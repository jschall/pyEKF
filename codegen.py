from sympy import *
from sympy.printing.ccode import *
from helpers import *
import math

# - max subexpression array size
# - state array size
# - state index enum
# - covariance prediction
# - mag x,y,z innovation variance, kalman gain or state update, covariance update
# - vel x,y,z innovation variance, kalman gain or state update, covariance update
# - pos x,y,z innovation variance, kalman gain or state update, covariance update
# - airspeed innovation variance, kalman gain or state update, covariance update
# - beta innovation variance, kalman gain or state update, covariance update
# - yaw312 innovation variance, kalman gain or state update, covariance update
# - yaw321 innovation variance, kalman gain or state update, covariance update
# - declination innovation variance, kalman gain or state update, covariance update

def generateCode(jsondict, cfile):
    # extract filter operations
    filterOps = {}
    for n,fn in jsondict.iteritems():
        operations = loadExprsFromJSON(fn)['operations']
        if len(operations) > 1:
            for i in range(len(operations)):
                filterOps['%s%s' % (n.upper(),i)] = operations[i]
        else:
            filterOps['%s' % (n.upper(),)] = operations[0]

    hashdefines = []
    hashdefines.extend(getConstants(filterOps))
    for k in sorted(filterOps.keys()):
        hashdefines.extend(getSnippetDefines(k,filterOps[k]))

    with open(cfile, 'w') as f:
        f.truncate()
        f.write(getHeader(hashdefines, 'EKF_'))
        print('Generated code saved to %s'%(cfile,))

def getConstants(filterOps):
    max_num_subx = 0
    for k,v in filterOps.iteritems():
        max_num_subx = max(len(v['subx']['ret']),max_num_subx)

    x = filterOps['COVPRED']['P']['params']['x']
    u = filterOps['COVPRED']['P']['params']['u']
    ret = []
    ret.append(('NUM_STATES', x.rows))
    ret.append(('NUM_CONTROL_INPUTS', u.rows))
    ret.append(('MAX_NUM_SUBX', max_num_subx))
    for r in range(x.rows):
        ret.append(('STATE_IDX_'+str(x[r,0]).upper(), r))
    for r in range(u.rows):
        ret.append(('U_IDX_'+str(u[r,0]).upper(), r))
    return ret

def getSnippetDefines(opname, funcs):
    ret = []

    for retName,func in funcs.iteritems():
        retParamName = '__RET_'+retName.upper()

        paramlist = []
        substitutions = []
        for paramname in sorted(func['params'].keys()):
            if not isinstance(func['params'][paramname], MatrixBase):
                func['params'][paramname] = Matrix([[func['params'][paramname]]])

            if func['params'][paramname].shape != (1,1) and func['params'][paramname].is_symmetric(simplify=False):
                func['params'][paramname] = symmetricMatrixToVec(func['params'][paramname])

            nr,nc = func['params'][paramname].shape

            paramlist.append('__'+paramname.upper())

            if (nr,nc) == (1,1):
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, [Symbol(paramlist[-1])]))
            elif nc == 1:
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, symbols(paramlist[-1]+'[0:%u]'%(nr,))))
            else:
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, symbols(paramlist[-1]+'[0:%u][0:%u]'%(nr,nc))))

        if 'retsymbols' in func:
            if not isinstance(func['retsymbols'], MatrixBase):
                func['retsymbols'] = Matrix([[func['retsymbols']]])

            nr,nc = func['retsymbols'].shape

            if (nr,nc) == (1,1):
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, [Symbol(retParamName)]))
            elif nc == 1:
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, symbols(retParamName+'[0:%u]'%(nr,))))
            else:
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, symbols(retParamName+'[0:%u][0:%u]'%(nr,nc))))

        func['ret'] = func['ret'].xreplace(dict(substitutions))

        defineName = '%s_CALC_%s(%s)' % (opname,retName.upper(),','.join(paramlist+[retParamName]))
        ret.append((defineName,getSnippet(retParamName,func['ret'])))

    return ret

def getSnippet(retParamName, outputMatrix):
    if not isinstance(outputMatrix, MatrixBase):
        outputMatrix = Matrix([[outputMatrix]])

    if outputMatrix.is_symmetric(simplify=False):
        outputMatrix = symmetricMatrixToVec(outputMatrix)

    nr,nc = outputMatrix.shape

    if (nr,nc) == (1,1):
        retMatrix = Matrix(nr,nc, [Symbol(retParamName)])
    elif nc == 1:
        retMatrix = Matrix(nr,nc, symbols(retParamName+'[0:%u]'%(nr,)))
    else:
        retMatrix = Matrix(nr,nc, symbols(retParamName+'[0:%u][0:%u]'%(nr,nc)))

    ret = ''
    for assignee,expr in zip(retMatrix,outputMatrix):
        ret += ccode(expr, assign_to=assignee)+' '

    return str(ret)

def wrapstring(string, linemax, delim):
    strings = string.split(' ')
    lines = [strings[0]]
    for s in strings[1:]:
        if len(lines[-1]) != 0 and len(lines[-1])+len(s) > linemax:
            lines.append(s)
        else:
            lines[-1] += ' '+s
    return delim.join(lines)

def getHeader(hashdefines, prefix=''):
    ret = '/*\n'
    for (key,val) in hashdefines:
        if type(val) is str:
            ret += prefix+key+'\n'
    ret += '*/\n\n'
    for (key,val) in hashdefines:
        val = wrapstring(str(val),100,' \\\n')
        if '\n' in val:
            ret += '\n#define %s%s \\\n%s\n' % (prefix,key,val)
        else:
            ret += '#define %s%s %s\n' % (prefix,key,val)
    return ret
