from sympy import *
from sympy.printing.lambdarepr import lambdarepr
from helpers import *
import math

def generateCode(jsondict, pyfile):
    # extract filter operations
    filterOps = {}
    for n,fn in jsondict.iteritems():
        filterOps[n.upper()] = deserialize_exprs_in_structure(loadExprsFromJSON(fn)['funcs'])

    constants = getConstants(filterOps)
    functions = []
    for k in sorted(filterOps.keys()):
        functions.extend(getFunctions(k,filterOps[k]))

    fileCont = 'import numpy as np\nfrom math import *\n\n'
    for name, val in constants:
        fileCont += 'EKF_%s = %s\n' % (name,val)
    fileCont += '\n'
    for func in functions:
        fileCont += func+'\n'

    with open(pyfile, 'wb') as f:
        f.truncate()
        f.write(fileCont)
        print('Generated code saved to %s'%(pyfile,))

def getConstants(filterOps):
    max_num_subx = 0
    for k,v in filterOps.iteritems():
        if 'subx' in v:
            max_num_subx = max(len(v['subx']['ret']),max_num_subx)

    x = filterOps['PREDICTION']['cov']['params']['x']
    u = filterOps['PREDICTION']['cov']['params']['u']
    ret = []
    ret.append(('NUM_STATES', x.rows))
    ret.append(('NUM_CONTROL_INPUTS', u.rows))
    ret.append(('MAX_NUM_SUBX', max_num_subx))
    for r in range(x.rows):
        ret.append(('STATE_IDX_'+str(x[r,0]).upper(), r))
    for r in range(u.rows):
        ret.append(('U_IDX_'+str(u[r,0]).upper(), r))
    return ret

def getFunctions(opname, funcs):
    ret = []

    for retName,func in funcs.iteritems():
        retParamName = 'ret_'+retName

        paramlist = []
        substitutions = []
        for paramname in sorted(func['params'].keys()):
            if not isinstance(func['params'][paramname], MatrixBase):
                func['params'][paramname] = Matrix([[func['params'][paramname]]])

            nr,nc = func['params'][paramname].shape

            paramlist.append(paramname)

            if (nr,nc) == (1,1):
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, [Symbol(paramlist[-1])]))
            elif nc == 1:
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, symbols(paramlist[-1]+'[0:%u]'%(nr,))))
            else:
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, symbols(paramlist[-1]+'[0:%u\\,0:%u]'%(nr,nc))))

        if 'retsymbols' in func:
            if not isinstance(func['retsymbols'], MatrixBase):
                func['retsymbols'] = Matrix([[func['retsymbols']]])

            nr,nc = func['retsymbols'].shape

            if (nr,nc) == (1,1):
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, [Symbol(retParamName)]))
            elif nc == 1:
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, symbols(retParamName+'[0:%u]'%(nr,))))
            else:
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, symbols(retParamName+'[0:%u\\,0:%u]'%(nr,nc))))

        nr,nc = func['ret'].shape
        if (nr,nc) == (1,1):
            retMatrix = Matrix(nr,nc, [Symbol(retParamName)])
        elif nc == 1:
            retMatrix = Matrix(nr,nc, symbols(retParamName+'[0:%u]'%(nr,)))
        else:
            retMatrix = Matrix(nr,nc, symbols(retParamName+'[0:%u\\,0:%u]'%(nr,nc)))

        func['ret'] = func['ret'].xreplace(dict(substitutions))

        funcText = 'def EKF_%s_CALC_%s(%s):\n' % (opname,retName.upper(),', '.join(paramlist))

        if (nr,nc) == (1,1):
            funcText += '    '+retParamName+' = 0\n'
        elif nc == 1:
            funcText += '    '+retParamName+' = np.zeros(%u)\n' % func['ret'].shape[0]
        else:
            funcText += '    '+retParamName+' = np.zeros((%u,%u))\n' % func['ret'].shape

        for assignee,expr in zip(retMatrix,func['ret']):
            funcText += '    '+str(assignee)+' = '+lambdarepr(expr)+'\n'

        funcText += '\n    return %s\n' % retParamName

        ret.append(funcText)

    return ret

def getSnippet(retParamName, outputMatrix):
    if not isinstance(outputMatrix, MatrixBase):
        outputMatrix = Matrix([[outputMatrix]])

    outputMatrix = outputMatrix



    ret = ''
    for assignee,expr in zip(retMatrix,outputMatrix):
        ret += '    '+str(assignee)+' = '+lambdarepr(expr)+'\n'

    return str(ret)
