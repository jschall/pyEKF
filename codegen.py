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
        exprs = loadExprsFromJSON(fn)['exprs']
        if len(exprs) > 1:
            for i in range(len(exprs)):
                filterOps['%s%s' % (n.upper(),i)] = exprs[i]
        else:
            filterOps['%s' % (n.upper(),)] = exprs[0]

    hashdefines = []
    hashdefines.extend(getConstants(filterOps))
    for k,v in filterOps.iteritems():
        hashdefines.extend(getSnippets(k,v))

    with open(cfile, 'w') as f:
        f.truncate()
        f.write(getHeader(hashdefines, 'EKF_'))
        print('Generated code saved to %s'%(cfile,))

def getConstants(filterOps):
    max_num_subx = 0
    for k,v in filterOps.iteritems():
        max_num_subx = max(len(v['output']['subx']),max_num_subx)

    x = filterOps['PREDICTION']['input']['x']
    u = filterOps['PREDICTION']['input']['u']
    ret = []
    ret.append(('NUM_STATES', x.rows))
    ret.append(('NUM_CONTROL_INPUTS', u.rows))
    ret.append(('MAX_NUM_SUBX', max_num_subx))
    for r in range(x.rows):
        ret.append(('X_IDX_'+str(x[r,0]).upper(), r))
    for r in range(u.rows):
        ret.append(('U_IDX_'+str(u[r,0]).upper(), r))
    return ret

def getSnippets(name, op):
    ret = []

    params = []
    substitutions = []
    for k in sorted(op['input'].keys()):
        if not isinstance(op['input'][k], MatrixBase):
            op['input'][k] = Matrix([[op['input'][k]]])

        nr,nc = op['input'][k].shape

        params.append('__'+k.upper())

        if (nr,nc) == (1,1):
            substitutions += zip(op['input'][k],Matrix(nr,nc,[Symbol(params[-1])]))
        elif nc == 1:
            substitutions += zip(op['input'][k],Matrix(nr,nc,symbols(params[-1]+'[0:%u]'%(nr,))))
        else:
            substitutions += zip(op['input'][k],Matrix(nr,nc,symbols(params[-1]+'[0:%u][0:%u]'%(nr,nc))))

    for k in sorted(op['output'].keys()):
        retparam = '__RET_'+k.upper()

        paramlist = params[:]
        paramlist.append(retparam)
        paramlist = ','.join(paramlist)
        macroproto = '%s_CALC_%s(%s)'%(name,k.upper(),paramlist)

        #ret.append((macroproto, ))
        print macroproto

        #if not isinstance(op['output'][k], MatrixBase):
            #op['output'][k] = Matrix([[op['output'][k]]])

        #nr,nc = op['output'][k].shape

        #op['output'][k] = op['output'][k].xreplace(dict(substitutions))

        #print op['output'][k]

    return ret

def getSnippetsCovPred(predictionjson):
    stateVector, Pn_O, subx = loadExprsFromJSON(predictionjson, ['stateVector', 'Pn_O', 'subx'])
    Pn_O = symmetricMatrixToList(Pn_O)

    ret = []

    snippet = ''
    for e in subx:
        snippet += ccode(e[1], assign_to=e[0])+' '
    snippet = '{ '+snippet+' }'

    ret.append(('COV_PRED_CALC_SUBX(_STATE, _COV, _SUBX)', snippet))

    snippet = ''
    for i in range(len(Pn_O)):
        snippet += ccode(Pn_O[i], assign_to='_NEXTCOV[%u]'%(i,))+' '
    snippet = '{ '+snippet+' }'

    ret.append(('COV_PRED_CALC_NEXTCOV(_STATE, _COV, _SUBX, _NEXTCOV)', snippet))

    return ret

def getFusionSnippets(json,name):
    name = name.upper()
    nObs, S_O, K_O, Pn_O, subx = loadExprsFromJSON(json, ['nObs', 'S_O', 'K_O', 'Pn_O', 'subx'])
    Pn_O = map(symmetricMatrixToList,Pn_O)

    ret = []

    for obs in range(nObs):
        obsname = '%s%u' % (name,obs) if nObs > 1 else name

        snippet = ''
        for e in subx[obs]:
            snippet += ccode(e[1], assign_to=e[0])+' '
        ret.append(('%s_CALC_SUBX(_STATE, _COV, _SUBX)' % (obsname,), '{ '+snippet+' }'))

        snippet = ccode(S_O[obs], assign_to='_S')+' '
        ret.append(('%s_CALC_S(_STATE, _COV, _SUBX, _S)' % (obsname,), '{ '+snippet+' }'))

        snippet = ''
        for i in range(len(K_O[obs])):
            snippet += ccode(K_O[obs][i], assign_to='_K[%u]'%(i,))+' '
        ret.append(('%s_CALC_K(_STATE, _COV, _SUBX, _K)' % (obsname,), '{ '+snippet+' }'))

        snippet = ''
        for i in range(len(Pn_O[obs])):
            snippet += ccode(Pn_O[obs][i], assign_to='_NEXTCOV[%u]'%(i,))+' '
        ret.append(('%s_CALC_NEXTCOV(_STATE, _COV, _SUBX, _NEXTCOV)' % (obsname,), '{ '+snippet+' }'))

    return ret

def wrapstring(string, linemax, delim):
    import re
    return re.sub("(.{1,%u})\s(?!$)" % (linemax,), '\g<1>'+delim, string)

def getHeader(hashdefines, prefix=''):
    ret = '/*\n'
    for (key,val) in hashdefines:
        if type(val) is str:
            ret += prefix+key+'\n'
    ret += '*/\n\n'
    for (key,val) in hashdefines:
        val = wrapstring(str(val),120,' \\\n')
        if '\n' in val:
            ret += '\n#define %s%s \\\n%s\n' % (prefix,key,val)
        else:
            ret += '#define %s%s %s\n' % (prefix,key,val)
    return ret
