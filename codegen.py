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
    for k,v in jsondict.iteritems():
        exprs = loadExprsFromJSON(fn)['exprs']
        if len(exprs) > 1:
            for i in range(len(exprs)):
                filterOps['%s%s' % (k.upper(),i)] = exprs[i]
        else:
            filterOps['%s' % (k.upper(),)] = exprs[0]

    hashdefines = []
    hashdefines.extend(getConstants(jsondict))
    for k,v in filterOps.iteritems():
        hashdefines.extend(getSnippets(k,v))

    with open(cfile, 'w') as f:
        f.truncate()
        f.write(getHeader(hashdefines, 'EKF_'))
        print('Generated code saved to %s'%(cfile,))

def getConstants(filterOps):
    max_num_subx = 0
    for k,v in filterOps:
        max_num_subx = max(len(v['output']['subx']),max_num_subx)

    x = filterOps['prediction']['input']['x']
    constants = []
    constants.append(('NUM_STATES', x.rows))
    constants.append(('MAX_NUM_SUBX', max_num_subx))
    for r in range(x.rows):
        constants.append(('STATEIDX_'+str(x[r,0]).upper(), r))
    return constants

def getSnippets(name, op):
    pass

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
