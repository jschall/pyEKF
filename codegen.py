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

def generateCode(predictionjson, airspeedjson, cfile):
    assert checkFileVersionMatch([predictionjson, airspeedjson])

    hashdefines = []
    hashdefines.extend(getConstants(predictionjson))
    hashdefines.extend(getSnippetsCovPred(predictionjson))
    hashdefines.extend(getFusionSnippets(airspeedjson,'AIRSPEED'))

    with open(cfile, 'w') as f:
        f.truncate()
        f.write(getHeader(hashdefines))
        print('Generated code saved to %s'%(cfile,))

def getHeader(hashdefines):
    ret = ''
    for (key,val) in hashdefines:
        ret += '#define %s %s\n\n' % (key,val)
    return ret

def checkFileVersionMatch(fnames):
    stateVectors = map(lambda x: ImmutableDenseMatrix(loadExprsFromJSON(x, ['stateVector'])[0]),fnames)
    return len(set(stateVectors)) <= 1

def getConstants(predictionjson):
    stateVector,subx = loadExprsFromJSON(predictionjson, ['stateVector','subx'])
    constants = []
    constants.append(('NUM_STATES', stateVector.rows))
    constants.append(('MAX_NUM_SUBX', len(subx)))
    for r in range(stateVector.rows):
        constants.append(('STATEIDX_'+str(stateVector[r,0]).upper(), r))
    return constants

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
    nObs, S_O, K_O, Pn_O, subx = loadExprsFromJSON(json, ['nObs', 'S_O', 'K_O', 'NEXTP', 'subx'])
    Pn_O = map(symmetricMatrixToList,Pn_O)

    ret = []

    snippet = ''
    for e in subx:
        snippet += ccode(e[1], assign_to=e[0])+' '
    snippet = '{ '+snippet+' }'

    ret.append(('%s_CALC_SUBX(_STATE, _COV, _SUBX)' % (name,), snippet))

    for obs in range(nObs):
        obsname = '%s%u' % (name,obs) if nObs > 1 else name
        snippet = ccode(S_O[obs], assign_to='_S')+' '
        snippet = '{ '+snippet+' }'

        ret.append(('%s_CALC_S(_STATE, _COV, _SUBX, _S)' % (obsname,), snippet))

        snippet = ''
        for i in range(len(K_O[obs])):
            snippet += ccode(K_O[obs][i], assign_to='_K[%u]'%(i,))+' '
        snippet = '{ '+snippet+' }'

        ret.append(('%s_CALC_K(_STATE, _COV, _SUBX, _K)' % (obsname,), snippet))

        snippet = ''
        for i in range(len(Pn_O[obs])):
            snippet += ccode(Pn_O[obs][i], assign_to='_NEXTCOV[%u]'%(i,))+' '
        snippet = '{ '+snippet+' }'

        ret.append(('%s_CALC_NEXTCOV(_STATE, _COV, _SUBX, _NEXTCOV)' % (obsname,), snippet))

    return ret
