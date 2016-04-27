from sympy import *
from helpers import *
from sympy.printing.ccode import *
from sys import exit
import math

# Parameters
estQuat = toVec(symbols('q0 q1 q2 q3'))
dAngMeas = toVec(symbols('dax day daz'))
dVelMeas = toVec(symbols('dvx dvy dvz'))
dAngNoise = toVec(symbols('daxNoise dayNoise dazNoise'))
dVelNoise = toVec(symbols('dvxNoise dvyNoise dvzNoise'))
gravityNED = toVec(0.,0.,symbols('gravity'))
dt = Symbol('dt')
R_TAS, R_BETA, R_MAG = symbols('R_TAS R_BETA R_MAG')

# States
rotErr = toVec(symbols('rotErrX rotErrY rotErrZ'))
vn,ve,vd = symbols('vn ve vd')
velNED = toVec(vn,ve,vd)
posNED = toVec(symbols('pn pe pd'))
dAngBias = toVec(symbols('dax_b day_b daz_b'))
dAngScale = toVec(symbols('dax_s day_s daz_s'))
dVelBias = toVec(0,0,symbols('dvz_b'))
magBody = toVec(symbols('magX magY magZ'))
magNED = toVec(symbols('magN magE magD'))
vwn, vwe = symbols('vwn vwe')
windNED = toVec(vwn,vwe,0.)
stateVector = toVec(rotErr,velNED,posNED,dAngBias,dAngScale,dVelBias[2],magNED,magBody,vwn,vwe)
nStates = len(stateVector)

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P_symmetric = P
for r in range(P_symmetric.rows):
    for c in range(P_symmetric.cols):
        if r > c:
            P_symmetric[c,r] = P_symmetric[r,c]

# Common computations
errQuat = rot_vec_to_quat_approx(rotErr)
truthQuat = quat_multiply(estQuat, errQuat)
Tbn = quat_to_matrix(truthQuat)

def deriveCovariancePrediction(jsonfile):
    # Prediction step
    dAngTruth = dAngMeas.multiply_elementwise(dAngScale) - dAngBias
    deltaQuat = rot_vec_to_quat_approx(dAngTruth)
    truthQuatNew = quat_multiply(truthQuat, deltaQuat)
    estQuatNew = quat_multiply(estQuat, deltaQuat)
    errQuatNew = quat_multiply(quat_inverse(estQuatNew), truthQuatNew)
    rotErrNew = quat_to_rot_vec_approx(errQuatNew)
    dVelTruth = dVelMeas - dVelBias
    velNEDNew = velNED+gravityNED*dt + Tbn*dVelTruth
    posNEDNew = posNED+velNED*dt
    dAngBiasNew = dAngBias
    dAngScaleNew = dAngScale
    dVelBiasNew = dVelBias
    magBodyNew = magBody
    magNEDNew = magNED
    vwnNew = vwn
    vweNew = vwe

    stateVectorNew = toVec(rotErrNew,velNEDNew,posNEDNew,dAngBiasNew,dAngScaleNew,dVelBiasNew[2],magNEDNew,magBodyNew,vwnNew,vweNew)

    F = stateVectorNew.jacobian(stateVector)
    F = F.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])

    G = -stateVectorNew.jacobian(toVec(dAngMeas,dVelMeas))
    G = G.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])

    distVector = toVec(dAngNoise, dVelNoise)
    distMatrix = diag(*distVector.multiply_elementwise(distVector))
    Q = G*distMatrix*G.T

    PP = F*P*F.T + Q

    # Optimizations
    PP_O = PP

    # assume that the P matrix is symmetrical
    PP_O = PP_O.subs(zip(P, P_symmetric))

    # zero the lower off-diagonals
    for r in range(PP_O.rows):
        for c in range(PP_O.cols):
            if r > c:
                PP_O[r,c] = 0.

    PP_O,PP_S = extractSubexpressions(PP_O,'PP_S')

    saveExprsToJSON(jsonfile, {'PP':PP, 'PP_O':PP_O, 'PP_S':PP_S})

def deriveAirspeedFusion(jsonfile):
    VtasPred = Matrix([[sqrt((vn-vwn)**2 + (ve-vwe)**2 + vd**2)]])

    H = VtasPred.jacobian(stateVector)

    H_sym = Matrix(1,nStates,symbols('H[0:%u][0:%u]' % (1,nStates)))

    H = H.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])
    H_O, H_S = extractSubexpressions(H,'H_S')

    K = (P*H_sym.T)/(H_sym*P*H_sym.T + Matrix([[R_TAS]]))[0]

    K_O = K.subs(zip(H_sym, H_O)+zip(P,P_symmetric))
    K = K.subs(zip(H_sym,H))

    K_O, K_S = extractSubexpressions(K_O,'K_S')

    saveExprsToJSON(jsonfile, {'H':H, 'H_O':H_O, 'H_S':H_S, 'K':K, 'K_O':K_O, 'K_S':K_S})

def deriveBetaFusion(jsonfile):
    Vbw = Tbn.T*(velNED-windNED)
    betaPred = Matrix([[Vbw[1]/Vbw[0]]])

    H = betaPred.jacobian(stateVector)

    H_sym = Matrix(1,nStates,symbols('H[0:%u][0:%u]' % (1,nStates)))

    H = H.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])
    H_O, H_S = extractSubexpressions(H,'H_S')

    K = (P*H_sym.T)/(H_sym*P*H_sym.T + Matrix([[R_BETA]]))[0]
    K_O = K.subs(zip(H_sym, H_O)+zip(P,P_symmetric))
    K = K.subs(zip(H_sym,H))

    K_O, K_S = extractSubexpressions(K_O,'K_S')

    saveExprsToJSON(jsonfile, {'H':H, 'H_O':H_O, 'H_S':H_S, 'K':K, 'K_O':K_O, 'K_S':K_S})

def deriveMagFusion(jsonfile):
    magPred = Tbn.T*magNED+magBody
    H = magPred.jacobian(stateVector)

    H_sym = Matrix(1,nStates,symbols('H[0:%u][0:%u]' % (1,nStates)))

    H = H.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])

    H_O, H_S = extractSubexpressions(H,'H_S')

    K = (P*H_sym.T)/(H_sym*P*H_sym.T + Matrix([[R_MAG]]))[0]

    KX = K.subs(zip(H_sym,H.row(0)))
    KY = K.subs(zip(H_sym,H.row(1)))
    KZ = K.subs(zip(H_sym,H.row(2)))
    KX_O = K.subs(zip(H_sym,H_O.row(0))+zip(P,P_symmetric))
    KY_O = K.subs(zip(H_sym,H_O.row(1))+zip(P,P_symmetric))
    KZ_O = K.subs(zip(H_sym,H_O.row(2))+zip(P,P_symmetric))

    KX_O, KY_O, KZ_O, K_S = extractSubexpressions([KX_O,KY_O,KZ_O],'K_S')

    saveExprsToJSON(jsonfile, {'H':H, 'H_O':H_O, 'H_S':H_S, 'KX_O':KX, 'KY_O':KY, 'KZ_O':KZ, 'K_S':K_S})
