from sympy import *
from helpers import *
from sys import exit
import math

# Parameters
dAngNoise = toVec(symbols('daxNoise dayNoise dazNoise'))
dVelNoise = toVec(symbols('dvxNoise dvyNoise dvzNoise'))
gravityNED = toVec(0.,0.,symbols('gravity'))
dt = Symbol('dt')
ptd = Symbol('ptd')
R_TAS = toVec(Symbol('R_TAS'))
R_BETA = toVec(Symbol('R_BETA'))
R_MAG = toVec(Symbol('R_MAG'))
R_LOS = toVec(Symbol('R_LOS'))
R_YAW = toVec(Symbol('R_YAW'))
R_DECL = toVec(Symbol('R_DECL'))

# Inputs
dAngMeas = toVec(symbols('dax day daz'))
dVelMeas = toVec(symbols('dvx dvy dvz'))

# States
estQuat = toVec(symbols('q0 q1 q2 q3'))
rotErr = toVec(symbols('rex rey rez'))
vn,ve,vd = symbols('vn ve vd')
velNED = toVec(vn,ve,vd)
posNED = toVec(symbols('pn pe pd'))
dAngBias = toVec(symbols('dax_b day_b daz_b'))
dAngScale = toVec(symbols('dax_s day_s daz_s'))
dVelBias = toVec(0.,0.,symbols('dvz_b'))
magBody = toVec(symbols('magx magy magz'))
magNED = toVec(symbols('magn mage magd'))
vwn, vwe = symbols('vwn vwe')
windNED = toVec(vwn,vwe,0.)
stateVector = toVec(rotErr,velNED,posNED,dAngBias,dAngScale,dVelBias[2],magNED,magBody,vwn,vwe)
nStates = len(stateVector)
stateVectorIndexed = toVec(symbols('_STATE[0:%u]'%(nStates,)))

# Covariance matrix
P = Matrix(nStates,nStates,symbols('_COV[0:%u][0:%u]' % (nStates,nStates)))
P_symmetric = P
for r in range(P_symmetric.rows):
    for c in range(P_symmetric.cols):
        if r > c:
            P_symmetric[c,r] = P_symmetric[r,c]

# Common computations
# Quaternion from body frame at time k to earth frame
truthQuat = quat_rotate_approx(estQuat,rotErr)
# Rotation matrix from body frame at time k to earth frame
Tbn = quat_to_matrix(truthQuat)

def deriveCovariancePrediction(jsonfile):
    print('Beginning covariance prediction derivation')
    # The prediction step predicts the state at time k+1 as a function of the
    # state at time k and the control inputs. This attitude estimation EKF is
    # formulated with the IMU data as control inputs rather than observations.

    # Rotation vector from body frame at time k to body frame at time k+1
    dAngTruth = dAngMeas.multiply_elementwise(dAngScale) - dAngBias

    # Change in velocity from time k to time k+1 in body frame at time k+1
    dVelTruth = dVelMeas - dVelBias

    truthQuatNew = quat_rotate_approx(truthQuat,dAngTruth)
    errQuatNew = quat_multiply(quat_inverse(estQuat),truthQuatNew)

    # States at time k+1
    #estQuatNew = quat_rotate_approx(estQuat, dAngTruth)
    rotErrNew = quat_to_rot_vec_approx(errQuatNew)
    velNEDNew = velNED+gravityNED*dt + Tbn*dVelTruth
    posNEDNew = posNED+velNED*dt
    dAngBiasNew = dAngBias
    dAngScaleNew = dAngScale
    dVelBiasNew = dVelBias
    magBodyNew = magBody
    magNEDNew = magNED
    vwnNew = vwn
    vweNew = vwe

    # f: state-transtition model
    f = toVec(rotErrNew,velNEDNew,posNEDNew,dAngBiasNew,dAngScaleNew,dVelBiasNew[2],magNEDNew,magBodyNew,vwnNew,vweNew)
    assert f.shape == stateVector.shape

    # F: linearized state-transition model
    F = f.jacobian(stateVector)

    # u: control input vector
    u = toVec(dAngMeas,dVelMeas)

    # G: control-influence matrix, AKA "B" in literature
    G = f.jacobian(u)

    # w_u_sigma: additive noise on u
    w_u_sigma = toVec(dAngNoise, dVelNoise)

    # Q_u: covariance of additive noise on u
    Q_u = diag(*w_u_sigma.multiply_elementwise(w_u_sigma))

    # Q: covariance of additive noise on x
    Q = G*Q_u*G.T

    # Pn: covariance matrix at time k+1
    Pn = F*P*F.T + Q

    # Optimizations
    Pn_O = Pn

    # zero the lower off-diagonals before extracting subexpressions
    Pn_O = zero_lower_offdiagonals(Pn_O)

    # assume that the P matrix is symmetrical
    Pn_O = Pn_O.subs(zip(P, P_symmetric)+zip(stateVector,stateVectorIndexed))

    Pn_O,subx = extractSubexpressions(Pn_O,'_SUBX')

    # make Pn_O symmetric
    Pn_O = copy_upper_to_lower_offdiagonals(Pn_O)

    saveExprsToJSON(jsonfile, {'stateVector':stateVector,'Pn_O':Pn_O,'subx':subx})
    print('Covariance predicton derivation saved to %s' % (jsonfile))

    ops = count_ops(zero_lower_offdiagonals(Pn_O))+count_ops(subx)
    print('%u subexpressions, %u ops' % (len(subx),ops))

def deriveAirspeedFusion(jsonfile):
    measPred = Matrix([[sqrt((vn-vwn)**2 + (ve-vwe)**2 + vd**2)]])

    deriveFusionSequential('Airspeed',jsonfile,measPred,R_TAS)

def deriveBetaFusion(jsonfile):
    Vbw = Tbn.T*(velNED-windNED)
    measPred = Matrix([[Vbw[1]/Vbw[0]]])

    deriveFusionSequential('Beta',jsonfile,measPred,R_BETA)

def deriveMagFusion(jsonfile):
    measPred = Tbn.T*magNED+magBody

    deriveFusionSequential('Mag',jsonfile,measPred,R_MAG)

def deriveOptFlowFusion(jsonfile):
    velBody = Tbn.T*velNED
    rangeToGround = (ptd-posNED[2])/Tbn[2,2]
    measPred = toVec(velBody[1]/rangeToGround, -velBody[0]/rangeToGround)

    deriveFusionSequential('Flow',jsonfile,measPred,R_LOS)

def deriveYaw321Fusion(jsonfile):
    measPred = toVec(atan(Tbn[1,0]/Tbn[0,0]))

    deriveFusionSequential('Yaw321',jsonfile,measPred,R_YAW)

def deriveYaw312Fusion(jsonfile):
    measPred = toVec(atan(-Tbn[0,1]/Tbn[1,1]))

    deriveFusionSequential('Yaw312',jsonfile,measPred,R_YAW)

def deriveDeclinationFusion(jsonfile):
    measPred = toVec(atan(magNED[1]/magNED[0]))

    deriveFusionSequential('Declination',jsonfile,measPred,R_YAW)

def deriveFusionSequential(fusionName,jsonfile,measPred,R):
    assert isinstance(measPred,MatrixBase) and isinstance(R,MatrixBase) and R.shape[0] == R.shape[1] == 1
    print('Beginning %s fusion derivation' % (fusionName,))

    H = measPred.jacobian(stateVector)

    nObs = measPred.rows
    S_O = [None for _ in range(nObs)]
    K_O = [None for _ in range(nObs)]
    Pn_O = [None for _ in range(nObs)]

    for i in range(nObs):
        S_O[i] = (H.row(i)*P*H.row(i).T + R)[0]
        K_O[i] = (P*H.row(i).T)/S_O[i]
        Pn_O[i] = (eye(nStates)-K_O[i]*H.row(i))*P
        Pn_O[i] = zero_lower_offdiagonals(Pn_O[i])

        subs = zip(P, P_symmetric)+zip(stateVector,stateVectorIndexed)
        S_O[i] = S_O[i].subs(subs)
        K_O[i] = K_O[i].subs(subs)
        Pn_O[i] = Pn_O[i].subs(subs)

    result = extractSubexpressions(S_O+K_O+Pn_O,'_SUBX')
    S_O = list(result[nObs*0:nObs*1])
    K_O = list(result[nObs*1:nObs*2])
    Pn_O = list(result[nObs*2:nObs*3])
    subx = result[-1]

    Pn_O = map(copy_upper_to_lower_offdiagonals, Pn_O)

    saveExprsToJSON(jsonfile, {'stateVector':stateVector,'nObs':nObs,'S_O':S_O,'K_O':K_O,'Pn_O':Pn_O,'subx':subx})
    print('%s fusion derivation saved to %s' % (fusionName,jsonfile))
    ops = count_ops(S_O)+count_ops(K_O)+count_ops(map(zero_lower_offdiagonals, Pn_O))+count_ops(subx)
    print('%u subexpressions, %u ops' % (len(subx),ops))
