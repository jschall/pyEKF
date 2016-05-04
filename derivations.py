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

# Inputs
dAngMeas = toVec(symbols('dax day daz'))
dVelMeas = toVec(symbols('dvx dvy dvz'))

# States
estQuat = toVec(symbols('q0 q1 q2 q3'))
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
stateVector = toVec(estQuat,velNED,posNED,dAngBias,dAngScale,dVelBias[2],magNED,magBody,vwn,vwe)
nStates = len(stateVector)
stateVectorIndexed = toVec(symbols('x[0:%u]' % (nStates)))

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P_symmetric = P
for r in range(P_symmetric.rows):
    for c in range(P_symmetric.cols):
        if r > c:
            P_symmetric[c,r] = P_symmetric[r,c]

def deriveCovariancePrediction(jsonfile):
    # The prediction step predicts the state at time k+1 as a function of the
    # state at time k and the control inputs. This attitude estimation EKF is
    # formulated with the IMU data as control inputs rather than observations.

    # Rotation vector from body frame at time k to body frame at time k+1
    dAngTruth = dAngMeas.multiply_elementwise(dAngScale) - dAngBias

    # Rotation matrix from body frame at time k to earth frame
    Tbn = quat_to_matrix(estQuat)

    # Change in velocity from time k to time k+1 in body frame at time k+1
    dVelTruth = dVelMeas - dVelBias

    # States at time k+1
    estQuatNew = quat_rotate_approx(estQuat, dAngTruth)
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
    f = toVec(estQuatNew,velNEDNew,posNEDNew,dAngBiasNew,dAngScaleNew,dVelBiasNew[2],magNEDNew,magBodyNew,vwnNew,vweNew)

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

    # assume that the P matrix is symmetrical
    Pn_O = Pn_O.subs(zip(P, P_symmetric))

    # zero the lower off-diagonals before extracting subexpressions
    Pn_O = zero_lower_offdiagonals(Pn_O)

    Pn_O,subx = extractSubexpressions(Pn_O,'subx')

    # make Pn_O symmetric
    Pn_O = copy_upper_to_lower_offdiagonals(Pn_O)

    saveExprsToJSON(jsonfile, {'Pn_O':Pn_O, 'subx':subx})

def deriveFusionSequential(jsonfile,measPred,R):
    assert isinstance(measPred,MatrixBase) and isinstance(R,MatrixBase) and R.shape[0] == R.shape[1] == 1

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

        subs = zip(P, P_symmetric)
        S_O[i] = S_O[i].subs(subs)
        K_O[i] = K_O[i].subs(subs)
        Pn_O[i] = Pn_O[i].subs(subs)


    result = extractSubexpressions(S_O+K_O+Pn_O,'subx')
    S_O = list(result[nObs*0:nObs*1])
    K_O = list(result[nObs*1:nObs*2])
    Pn_O = list(result[nObs*2:nObs*3])
    subx = result[-1]

    Pn_O = map(copy_upper_to_lower_offdiagonals, Pn_O)

    saveExprsToJSON(jsonfile, {'S_O':S_O,'K_O':K_O,'Pn_O':Pn_O,'subx':subx})

def deriveAirspeedFusion(jsonfile):
    measPred = Matrix([[sqrt((vn-vwn)**2 + (ve-vwe)**2 + vd**2)]])

    deriveFusionSequential(jsonfile,measPred,R_TAS)

def deriveBetaFusion(jsonfile):
    Tbn = quat_to_matrix(estQuat)
    Vbw = Tbn.T*(velNED-windNED)
    measPred = Matrix([[Vbw[1]/Vbw[0]]])

    deriveFusionSequential(jsonfile,measPred,R_BETA)

def deriveMagFusion(jsonfile):
    Tbn = quat_to_matrix(estQuat)
    measPred = Tbn.T*magNED+magBody

    deriveFusionSequential(jsonfile,measPred,R_MAG)

def deriveOptFlowFusion(jsonfile):
    Tbn = quat_to_matrix(estQuat)
    velBody = Tbn.T*velNED
    rangeToGroud = ((ptd-pd)/Tbn[2,2])
    measPred = toVec(velBody[1]/rangeToGround, -velBody[0]/rangeToGround)

    deriveFusionSequential(jsonfile,H,R_LOS)
