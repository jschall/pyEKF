from sympy import *
from helpers import *
from sys import exit
import math

# Parameters
dAngNoise = toVec(symbols('daxNoise dayNoise dazNoise'))
dVelNoise = toVec(symbols('dvxNoise dvyNoise dvzNoise'))
gravity = Symbol('gravity')
gravityNED = toVec(0.,0.,gravity)
dt = Symbol('dt')
ptd = Symbol('ptd')

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

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P_symmetric = copy_upper_to_lower_offdiagonals(P)

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

    # zero the lower off-diagonals before extracting subexpressions
    Pn = zero_lower_offdiagonals(Pn)

    # assume that the P matrix is symmetrical
    Pn = Pn.xreplace(dict(zip(P, P_symmetric)))

    Pn,subx = extractSubexpressions(Pn,'_SUBX')

    # make Pn symmetric
    Pn = copy_upper_to_lower_offdiagonals(Pn)

    funcParams = {'quat':estQuat,'x':stateVector,'P':P_symmetric,'u':u,'w_u_sigma':w_u_sigma,'gravity':gravity,'dt':dt}

    funcs = {}

    funcs['subx'] = {}
    funcs['subx']['params'] = funcParams
    funcs['subx']['ret'] = toVec([x[1] for x in subx])
    funcs['subx']['retsymbols'] = toVec([x[0] for x in subx])

    funcParams = funcParams.copy()
    funcParams['subx'] = toVec([x[0] for x in subx])

    funcs['P'] = {}
    funcs['P']['params'] = funcParams
    funcs['P']['ret'] = Pn

    check_funcs(funcs)

    saveExprsToJSON(jsonfile, {'operations':[funcs]})

    print('Covariance predicton: derivation saved to %s.' % (jsonfile,))

def derivePosFusion(jsonfile):
    measPred = posNED

    deriveFusionSequential('Pos',jsonfile,measPred)

def deriveVelFusion(jsonfile):
    measPred = velNED

    deriveFusionSequential('Vel',jsonfile,measPred)

def deriveAirspeedFusion(jsonfile):
    measPred = toVec(sqrt((vn-vwn)**2 + (ve-vwe)**2 + vd**2))

    deriveFusionSequential('Airspeed',jsonfile,measPred)

def deriveBetaFusion(jsonfile):
    Vbw = Tbn.T*(velNED-windNED)
    measPred = toVec(Vbw[1]/Vbw[0])

    deriveFusionSequential('Beta',jsonfile,measPred)

def deriveMagFusion(jsonfile):
    measPred = Tbn.T*magNED+magBody

    deriveFusionSequential('Mag',jsonfile,measPred)

def deriveOptFlowFusion(jsonfile):
    velBody = Tbn.T*velNED
    rangeToGround = (ptd-posNED[2])/Tbn[2,2]
    measPred = toVec(velBody[1]/rangeToGround, -velBody[0]/rangeToGround)

    deriveFusionSequential('Flow',jsonfile,measPred,{'ptd':ptd})

def deriveYaw321Fusion(jsonfile):
    measPred = toVec(atan(Tbn[1,0]/Tbn[0,0]))

    deriveFusionSequential('Yaw321',jsonfile,measPred)

def deriveYaw312Fusion(jsonfile):
    measPred = toVec(atan(-Tbn[0,1]/Tbn[1,1]))

    deriveFusionSequential('Yaw312',jsonfile,measPred)

def deriveDeclinationFusion(jsonfile):
    measPred = toVec(atan(magNED[1]/magNED[0]))

    deriveFusionSequential('Declination',jsonfile,measPred)

def deriveFusionSimultaneous(fusionName,jsonfile,measPred,additionalinputs={}):
    assert isinstance(measPred,MatrixBase) and measPred.cols == 1
    print('Beginning %s fusion derivation' % (fusionName,))

    nObs = measPred.rows
    substitute = dict(zip(P, P_symmetric))

    R = diag(*symbols('R[0:%u]' % (nObs,)))

    H = measPred.jacobian(stateVector)
    I = eye(nStates)

    funcs = {}

    S = (H*P*H.T + R)
    K = P*H.T*quickinv_sym(S)
    K = K.xreplace(substitute)
    K, K_subx = extractSubexpressions([K],'subx')

    Pn = (I-K*H)*P*(I-K*H).T+K*R*K.T # Joseph form for numerical reasons if desired

    Pn = zero_lower_offdiagonals(Pn)
    Pn = Pn.xreplace(substitute)

    Pn, P_subx = extractSubexpressions([Pn],'subx')

    ops = [count_ops(K), count_ops(K_subx), count_ops(Pn), count_ops(P_subx)]



def deriveFusionSequential(fusionName,jsonfile,measPred,additionalinputs={}):
    assert isinstance(measPred,MatrixBase) and measPred.cols == 1
    print('Beginning %s fusion derivation' % (fusionName,))
    R = toVec(Symbol('R'))

    nObs = measPred.rows

    H_simul = measPred.jacobian(stateVector)
    I = eye(nStates)

    operations = []

    for i in range(nObs):
        funcs = {}

        H = H_simul.row(i)

        S = (H*P*H.T + R)[0]
        K = (P*H.T)/S
        Pn = (I-K*H)*P #*(I-K*H).T+K*R*K.T # Joseph form for numerical reasons if desired

        substitute = dict(zip(P, P_symmetric))
        S = S.xreplace(substitute)
        K = K.xreplace(substitute)
        Pn = zero_lower_offdiagonals(Pn)
        Pn = Pn.xreplace(substitute)

        S, K, Pn, subx = extractSubexpressions([S,K,Pn],'subx')

        ops = [count_ops(S), count_ops(K), count_ops(Pn), count_ops(subx)]
        print ops, sum(ops), len(subx)

        Pn = copy_upper_to_lower_offdiagonals(Pn)

        funcParams = {'quat':estQuat,'x':stateVector,'P':P_symmetric,'R':R}
        funcParams.update(additionalinputs)


        funcs['subx'] = {}
        funcs['subx']['params'] = funcParams
        funcs['subx']['ret'] = toVec([x[1] for x in subx])
        funcs['subx']['retsymbols'] = toVec([x[0] for x in subx])

        funcParams = funcParams.copy()
        funcParams['subx'] = toVec([x[0] for x in subx])

        funcs['S'] = {}
        funcs['S']['params'] = funcParams
        funcs['S']['ret'] = S

        funcs['K'] = {}
        funcs['K']['params'] = funcParams
        funcs['K']['ret'] = K

        funcs['P'] = {}
        funcs['P']['params'] = funcParams
        funcs['P']['ret'] = Pn
        operations.append(funcs)

    for funcs in operations:
        check_funcs(funcs)

    saveExprsToJSON(jsonfile, {'operations':operations})

    print('%s fusion: derivation saved to %s.' % (fusionName,jsonfile))
