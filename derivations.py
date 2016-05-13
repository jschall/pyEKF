from sympy import *
from helpers import *
from sys import exit
import math

# Parameters
dAngNoise = toVec(symbols('daxNoise dayNoise dazNoise'))
dVelNoise = toVec(symbols('dvxNoise dvyNoise dvzNoise'))
gravity = Symbol('gravity')
gravityNED = toVec(0,0,gravity)
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
dVelBias = toVec(0,0,symbols('dvz_b'))
magBody = toVec(symbols('magx magy magz'))
magNED = toVec(symbols('magn mage magd'))
vwn, vwe = symbols('vwn vwe')
windNED = toVec(vwn,vwe,0)
stateVector = toVec(rotErr,velNED,posNED,dAngBias,dAngScale,dVelBias[2],magNED,magBody,vwn,vwe)
nStates = len(stateVector)

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P = copy_upper_to_lower_offdiagonals(P)

# Common computations
# Quaternion from body frame at time k to earth frame
truthQuat = quat_rotate_approx(estQuat,rotErr)
# Rotation matrix from body frame at time k to earth frame
Tbn = quat_to_matrix(truthQuat)

def deriveCovariancePrediction(jsonfile):
    print('Beginning covariance prediction derivation')
    t1 = datetime.datetime.now()

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

    # P_n: covariance matrix at time k+1
    P_n = F*P*F.T + Q

    # Optimizations
    P_n = upperTriangularToVec(P_n)
    P_n,subx = extractSubexpressions(P_n,'subx',threshold=10)

    # Output generation
    funcParams = {'quat':estQuat,'x':stateVector,'P':upperTriangularToVec(P),'u':u,'w_u_sigma':w_u_sigma,'gravity':gravity,'dt':dt}

    funcs = {}

    funcs['subx'] = {}
    funcs['subx']['params'] = funcParams
    funcs['subx']['ret'] = toVec([x[1] for x in subx])
    funcs['subx']['retsymbols'] = toVec([x[0] for x in subx])

    funcParams = funcParams.copy()
    funcParams['subx'] = toVec([x[0] for x in subx])

    funcs['cov'] = {}
    funcs['cov']['params'] = funcParams
    funcs['cov']['ret'] = P_n

    check_funcs(funcs)

    saveExprsToJSON(jsonfile, {'funcs':funcs})

    op_count, subx_count = getOpStats(funcs)
    t2 = datetime.datetime.now()
    print('%s Covariance predicton: derivation saved to %s. %u ops, %u subexpressions.' % (t2-t1,jsonfile,op_count,subx_count))

def derivePosNEFusion(jsonfile):
    measPred = posNED[0:2,:]

    deriveFusionSimultaneous('posNE',jsonfile,measPred)

def derivePosDFusion(jsonfile):
    measPred = posNED[2:3,:]

    deriveFusionSimultaneous('posD',jsonfile,measPred)

def deriveVelNEFusion(jsonfile):
    measPred = velNED[0:2,:]

    deriveFusionSimultaneous('velNE',jsonfile,measPred)

def deriveVelDFusion(jsonfile):
    measPred = velNED[2:3,:]

    deriveFusionSimultaneous('velD',jsonfile,measPred)

def deriveAirspeedFusion(jsonfile):
    measPred = toVec(sqrt((vn-vwn)**2 + (ve-vwe)**2 + vd**2))

    deriveFusionSimultaneous('airspeed',jsonfile,measPred)

def deriveBetaFusion(jsonfile):
    Vbw = Tbn.T*(velNED-windNED)
    measPred = toVec(Vbw[1]/Vbw[0])

    deriveFusionSimultaneous('beta',jsonfile,measPred,{'quat':estQuat})

def deriveMagFusion(jsonfile):
    measPred = Tbn.T*magNED+magBody

    deriveFusionSimultaneous('mag',jsonfile,measPred,{'quat':estQuat})

def deriveOptFlowFusion(jsonfile):
    velBody = Tbn.T*velNED
    rangeToGround = (ptd-posNED[2])/Tbn[2,2]
    measPred = toVec(velBody[1]/rangeToGround, -velBody[0]/rangeToGround)

    deriveFusionSimultaneous('flow',jsonfile,measPred,{'ptd':ptd,'quat':estQuat})

def deriveYaw321Fusion(jsonfile):
    measPred = toVec(atan(Tbn[1,0]/Tbn[0,0]))

    deriveFusionSimultaneous('yaw321',jsonfile,measPred,{'quat':estQuat})

def deriveYaw312Fusion(jsonfile):
    measPred = toVec(atan(-Tbn[0,1]/Tbn[1,1]))

    deriveFusionSimultaneous('yaw312',jsonfile,measPred,{'quat':estQuat})

def deriveDeclinationFusion(jsonfile):
    measPred = toVec(atan(magNED[1]/magNED[0]))

    deriveFusionSimultaneous('declination',jsonfile,measPred)

def deriveFusionSimultaneous(fusionName,jsonfile,measPred,additionalinputs={},R_type='scalar'):
    assert isinstance(measPred,MatrixBase) and measPred.cols == 1
    print('Beginning %s fusion derivation' % (fusionName,))
    t1 = datetime.datetime.now()

    nObs = measPred.rows
    I = eye(nStates)

    # Define symbols
    z = toVec(symbols('z[0:%u]' % (nObs,))) # Measurement
    if R_type == 'matrix':
        R_param = Matrix(nObs,nObs, symbols('R[0:%u][0:%u]' % (nObs,nObs)))
        R = R_param
    elif R_type == 'vector':
        R_param = toVec(symbols('R[0:%u]' % (nObs,)))
        R = diag(*R_param)
    elif R_type == 'scalar':
        R_param = Symbol('R')
        R = eye(nObs)*R_param

    # Intermediates
    y = z-measPred                       # Innovation
    H = measPred.jacobian(stateVector)   # Obervation sensitivity matrix
    S = H*P*H.T + R                      # Innovation covariance
    S_I = quickinv_sym(S)                # Innovation covariance inverse
    K = P*H.T*S_I                        # Near-optimal Kalman gain

    y,H,S_I,K,temp_subx = extractSubexpressions([y,H,S_I,K],'temp')

    # Outputs

    # NOTE: The covariance update involves subtraction and can result in loss
    # of symmetry and positive definiteness due to rounding errors. Joseph's
    # form covariance update avoids this at the expense of computation burden.

    NIS = y.T*S_I*y                      # Normalized innovation squared
    x_n = stateVector+K*y                # Updated state vector
    P_n = (I-K*H)*P*(I-K*H).T+K*R*K.T    # Updated covariance matrix

    # Optimizations
    P_n = upperTriangularToVec(P_n)
    y, NIS, x_n, P_n, subx = extractSubexpressions([y,NIS,x_n,P_n],'subx',threshold=10,prev_subx=temp_subx)

    funcParams = {'x':stateVector,'P':upperTriangularToVec(P),'R':R_param,'z':z}
    funcParams.update(additionalinputs)

    funcs = {}
    funcs['subx'] = {}
    funcs['subx']['params'] = funcParams
    funcs['subx']['ret'] = toVec([x[1] for x in subx])
    funcs['subx']['retsymbols'] = toVec([x[0] for x in subx])

    funcParams = funcParams.copy()
    funcParams['subx'] = toVec([x[0] for x in subx])

    funcs['innov'] = {}
    funcs['innov']['params'] = funcParams
    funcs['innov']['ret'] = y

    funcs['NIS'] = {}
    funcs['NIS']['params'] = funcParams
    funcs['NIS']['ret'] = NIS

    funcs['state'] = {}
    funcs['state']['params'] = funcParams
    funcs['state']['ret'] = x_n

    funcs['cov'] = {}
    funcs['cov']['params'] = funcParams
    funcs['cov']['ret'] = P_n

    saveExprsToJSON(jsonfile, {'funcs':funcs})

    op_count, subx_count = getOpStats(funcs)
    t2 = datetime.datetime.now()
    print('%s %s fusion: derivation saved to %s. %u ops, %u subexpressions.' % (t2-t1,fusionName,jsonfile,op_count,subx_count))

def getOpStats(funcs):
    op_count = sum([count_ops(x['ret']) for x in funcs.values()])
    subx_count = len(funcs['subx']['ret']) if 'subx' in funcs else 0
    return op_count, subx_count
