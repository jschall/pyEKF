from sympy import *
from sympy.solvers import solve
from helpers import *
from sys import exit

# EKF to estimate target position and velocity relative to vehicle

# Goals:
# - Provide target height estimation (i.e. depth-from-motion) for cases where
#   there is no range finder or the target is on an elevated platform (like the
#   roof of a car)
# - Decouple precision landing performance from navigation performance and
#   potentially allow precision landings indoors.
# - Allow fusion of multiple sources of data relating to the target location
#   (e.g. if your computer vision algorithm provides an estimate of distance to
#   target based on its size in the frame)

# Parameters
dt = Symbol('dt')
Tbn = Matrix(3,3,symbols('Tbn[0:3][0:3]'))

vehicleVelNED_R = toVec(symbols('vv_n_R vv_e_R vv_d_R'))
targetCameraPos_R_sca = Symbol('cam_pos_R')
terrainDistD_R_sca = Symbol('gnd_dist_d_R')
targetDist_R_sca = Symbol("target_dist_R")
initFocLen = Symbol("foc_len_init")
initFocLen_R = Symbol("foc_len_init_R")

terrainDistD_R = toVec(terrainDistD_R_sca)
targetDist_R = toVec(targetDist_R_sca)
targetCameraPos_R = toVec(targetCameraPos_R_sca, targetCameraPos_R_sca)

vehicleDelVelNED_noise = toVec(symbols('vdv_n_noise vdv_e_noise vdv_d_noise'))

# Observations
terrainDistD = Symbol('gnd_dist_d')
targetCameraPos = toVec(symbols('cam_pos_x cam_pos_y'))
vehicleVelNED = toVec(symbols('vv_n vv_e vv_d'))
targetDist = Symbol("target_dist")

# Inputs
vehicleDeltaVelocityNED = toVec(symbols('dvv_n dvv_e dvv_d'))

# States
targetPosNED = toVec(symbols('pt_n pt_e pt_d'))
targetVelNED = toVec(symbols('vt_n vt_e vt_d'))
focLen = Symbol('foc_len')
#delVelBiasNE = toVec(symbols('dvv_n_b dvv_e_b'))
stateVector = toVec(targetPosNED, targetVelNED, focLen)
nStates = len(stateVector)

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P = copy_upper_to_lower_offdiagonals(P)

def deriveInitialization(jsonfile):
    print('Beginning initialization derivation')
    t1 = datetime.datetime.now()

    # Vector from rotation vector
    vecToTargetBody = toVec(targetCameraPos[0], targetCameraPos[1], 1.)
    vecToTargetNED = Tbn*vecToTargetBody

    # Derive initial state and covariance from observations
    initTargetPosNED = vecToTargetNED*terrainDistD/vecToTargetNED[2]
    initTargetVelNED = -vehicleVelNED

    # x_n: initial state
    x_n = toVec(initTargetPosNED, initTargetVelNED, initFocLen)

    assert x_n.shape == stateVector.shape

    # z: initialization measurement vector
    z = toVec(terrainDistD, targetCameraPos, vehicleVelNED, initFocLen)

    # R: covariance of additive noise on z
    R = diag(*toVec(terrainDistD_R, targetCameraPos_R, vehicleVelNED_R, initFocLen_R))

    # H: initialization measurement influence matrix
    H = x_n.jacobian(z)

    # P_n: initial covariance
    P_n = H*R*H.T

    assert P_n.shape == P.shape

    # Optimizations
    P_n = upperTriangularToVec(P_n)
    x_n,P_n,subx = extractSubexpressions([x_n,P_n],'subx',threshold=1)

    # Output generation
    funcParams = {'Tbn': Tbn, 'hgt': terrainDistD, 'hgt_R': terrainDistD_R, 'cam_pos': targetCameraPos, 'cam_pos_R': targetCameraPos_R_sca, 'init_foc_len': initFocLen, 'init_foc_len_R': initFocLen_R, 'vel': vehicleVelNED, 'vel_R': vehicleVelNED_R}

    funcs = {}

    funcs['subx'] = {}
    funcs['subx']['params'] = funcParams
    funcs['subx']['ret'] = toVec([x[1] for x in subx])
    funcs['subx']['retsymbols'] = toVec([x[0] for x in subx])

    funcParams = funcParams.copy()
    funcParams['subx'] = toVec([x[0] for x in subx])

    funcs['state'] = {}
    funcs['state']['params'] = funcParams
    funcs['state']['ret'] = x_n

    funcs['cov'] = {}
    funcs['cov']['params'] = funcParams
    funcs['cov']['ret'] = P_n

    check_funcs(funcs)

    saveExprsToJSON(jsonfile, {'funcs':serialize_exprs_in_structure(funcs.copy())})

    op_count, subx_count = getOpStats(funcs)
    t2 = datetime.datetime.now()
    print('%s initialization: derivation saved to %s. %u ops, %u subexpressions.' % (t2-t1,jsonfile,op_count,subx_count))

def derivePrediction(jsonfile):
    print('Beginning prediction derivation')
    t1 = datetime.datetime.now()

    # States at time k+1
    targetPosNEDNew = targetPosNED+targetVelNED*dt
    targetVelNEDNew = targetVelNED-vehicleDeltaVelocityNED
    focLenNew = focLen

    # f: state-transtition model
    f = toVec(targetPosNEDNew, targetVelNEDNew, focLenNew)
    assert f.shape == stateVector.shape

    # F: linearized state-transition model
    F = f.jacobian(stateVector)

    # u: control input vector
    u = toVec(vehicleDeltaVelocityNED)

    # G: control-influence matrix, AKA "B" in literature
    G = f.jacobian(u)

    # w_u_sigma: additive noise on u
    w_u_sigma = toVec(vehicleDelVelNED_noise)

    # Q_u: covariance of additive noise on u
    Q_u = diag(*w_u_sigma.multiply_elementwise(w_u_sigma))

    # Q: covariance of additive noise on x
    Q = G*Q_u*G.T

    # P_n: covariance matrix at time k+1
    P_n = F*P*F.T + Q
    assert P_n.shape == P.shape

    x_n = f

    # Optimizations
    P_n = upperTriangularToVec(P_n)
    x_n,P_n,subx = extractSubexpressions([x_n,P_n],'subx',threshold=10)

    # Output generation
    funcParams = {'x':stateVector,'P':upperTriangularToVec(P),'u':u,'w_u_sigma':w_u_sigma,'dt':dt}

    funcs = {}

    funcs['subx'] = {}
    funcs['subx']['params'] = funcParams
    funcs['subx']['ret'] = toVec([x[1] for x in subx])
    funcs['subx']['retsymbols'] = toVec([x[0] for x in subx])

    funcParams = funcParams.copy()
    funcParams['subx'] = toVec([x[0] for x in subx])

    funcs['state'] = {}
    funcs['state']['params'] = funcParams
    funcs['state']['ret'] = x_n

    funcs['cov'] = {}
    funcs['cov']['params'] = funcParams
    funcs['cov']['ret'] = P_n

    check_funcs(funcs)

    saveExprsToJSON(jsonfile, {'funcs':serialize_exprs_in_structure(funcs.copy())})

    op_count, subx_count = getOpStats(funcs)
    t2 = datetime.datetime.now()
    print('%s prediction: derivation saved to %s. %u ops, %u subexpressions.' % (t2-t1,jsonfile,op_count,subx_count))

def deriveCameraFusion(jsonfile):
    targetPosBody = Tbn.T*targetPosNED
    measPred = toVec(targetPosBody[0]/targetPosBody[2], targetPosBody[1]/targetPosBody[2])*focLen

    deriveFusionSimultaneous('camera', jsonfile, measPred, additionalinputs={'Tbn':Tbn})

def deriveVelNEFusion(jsonfile):
    measPred = -targetVelNED[0:2,0]

    deriveFusionSimultaneous('velNE',jsonfile,measPred)

def deriveVelDFusion(jsonfile):
    measPred = -targetVelNED[2:3,0]

    deriveFusionSimultaneous('velD',jsonfile,measPred)

def deriveHeightFusion(jsonfile):
    measPred = targetPosNED[2:3,0]

    deriveFusionSimultaneous('height',jsonfile,measPred)

def deriveFusionSimultaneous(fusionName,jsonfile,measPred,additionalinputs={},subs={},R_type='scalar'):
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

    # Apply specified substitutions
    y = y.xreplace(subs)
    NIS = NIS.xreplace(subs)
    x_n = x_n.xreplace(subs)
    P_n = P_n.xreplace(subs)
    temp_subx = [(x[0], x[1].xreplace(subs)) for x in temp_subx]

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


    saveExprsToJSON(jsonfile, {'funcs':serialize_exprs_in_structure(funcs.copy())})

    op_count, subx_count = getOpStats(funcs)
    t2 = datetime.datetime.now()
    print('%s %s fusion: derivation saved to %s. %u ops, %u subexpressions.' % (t2-t1,fusionName,jsonfile,op_count,subx_count))

def getOpStats(funcs):
    op_count = sum([count_ops(x['ret']) for x in funcs.values()])
    subx_count = len(funcs['subx']['ret']) if 'subx' in funcs else 0
    return op_count, subx_count
