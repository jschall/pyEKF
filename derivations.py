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
# Parameterization: Target position is encoded as posNED = normalize(p_n, p_e, 1)/range_inv, and inverse range.
targetPos = toVec(symbols('p_n p_e range_inv'))
targetVelNED = toVec(symbols('vt_n vt_e vt_d'))
stateVector = toVec(targetPos, targetVelNED)
nStates = len(stateVector)

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P = copy_upper_to_lower_offdiagonals(P)

# Common computations
def posStateToPosNED(posState):
    posStateSym = toVec(symbols('_pos_state[0:3]'))

    ret = toVec(posStateSym[0], posStateSym[1], 1)
    ret /= vec_norm(ret)*posStateSym[2]

    return simplify(ret).xreplace(dict(zip(posStateSym, posState)))

def posNEDToPosState(posNED):
    posNEDSym = toVec(symbols('_pos_ned[0:3]'))
    posStateSym = toVec(symbols('_pos_state[0:3]'))

    posNEDEqn = posStateToPosNED(posStateSym)

    return simplify(toVec(solve(posNEDEqn-posNEDSym, posStateSym)[0])).xreplace(dict(zip(posNEDSym, posNED)))

def velNEDToPosStateDeriv(velNED, posState):
    t = Symbol('_t')
    velNEDSym = toVec(symbols('_vel_ned[0:3]'))
    posStateSym = toVec(symbols('_pos_state[0:3]'))

    return simplify(posNEDToPosState(posStateToPosNED(posStateSym)+velNEDSym*t).diff(t).subs(t,0)).xreplace(dict(zip(velNEDSym, velNED))).xreplace(dict(zip(posStateSym, posState)))

def posStateDerivToVelNED(posStateDeriv, posState):
    t = Symbol('_t')
    posStateDerivSym = toVec(symbols('_pos_state_deriv[0:3]'))
    posStateSym = toVec(symbols('_pos_state[0:3]'))

    return simplify(posStateToPosNED(posStateSym+posStateDeriv*t).diff(t).subs(t,0)).xreplace(dict(zip(posStateDerivSym, posStateDeriv))).xreplace(dict(zip(posStateSym, posState)))

targetPosNED = posStateToPosNED(targetPos)
#targetVelNED = posStateDerivToVelNED(targetPosDeriv, targetPos)

def deriveConversions(jsonfile):
    print('Beginning conversions derivation')
    t1 = datetime.datetime.now()

    posNEDSym = toVec(symbols('_pos_ned[0:3]'))
    posStateSym = toVec(symbols('_pos_state[0:3]'))

    z = targetPos
    x = posStateToPosNED(z)
    R = P[0:3,0:3]
    H = x.jacobian(z)
    cov = H*R*H.T

    funcs = {}

    funcs['pos_ned'] = {}
    funcs['pos_ned']['params'] = {'x': stateVector}
    funcs['pos_ned']['ret'] = posStateToPosNED(targetPos)

    funcs['pos_ned_var'] = {}
    funcs['pos_ned_var']['params'] = {'x': stateVector, 'P':upperTriangularToVec(P)}
    funcs['pos_ned_var']['ret'] = toVec([cov[i,i] for i in range(3)])


    check_funcs(funcs)

    saveExprsToJSON(jsonfile, {'funcs':serialize_exprs_in_structure(funcs.copy())})

    op_count, subx_count = getOpStats(funcs)
    t2 = datetime.datetime.now()
    print('%s conversions: derivation saved to %s. %u ops, %u subexpressions.' % (t2-t1,jsonfile,op_count,subx_count))

def deriveInitialization(jsonfile):
    print('Beginning initialization derivation')
    t1 = datetime.datetime.now()

    depthInvInit, depthInvInitSigma = symbols('depth_inv_init depth_inv_init_sigma')

    # Vector from rotation vector
    unitvecToTargetBody = toVec(targetCameraPos[0], targetCameraPos[1], 1.)
    unitvecToTargetBody /= vec_norm(unitvecToTargetBody)
    unitvecToTargetNED = Tbn*unitvecToTargetBody

    # Derive initial state and covariance from observations
    initTargetPosNED = unitvecToTargetNED/depthInvInit
    initTargetPos = posNEDToPosState(initTargetPosNED)
    initTargetVelNED = -vehicleVelNED
    #initTargetPosDeriv = velNEDToPosStateDeriv(initTargetVelNED, initTargetPos)

    # x_n: initial state
    x_n = toVec(initTargetPos, initTargetVelNED)

    assert x_n.shape == stateVector.shape

    # z: initialization measurement vector
    z = toVec(depthInvInit, targetCameraPos, vehicleVelNED)

    # R: covariance of additive noise on z
    R = diag(*toVec(depthInvInitSigma**2, targetCameraPos_R, vehicleVelNED_R))

    # H: initialization measurement influence matrix
    H = x_n.jacobian(z)

    # P_n: initial covariance
    P_n = H*R*H.T

    assert P_n.shape == P.shape

    # Compute depth initialization
    # This puts the 95% confidence interval between 0.2m and infinity
    subs = solve((1/(depthInvInit+2*depthInvInitSigma)-0.2, depthInvInit-2*depthInvInitSigma), (depthInvInit, depthInvInitSigma))
    pprint(subs)
    x_n = x_n.xreplace(subs)
    P_n = P_n.xreplace(subs)

    # Optimizations
    P_n = upperTriangularToVec(P_n)
    x_n,P_n,subx = extractSubexpressions([x_n,P_n],'subx',threshold=10)

    # Output generation
    funcParams = {'Tbn': Tbn, 'cam_pos': targetCameraPos, 'cam_pos_R': targetCameraPos_R_sca, 'vel': vehicleVelNED, 'vel_R': vehicleVelNED_R}

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
    #targetPosDeriv = velNEDToPosStateDeriv(targetVelNED, targetPos)
    targetPosNEDNew = targetPosNED+targetVelNED*dt
    targetPosNew = posNEDToPosState(targetPosNEDNew)
    #targetPosNew = targetPos+targetPosDeriv*dt
    targetVelNEDNew = targetVelNED-vehicleDeltaVelocityNED
    #targetPosDerivNew = velNEDToPosStateDeriv(targetVelNEDNew, targetPos)

    # f: state-transtition model
    f = simplify(toVec(targetPosNew, targetVelNEDNew))
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
    measPred = toVec(targetPosBody[0]/targetPosBody[2], targetPosBody[1]/targetPosBody[2])

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
