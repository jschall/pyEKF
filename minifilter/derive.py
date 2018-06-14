import sys
sys.path.append('..')
from sympy import *
from helpers import *
import math

from symbols import *

P_symmetric = P
for r in range(P_symmetric.rows):
    for c in range(P_symmetric.cols):
        if r > c:
            P_symmetric[c,r] = P_symmetric[r,c]

def deriveZeroRotErr(jsonfile):
    # This represents the rotation from q to q_n - we will substitute in the rotation error states later
    refIncrRot = Matrix(symbols('refIncr[0:3]'))
    T_qn_q = rot_vec_to_matrix(refIncrRot)

    # f: state transtition model
    f = x[:,:]
    f[0:3,0] = T_qn_q.T * rotErr

    # F: linearized state-transition model
    F = f.jacobian(x)

    # q_n: reference quaternion after step
    q_n = quat_normalize(quat_multiply(q, rot_vec_to_quat(refIncrRot)))

    # x_n: state matrix after step
    x_n = x[:,:]
    x_n[0:3,0] -= refIncrRot

    # P_n: covariance matrix after step
    P_n = F*P*F.T

    # Substitute rotErr for refIncrRot - this will cause rotErr to be reset to zero
    subs = dict(zip(refIncrRot, rotErr))
    q_n = q_n.subs(subs)
    x_n = x_n.subs(subs)
    P_n = P_n.subs(subs)

    x_n, P_n, q_n, subexp = extractSubexpressions([x_n,P_n,q_n],'subexp')


    saveExprsToJSON(jsonfile, {'x_n': x_n, 'P_n':P_n, 'q_n': q_n, 'subexp':subexp})


def derivePrediction(jsonfile):
    errQuat = rot_vec_to_quat_approx(rotErr)
    truthQuat = quat_multiply(q, errQuat)
    Tbn = quat_to_matrix(truthQuat)
    velNEDNew = velNED+gravityNED*dt+Tbn*dVelMeas

    # f: state transtition model
    f = toVec(rotErr+dAngMeas-gyroBias*dt, velNEDNew, gyroBias)

    # F: linearized state-transition model
    F = f.jacobian(x)

    # G: control-influence matrix, AKA "B" in literature
    G = f.jacobian(u)

    # covariance of additive noise on control inputs (u)
    Q_u = diag(*w_u.multiply_elementwise(w_u))


    # covariance of additive noise on states (x)
    Q = G*Q_u*G.T

    pnoise = zeros(nStates,1)
    pnoise[6:9,0] = ones(3,1) * 0.01 * dt
    Q_pnoise = diag(*pnoise.multiply_elementwise(pnoise))

    # x_n: state vector after prediction step
    x_n = f[:,:]

    # P_n: covariance matrix after prediction step
    P_n = F*P*F.T + Q + Q_pnoise

    # q_n: reference quaternion after prediction step
    q_n = q

    # assume symmetry
    P_n = P_n.subs(zip(P, P_symmetric))

    # rotErr is known to always be zero prior to this step, substutute zero for computational optimization
    subs = dict(zip(rotErr, zeros(3,1)))
    q_n = q_n.subs(subs)
    x_n = x_n.subs(subs)
    P_n = P_n.subs(subs)

    x_n, P_n, q_n, subexp = extractSubexpressions([x_n,P_n,q_n],'subexp')

    saveExprsToJSON(jsonfile, {'Tbn': Tbn, 'x_n': x_n, 'P_n':P_n, 'q_n': q_n, 'subexp':subexp})

def deriveUpdate(jsonfile):
    # P_n: covariance matrix after update step
    P_n = P

    # x_n: state vector after update step
    x_n = x

    # h: predicted measurement
    h = velNED

    # y: measurement innovation
    y = velNEDMeas-h

    # sequential fusion
    for i in range(len(h)):
        # H: observation sensitivity matrix
        H = Matrix([[h[i]]]).jacobian(x)

        # K: near-optimal kalman gain
        K = (P*H.T)/(H*P*H.T + Matrix([[R_VEL]]))[0]

        x_n = x_n + K*y[i]

        P_n = (eye(nStates)-K*H)*P_n

    # q_n: reference quaternion after update step
    q_n = q

    # rotErr is known to always be zero prior to this step, substutute zero for computational optimization
    subs = dict(zip(rotErr, zeros(3,1)))
    q_n = q_n.subs(subs)
    x_n = x_n.subs(subs)
    P_n = P_n.subs(subs)

    x_n, P_n, q_n, subexp = extractSubexpressions([x_n,P_n,q_n],'subexp')

    saveExprsToJSON(jsonfile, {'x_n': x_n, 'P_n':P_n, 'q_n': q_n, 'subexp':subexp})

zerojson = 'minizero.json'
predictjson = 'minipredict.json'
updatejson = 'miniupdate.json'

derivePrediction(predictjson)
deriveUpdate(updatejson)
deriveZeroRotErr(zerojson)
