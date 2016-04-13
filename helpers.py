from sympy import *
import itertools

def toVec(*args):
    ret = Matrix(map(lambda x: Matrix([x]), args)).vec()
    return ret

def rot_vec_to_quat(_v):
    v = toVec(_v)
    theta = v.norm()
    axis = v/theta
    return toVec(cos(theta/2.), sin(theta/2.) * axis[0], sin(theta/2.) * axis[1], sin(theta/2.) * axis[2])

def rot_vec_to_quat_approx(_v):
    v = toVec(_v)
    return toVec(1.,v*0.5)

def quat_to_rot_vec(_q):
    q = toVec(_q)
    return toVec(q[1],q[2],q[3])*q[0]/toVec(q[1],q[2],q[3]).norm()

def quat_to_rot_vec_approx(_q):
    q = toVec(_q)
    return 2.*toVec(q[1],q[2],q[3])

def quat_to_matrix(_q):
    q = toVec(_q)
    return Matrix([[ 1.-2.*(q[2]**2+q[3]**2), 2.*(q[1]*q[2]-q[0]*q[3]), 2.*(q[1]*q[3]+q[0]*q[2])],
                   [2.*(q[1]*q[2]+q[0]*q[3]),  1.-2.*(q[1]**2+q[3]**2), 2.*(q[2]*q[3]-q[0]*q[1])],
                   [2.*(q[1]*q[3]-q[0]*q[2]), 2.*(q[2]*q[3]+q[0]*q[1]),  1.-2.*(q[1]**2+q[2]**2)]])

def quat_inverse(_q):
    q = toVec(_q)
    q[1] = -q[1]
    q[2] = -q[2]
    q[3] = -q[3]
    return q

def quat_multiply(_q1, _q2):
    q1 = toVec(_q1)
    q2 = toVec(_q2)
    return toVec(q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
                 q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2],
                 q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1],
                 q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0])

def optimizeAlgebra(inexpr, prefix='X'):
    def nameGen():
        n=0
        while True:
            yield Symbol('%s[%u]' % (prefix,n))
            n+=1

    inexpr = Matrix(simplify(inexpr))
    subexpr, outexpr = cse(inexpr, nameGen(), optimizations='basic')
    if len(outexpr) == 1:
        outexpr = outexpr[0]
    return outexpr, subexpr#[x[1] for x in subexpr]
