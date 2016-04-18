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
    return Matrix([[ q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2,             2.*(q[1]*q[2] - q[0]*q[3]),              2.*(q[1]*q[3] + q[0]*q[2])],
                   [            2.*(q[1]*q[2] + q[0]*q[3]),  q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2,              2.*(q[2]*q[3] - q[0]*q[1])],
                   [              2.*(q[1]*q[3]-q[0]*q[2]),             2.*(q[2]*q[3] + q[0]*q[1]),  q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2]])

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

def optimizeAlgebra(inexpr, prefix='X',threshold=10):
    debug=True
    #inexpr = Matrix(simplify(inexpr))
    subexprs, outexprs = cse(inexpr)
    if len(outexprs) == 1:
        outexprs = outexprs[0]

    if debug:
        print "total ops before cse:",sum(map(lambda x: x[1].count_ops(), subexprs))+sum(map(lambda x: x.count_ops(), inexpr))
        total_ops_before = sum(map(lambda x: x[1].count_ops(), subexprs))+sum(map(lambda x: x.count_ops(), outexprs))
        print "total ops after cse:", total_ops_before

    # eliminate subexpressions that are below the threshold
    i = 0
    while i < len(subexprs):
        subexprs[i] = simplify(subexprs[i])
        subexpr_occurances = sum(map(lambda x: x[1].count(subexprs[i][0]), subexprs))
        outexpr_occurances = sum(map(lambda x: x.count(subexprs[i][0]), outexprs))
        occurances = subexpr_occurances+outexpr_occurances
        ops = subexprs[i][1].count_ops()
        ops_saved = occurances*ops-ops
        if ops_saved < threshold:
            sub = subexprs.pop(i)
            subexprs = map(lambda x: (x[0],x[1].subs(*sub)), subexprs)
            outexprs = outexprs.subs(*sub)
            continue
        i += 1

    # TODO: could search for discarded subexpressions that now exceed threshold and restore them

    for i in range(len(subexprs)):
        newSym = Symbol('%s[%u]' % (prefix,i))
        sub = (subexprs[i][0],newSym)
        subexprs[i] = (newSym,subexprs[i][1])
        subexprs = map(lambda x: (x[0],x[1].subs(*sub)), subexprs)
        outexprs = outexprs.subs(*sub)

    subexprs = map(lambda x: (x[0],simplify(x[1])), subexprs)
    outexprs = simplify(outexprs)

    if debug:
        total_ops_after = sum(map(lambda x: x[1].count_ops(), subexprs))+sum(map(lambda x: x.count_ops(), outexprs))

        print "total ops after cse elimination",total_ops_after

        ops_saved_total = total_ops_before-total_ops_after

        print "total ops: %u" % (sum(map(lambda x: x[1].count_ops(), subexprs))+sum(map(lambda x: x.count_ops(), outexprs)),)
        print "subexpressions: %u, total saved ops: %u" % (len(subexprs),ops_saved_total)

    return outexprs, subexprs

def pow_to_sq(string):
    import re
    return re.sub(r"pow\(([^,]+),\s*2\s*\)", 'sq(\g<1>)', string)

def row_column_to_linear_symmetric(N,r,c):
    return (c-r)+24*r-r*(r-1)/2

