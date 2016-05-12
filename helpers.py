from sympy import *
import math
import itertools

def toVec(*args):
    ret = Matrix(map(lambda x: Matrix([x]), args)).vec()
    return ret

def skew(_v):
    v = toVec(_v)
    assert v.rows == 3

    return Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def rot_vec_to_quat(_v):
    v = toVec(_v)
    assert v.rows == 3

    theta = sqrt(v[0]**2+v[1]**2+v[2]**2)
    axis = v/theta
    return toVec(cos(theta/2.), sin(theta/2.) * axis[0], sin(theta/2.) * axis[1], sin(theta/2.) * axis[2])

def rot_vec_to_quat_approx(_v):
    v = toVec(_v)
    assert v.rows == 3

    return toVec(1.,v*0.5)

def quat_to_rot_vec_approx(_q):
    q = toVec(_q)
    assert q.rows == 4

    return 2.*toVec(q[1],q[2],q[3])

def quat_rotate_approx(_q, _v):
    return quat_multiply(_q,rot_vec_to_quat_approx(_v))

def quat_to_rot_vec(_q):
    q = toVec(_q)
    assert q.rows == 4

    theta = 2.*acos(q[0])
    axis = toVec(q[1],q[2],q[3])/sqrt(q[1]**2+q[2]**2+q[3]**2)

    return theta*axis

def quat_inverse(_q):
    q = toVec(_q)
    assert q.rows == 4

    q[1] = -q[1]
    q[2] = -q[2]
    q[3] = -q[3]
    return q

def quat_multiply(_q1, _q2):
    q1 = toVec(_q1)
    q2 = toVec(_q2)
    assert q1.rows == 4 and q2.rows == 4

    return toVec(q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
                 q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2],
                 q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1],
                 q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0])

def quat_normalize(_q):
    q = toVec(_q)
    assert q.rows == 4

    return q/sqrt(q[0]**2+q[1]**2+q[2]**2+q[3]**2)

def quat_to_matrix(_q):
    q = toVec(_q)
    assert q.rows == 4

    return (q[0]**2-(q[1:,0].T*q[1:,0])[0])*eye(3) + 2.*(q[1:,0]*q[1:,0].T) + 2.*q[0]*skew(q[1:,0])

def rot_vec_to_matrix(_v):
    v = toVec(_v)
    assert v.rows == 3

    theta = sqrt(v[0]**2+v[1]**2+v[2]**2)
    axis = v/theta
    return eye(3)*cos(theta)+(1-cos(theta))*axis*axis.T+skew(axis)*sin(theta)

def quickinv_sym(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols
    n = M.rows
    A = Matrix(n,n,symbols('_X[0:%u][0:%u]' % (n,n)))
    A = copy_upper_to_lower_offdiagonals(A)
    B = simplify(A.inv())
    return B.xreplace(dict(zip(A,M)))

def zero_lower_offdiagonals(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols

    ret = M[:,:]

    for r in range(ret.rows):
        for c in range(ret.cols):
            if r > c:
                ret[r,c] = 0.
    return ret

def copy_upper_to_lower_offdiagonals(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols

    ret = M[:,:]

    for r in range(ret.rows):
        for c in range(ret.cols):
            if r > c:
                ret[r,c] = ret[c,r]
    return ret

def count_subexpression(subexpr, expr):
    if hasattr(expr, "__getitem__"):
        return sum(map(lambda x: count_subexpression(subexpr, x), expr))
    else:
        return expr.count(subexpr)

def extractSubexpressions(inexprs, prefix='X',threshold=10):
    subexprs, outexprs = cse(inexprs)

    for i in reversed(range(len(subexprs))):
        ops_saved = (count_subexpression(subexprs[i][0], [[x[1] for x in subexprs], outexprs])-1)*subexprs[i][1].count_ops()
        if ops_saved < threshold:
            sub = dict([subexprs.pop(i)])
            subexprs = map(lambda x: (x[0],x[1].xreplace(sub)), subexprs)
            outexprs = map(lambda x: x.xreplace(sub), outexprs)

    for i in range(len(subexprs)):
        newSym = Symbol('%s[%u]' % (prefix,i))
        sub = {subexprs[i][0]:newSym}
        subexprs[i] = (newSym,subexprs[i][1])
        subexprs = map(lambda x: (x[0],x[1].xreplace(sub)), subexprs)
        outexprs = map(lambda x: x.xreplace(sub), outexprs)

    outexprs = map(lambda x: Matrix(x) if type(x) is ImmutableDenseMatrix else x, outexprs)

    return tuple(outexprs+[subexprs])

def pow_to_sq(string):
    import re
    return re.sub(r"pow\(([^,]+),\s*2\s*\)", 'sq(\g<1>)', string)

def serialize_exprs_in_structure(obj):
    if isinstance(obj, dict):
        return {k: serialize_exprs_in_structure(v) for k, v in obj.iteritems()}
    elif isinstance(obj, list):
        return [serialize_exprs_in_structure(x) for x in obj]
    else:
        return srepr(obj)

def deserialize_exprs_in_structure(obj):
    if isinstance(obj, dict):
        return {k: deserialize_exprs_in_structure(v) for k, v in obj.iteritems()}
    elif isinstance(obj, list):
        return [deserialize_exprs_in_structure(x) for x in obj]
    else:
        return sympify(obj)

def loadExprsFromJSON(fname):
    with open(fname, 'r') as f:
        import json
        imported = json.load(f)
        imported['operations'] = deserialize_exprs_in_structure(imported['operations'])
        return imported

def saveExprsToJSON(fname, input_dict):
    with open(fname, 'w') as f:
        f.truncate()
        import json
        output_dict = input_dict.copy()
        output_dict['operations'] = serialize_exprs_in_structure(output_dict['operations'])

        json.dump(output_dict, f)

def symmetricMatrixToVec(M):
    assert M.rows == M.cols

    N = M.rows
    r = lambda k: int(math.floor((2*N+1-math.sqrt((2*N+1)*(2*N+1)-8*k))/2))
    c = lambda k: int(k - N*r(k) + r(k)*(r(k)-1)/2 + r(k))
    return toVec([M[r(k),c(k)] for k in range((N**2-N)/2+N)])

def listSymbols(expr):
    if hasattr(expr, "__getitem__"):
        return set([item for sublist in map(lambda x: listSymbols(x), expr) for item in sublist])
    else:
        return expr.atoms(Symbol)

def check_funcs(funcs):
    for v in funcs.values():
        insymbols = listSymbols(v['params'].values())
        if 'retsymbols' in v:
            insymbols = insymbols.union(listSymbols(v['retsymbols']))
        funcsymbols = listSymbols(v['ret'])
        straysymbols = funcsymbols-insymbols
        assert not straysymbols, 'stray symbols: %s' % (str(straysymbols),)
