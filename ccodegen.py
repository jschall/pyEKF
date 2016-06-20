from sympy import *
from sympy.printing.ccode import *
from helpers import *
import math

class MyPrinter(CCodePrinter):
    def __init__(self,settings={}):
        CCodePrinter.__init__(self, settings)
        self.known_functions = {
            "Abs": [(lambda x: not x.is_integer, "fabsf")],
            "gamma": "tgammaf",
            "sin": "sinf",
            "cos": "cosf",
            "tan": "tanf",
            "asin": "asinf",
            "acos": "acosf",
            "atan": "atanf",
            "atan2": "atan2f",
            "exp": "expf",
            "log": "logf",
            "erf": "erff",
            "sinh": "sinhf",
            "cosh": "coshf",
            "tanh": "tanhf",
            "asinh": "asinhf",
            "acosh": "acoshf",
            "atanh": "atanhf",
            "floor": "floorf",
            "ceiling": "ceilf",
        }

    def _print_Pow(self, expr):
        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)

        if expr.exp == -1:
            return '1.0/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp < 0:
            expr = 1/expr

        if expr.exp == 0.5:
            return 'sqrtf(%s)' % self._print(expr.base)
        elif expr.exp.is_integer and expr.exp <= 4:
            return '*'.join([self._print(expr.base) for _ in range(expr.exp)])
        else:
            return 'powf(%s, %s)' % (self._print(expr.base),
                                 self._print(expr.exp))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' % (p, q)

def generateCode(jsondict, cfile):
    # extract filter operations
    filterOps = {}
    for n,fn in jsondict.iteritems():
        filterOps[n.upper()] = deserialize_exprs_in_structure(loadExprsFromJSON(fn)['funcs'])

    hashdefines = []
    hashdefines.extend(getConstants(filterOps))
    for k in sorted(filterOps.keys()):
        hashdefines.extend(getSnippetDefines(k,filterOps[k]))

    with open(cfile, 'w') as f:
        f.truncate()
        f.write(getHeader(hashdefines, 'EKF_'))
        print('Generated code saved to %s'%(cfile,))

def getConstants(filterOps):
    max_num_subx = 0
    for k,v in filterOps.iteritems():
        if 'subx' in v:
            max_num_subx = max(len(v['subx']['ret']),max_num_subx)

    x = filterOps['PREDICTION']['cov']['params']['x']
    u = filterOps['PREDICTION']['cov']['params']['u']
    ret = []
    ret.append(('NUM_STATES', x.rows))
    ret.append(('NUM_CONTROL_INPUTS', u.rows))
    ret.append(('MAX_NUM_SUBX', max_num_subx))
    for r in range(x.rows):
        ret.append(('STATE_IDX_'+str(x[r,0]).upper(), r))
    for r in range(u.rows):
        ret.append(('U_IDX_'+str(u[r,0]).upper(), r))
    return ret

def getSnippetDefines(opname, funcs):
    ret = []

    for retName,func in funcs.iteritems():
        retParamName = '__RET_'+retName.upper()

        paramlist = []
        substitutions = []
        for paramname in sorted(func['params'].keys()):
            if not isinstance(func['params'][paramname], MatrixBase):
                func['params'][paramname] = Matrix([[func['params'][paramname]]])

            nr,nc = func['params'][paramname].shape

            paramlist.append('__'+paramname.upper())

            if (nr,nc) == (1,1):
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, [Symbol(paramlist[-1])]))
            elif nc == 1:
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, symbols(paramlist[-1]+'[0:%u]'%(nr,))))
            else:
                substitutions += zip(func['params'][paramname], Matrix(nr,nc, symbols(paramlist[-1]+'[0:%u][0:%u]'%(nr,nc))))

        if 'retsymbols' in func:
            if not isinstance(func['retsymbols'], MatrixBase):
                func['retsymbols'] = Matrix([[func['retsymbols']]])

            nr,nc = func['retsymbols'].shape

            if (nr,nc) == (1,1):
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, [Symbol(retParamName)]))
            elif nc == 1:
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, symbols(retParamName+'[0:%u]'%(nr,))))
            else:
                substitutions += zip(func['retsymbols'], Matrix(nr,nc, symbols(retParamName+'[0:%u][0:%u]'%(nr,nc))))

        func['ret'] = func['ret'].xreplace(dict(substitutions))

        defineName = '%s_CALC_%s(%s)' % (opname,retName.upper(),', '.join(paramlist+[retParamName]))
        ret.append((defineName,getSnippet(retParamName,func['ret'])))

    return ret

def getSnippet(retParamName, outputMatrix):
    if not isinstance(outputMatrix, MatrixBase):
        outputMatrix = Matrix([[outputMatrix]])

    outputMatrix = outputMatrix

    nr,nc = outputMatrix.shape

    if (nr,nc) == (1,1):
        retMatrix = Matrix(nr,nc, [Symbol(retParamName)])
    elif nc == 1:
        retMatrix = Matrix(nr,nc, symbols(retParamName+'[0:%u]'%(nr,)))
    else:
        retMatrix = Matrix(nr,nc, symbols(retParamName+'[0:%u][0:%u]'%(nr,nc)))

    ret = ''
    for assignee,expr in zip(retMatrix,outputMatrix):
        ret += double2float(MyPrinter().doprint(expr, assignee))+' '

    return str(ret).replace('\n','')

def double2float(string):
    import re
    string = re.sub(r"[0-9]+\.[0-9]+", '\g<0>f', string)

    return string

def wrapstring(string, linemax, delim):
    strings = string.split(' ')
    lines = [strings[0]]
    for s in strings[1:]:
        if len(lines[-1]) != 0 and len(lines[-1])+len(s) > linemax:
            lines.append(s)
        else:
            lines[-1] += ' '+s
    return delim.join(lines)

def getHeader(hashdefines, prefix=''):
    ret = '/*\n'
    for (key,val) in hashdefines:
        if type(val) is str:
            ret += prefix+key+'\n'
    ret += '*/\n\n'
    for (key,val) in hashdefines:
        if type(val) is str:
            val = wrapstring(str(val),100,' \\\n')
            ret += '\n#define %s%s \\\n%s\n' % (prefix,key,val)
        else:
            ret += '#define %s%s %s\n' % (prefix,key,val)
    return ret
