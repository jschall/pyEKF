import sys
sys.path.append('..')

from sympy import *
from helpers import *

import symbols
from random import gauss

class Filter:
    def __init__(self, predictjson, updatejson, zeroRotErrjson):
        self.gravity = 9.81
        self.w_u = toVec(.0003,.0003,.0003,.0025,.0025,.0025)
        self.q = toVec(1.,0.,0.,0.)
        self.x = toVec(0.,0.,0.,0.,0.,0.)
        self.P = diag(.1**2,.1**2,.1**2,.1**2,.1**2,.1**2)
        self.R_VEL = 0.5

        self.predict_x, self.predict_P, self.predict_q, self.predict_subexp = loadExprsFromJSON(predictjson, ('x_n', 'P_n', 'q_n', 'subexp'))
        self.predict_subexp_symbols = [x[0] for x in self.predict_subexp]
        self.predict_args = [self.predict_subexp_symbols, symbols.gravity, symbols.w_u, symbols.q, symbols.x, symbols.P, symbols.dt, symbols.u]

        self.predict_subexp = [lambdify(self.predict_args,x[1]) for x in self.predict_subexp]
        self.predict_x = lambdify(self.predict_args,self.predict_x)
        self.predict_P = lambdify(self.predict_args,self.predict_P)
        self.predict_q = lambdify(self.predict_args,self.predict_q)

        self.update_x, self.update_P, self.update_q, self.update_subexp = loadExprsFromJSON(updatejson, ('x_n', 'P_n', 'q_n', 'subexp'))
        self.update_subexp_symbols = [x[0] for x in self.update_subexp]
        self.update_args = [self.update_subexp_symbols, symbols.q, symbols.x, symbols.P, symbols.R_VEL, symbols.velNEDMeas]

        self.update_subexp = [lambdify(self.update_args,x[1]) for x in self.update_subexp]
        self.update_x = lambdify(self.update_args,self.update_x)
        self.update_P = lambdify(self.update_args,self.update_P)
        self.update_q = lambdify(self.update_args,self.update_q)

        self.zeroRotErr_x, self.zeroRotErr_P, self.zeroRotErr_q, self.zeroRotErr_subexp = loadExprsFromJSON(zeroRotErrjson, ('x_n', 'P_n', 'q_n', 'subexp'))
        self.zeroRotErr_subexp_symbols = [x[0] for x in self.zeroRotErr_subexp]
        self.zeroRotErr_args = [self.zeroRotErr_subexp_symbols, symbols.q, symbols.x, symbols.P]

        self.zeroRotErr_subexp = [lambdify(self.zeroRotErr_args,x[1]) for x in self.zeroRotErr_subexp]
        self.zeroRotErr_x = lambdify(self.zeroRotErr_args,self.zeroRotErr_x)
        self.zeroRotErr_P = lambdify(self.zeroRotErr_args,self.zeroRotErr_P)
        self.zeroRotErr_q = lambdify(self.zeroRotErr_args,self.zeroRotErr_q)

    def zeroRotErr(self):
        subexprs = [0. for _ in self.zeroRotErr_subexp]
        for i in range(len(subexprs)):
            subexprs[i] = self.zeroRotErr_subexp[i](subexprs,self.q,self.x,self.P)

        x = Matrix(self.zeroRotErr_x(subexprs,self.q,self.x,self.P))
        P = Matrix(self.zeroRotErr_P(subexprs,self.q,self.x,self.P))
        q = Matrix(self.zeroRotErr_q(subexprs,self.q,self.x,self.P))

        self.x = x
        self.P = P
        self.q = q

    def predict(self, dt, u):
        u = toVec(u)
        # compute subexpressions
        subexprs = [0. for _ in self.predict_subexp]
        for i in range(len(subexprs)):
            subexprs[i] = self.predict_subexp[i](subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u)

        x = Matrix(self.predict_x(subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u))
        P = Matrix(self.predict_P(subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u))
        #q = Matrix(self.predict_q(subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u))

        self.x = x
        self.P = P
        #self.q = q

    def update(self, velNEDMeas):
        velNEDMeas = toVec(velNEDMeas)

        # compute subexpressions
        subexprs = [0. for _ in self.update_subexp]
        for i in range(len(subexprs)):
            subexprs[i] = self.update_subexp[i](subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas)

        x = Matrix(self.update_x(subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas))
        P = Matrix(self.update_P(subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas))
        #q = Matrix(self.update_q(subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas))

        self.x = x
        self.P = P
        #self.q = q

f = Filter('minipredict.json', 'miniupdate.json', 'minizero.json')

#from visual import *

for i in range(1000):
    imunoise = toVec(map(lambda x: gauss(0.,x),toVec(.003,.003,.003,.025,.025,.025)))
    f.predict(0.01, toVec(0.0,0.0,0.3,0.0,0.0,-0.0981)+imunoise)
    f.zeroRotErr()
    velnoise = toVec(gauss(0.,2),gauss(0.,2),gauss(0.,2))
    f.update(velnoise)
    f.zeroRotErr()
    pprint(f.x)
    pprint(f.P)
    pprint(f.q)
    print "\n\n"
