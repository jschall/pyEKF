import sys
sys.path.append('..')

from sympy import *
from helpers import *

import symbols
from random import gauss
import math

class Filter:
    def __init__(self, predictjson, updatejson, w_u):
        self.gravity = 9.81
        self.w_u = w_u
        self.q = toVec(1.,0.,0.,0.)
        self.x = toVec(0.0000001,0.,0.,0.,0.,0.)
        self.P = diag(.2**2,.2**2,.2**2,.1**2,.1**2,.1**2)
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

    def predict(self, dt, u):
        u = toVec(u)
        # compute subexpressions
        subexprs = [0. for _ in self.predict_subexp]
        for i in range(len(subexprs)):
            subexprs[i] = self.predict_subexp[i](subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u)

        self.x = Matrix(self.predict_x(subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u))
        self.P = Matrix(self.predict_P(subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u))
        self.q = Matrix(self.predict_q(subexprs,self.gravity,self.w_u,self.q,self.x,self.P,dt,u))

    def update(self, velNEDMeas):
        velNEDMeas = toVec(velNEDMeas)

        # compute subexpressions
        subexprs = [0. for _ in self.update_subexp]
        for i in range(len(subexprs)):
            subexprs[i] = self.update_subexp[i](subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas)

        self.x = Matrix(self.update_x(subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas))
        self.P = Matrix(self.update_P(subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas))
        self.q = Matrix(self.update_q(subexprs,self.q,self.x,self.P,self.R_VEL,velNEDMeas))

dt = 0.01
w_u = toVec(.03,.03,.03,.25,.25,.25)*dt
f = Filter('minipredict.json', 'miniupdate.json', w_u)

frames = []

t = 0.
while t < 20:
    t += dt
    imunoise = toVec(map(lambda x: gauss(0.,x),w_u))
    f.predict(dt, toVec(0.0,0.0,math.radians(2000)*dt,0.0,0.0,-9.81*dt)+imunoise)
    velnoise = toVec(gauss(0.,0.5),gauss(0.,0.5),gauss(0.,0.5))
    f.update(velnoise)
    print t
    pprint(f.x)
    pprint(f.P)
    pprint(f.q)
    frames.append((t,f.x,f.P,f.q))
    print "\n\n"

import matplotlib.pyplot as plt

plt.plot([x[0] for x in frames], [x[1][3] for x in frames])
plt.plot([x[0] for x in frames], [x[1][4] for x in frames])
plt.plot([x[0] for x in frames], [x[1][5] for x in frames])
plt.show()
