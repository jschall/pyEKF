from sympy import *
from helpers import *
from sympy.printing.ccode import *
from sys import exit
import math
horizontalAccelBiases = False

# Parameters
estQuat = toVec(symbols('q0 q1 q2 q3'))
dAngMeas = toVec(symbols('dax day daz'))
dVelMeas = toVec(symbols('dvx dvy dvz'))
dAngNoise = toVec(symbols('daxNoise dayNoise dazNoise'))
dVelNoise = toVec(symbols('dvxNoise dvyNoise dvzNoise'))
gravityNED = toVec(0.,0.,symbols('gravity'))
dt = Symbol('dt')

# States
rotErr = toVec(symbols('rotErrX rotErrY rotErrZ'))
vel = toVec(symbols('vn ve vd'))
pos = toVec(symbols('pn pe pd'))
dAngBias = toVec(symbols('dax_b day_b daz_b'))
dAngScale = toVec(symbols('dax_s day_s daz_s'))
if horizontalAccelBiases:
    dVelBias = toVec(symbols('dvx_b dvy_b dvz_b'))
else:
    dVelBias = toVec(0,0,symbols('dvz_b'))
magBody = toVec(symbols('magX magY magZ'))
magNED = toVec(symbols('magN magE magD'))
vwn, vwe = symbols('vwn vwe')
stateVector = toVec(rotErr,vel,pos,dAngBias,dAngScale,dVelBias if horizontalAccelBiases else dVelBias[2],magNED,magBody,vwn,vwe)
nStates = len(stateVector)
print "nStates=%u" % (nStates,)

# Covariance matrix
P = Matrix(nStates,nStates,symbols('P[0:%u][0:%u]' % (nStates,nStates)))
P_sym = P
for r in range(P_sym.rows):
    for c in range(P_sym.cols):
        if r > c:
            P_sym[c,r] = P_sym[r,c]

# Prediction step
errQuat = rot_vec_to_quat_approx(rotErr)
truthQuat = quat_multiply(estQuat, errQuat)
Tbn = quat_to_matrix(truthQuat)
dAngTruth = dAngMeas.multiply_elementwise(dAngScale) - dAngBias - dAngNoise
deltaQuat = rot_vec_to_quat_approx(dAngTruth)
truthQuatNew = quat_multiply(truthQuat, deltaQuat)
errQuatNew = quat_multiply(quat_inverse(estQuat), truthQuatNew)
rotErrNew = quat_to_rot_vec_approx(errQuatNew)
dVelTruth = dVelMeas - dVelBias - dVelNoise
velNew = vel+gravityNED*dt + Tbn*dVelTruth
posNew = pos+vel*dt
dAngBiasNew = dAngBias
dAngScaleNew = dAngScale
dVelBiasNew = dVelBias
magBodyNew = magBody
magNEDNew = magNED
vwnNew = vwn
vweNew = vwe

stateVectorNew = toVec(rotErrNew,velNew,posNew,dAngBiasNew,dAngScaleNew,dVelBiasNew if horizontalAccelBiases else dVelBiasNew[2], magNEDNew, magBodyNew, vwnNew, vweNew)

F = stateVectorNew.jacobian(stateVector)
F = F.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])
distVector = toVec(dAngNoise, dVelNoise)

G = stateVectorNew.jacobian(distVector)
G = G.subs([(rotErr[0], 0.),(rotErr[1], 0.),(rotErr[2], 0.)])

distMatrix = diag(*distVector.multiply_elementwise(distVector))
Q = G*distMatrix*G.T

# NOTE: PP has been analytically verified to be equal to PP in the priseoborough/InertialNav matlab code
PP = F*P*F.T+Q

# Optimizations

# assume that the P matrix is symmetrical
for r in range(PP.rows):
    for c in range(PP.cols):
        PP = PP.subs(P[r,c], P_sym[r,c])

# zero the lower off-diagonals
for r in range(PP.rows):
    for c in range(PP.cols):
        if r > c:
            PP[r,c] = 0.

PP,SPP = optimizeAlgebra(PP,'SPP')

# Code generation

PPc = ''
PPc += 'float SPP[%u];\n\n' % (len(SPP),)

for e in SPP:
    PPc += ccode(e[1], assign_to=e[0])+'\n'

N = PP.rows
PPc += '\nfloat nextP[%u];\n\n' % ((N**2-N)/2+N,)

PPidx = [0 for i in range((N**2-N)/2+N)]

for k in range(len(PPidx)):
    r = int(math.floor((2*N+1-math.sqrt((2*N+1)*(2*N+1)-8*k))/2))
    c = int(k - N*r + r*(r-1)/2 + r)
    PPidx[k] = (r,c,k)

for (r,c,k) in PPidx:
    if r <= 15 and c <= 15:
        PPc += ccode(PP[r,c], assign_to='nextP[%u]' % (k,))+'\n'

PPc += '\n'

for (r,c,k) in PPidx:
    if (r > 15 or c > 15) and (r <= 21 and c <= 21):
        PPc += ccode(PP[r,c], assign_to='nextP[%u]' % (k,))+'\n'

PPc += '\n'

for (r,c,k) in PPidx:
    if r > 21 or c > 21:
        PPc += ccode(PP[r,c], assign_to='nextP[%u]' % (k,))+'\n'

PPc += '\n'

# change pow(x,2) to sq(x)
PPc = pow_to_sq(PPc)

f = open('PP.c', 'w')
f.write(PPc)
f.close()
