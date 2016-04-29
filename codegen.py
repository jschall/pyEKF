from sympy import *
from sympy.printing.ccode import *
from helpers import *
import math

def generateCovariance(jsonfile, cfile):
    PP_O, PP_S = loadExprsFromJSON(jsonfile, ['PP_O', 'PP_S'])

    N = PP_O.rows
    PPidx = []
    for k in range((N**2-N)/2+N):
        r = int(math.floor((2*N+1-math.sqrt((2*N+1)*(2*N+1)-8*k))/2))
        c = int(k - N*r + r*(r-1)/2 + r)
        PPidx.append((r,c,k))

    PPc = ''

    for e in PP_S:
        PPc += '    '+ccode(e[1], assign_to=e[0])+'\n'

    PPc += '\n'

    for (r,c,k) in PPidx:
        if r <= 15 and c <= 15:
            PPc += '    '+ccode(PP_O[r,c], assign_to='nextP[%u]' % (k,))+'\n'

    PPc += '    if (stateIndexLim > 15) {\n'

    for (r,c,k) in PPidx:
        if (r > 15 or c > 15) and (r <= 21 and c <= 21):
            PPc += '        '+ccode(PP_O[r,c], assign_to='nextP[%u]' % (k,))+'\n'

    PPc += '        if (stateIndexLim > 21) {\n'

    for (r,c,k) in PPidx:
        if r > 21 or c > 21:
            PPc += '            '+ccode(PP_O[r,c], assign_to='nextP[%u]' % (k,))+'\n'

    PPc += '        }\n    }\n'

    f = open(cfile, 'w')
    f.write(PPc)
    f.close()
