from derivations import *
from codegen import *
import os
import argparse
from multiprocessing import Process
import datetime

print 'Begin %s' % (datetime.datetime.now().time().isoformat(),)

parser = argparse.ArgumentParser()
parser.add_argument('--outputdir', nargs=1, dest='outdir', default=['output'])
parser.add_argument('--derive', nargs='*', dest='derive', default=[])
parser.add_argument('--codegen', dest='codegen', action='store_true')
args = parser.parse_args()

outdir = args.outdir[0]

if not os.path.exists(outdir):
    os.makedirs(outdir)

predictionjson = os.path.join(outdir, 'covariancePrediction.json')
posnejson = os.path.join(outdir, 'posNEFusion.json')
posdjson = os.path.join(outdir, 'posDFusion.json')
velnejson = os.path.join(outdir, 'velNEFusion.json')
veldjson = os.path.join(outdir, 'velDFusion.json')
airspeedjson = os.path.join(outdir, 'airspeedFusion.json')
betajson = os.path.join(outdir, 'betaFusion.json')
magjson = os.path.join(outdir, 'magFusion.json')
yaw312json = os.path.join(outdir, 'yaw312Fusion.json')
yaw321json = os.path.join(outdir, 'yaw321Fusion.json')
flowjson = os.path.join(outdir, 'flowFusion.json')
declinationjson = os.path.join(outdir, 'declinationFusion.json')
c_header = os.path.join(outdir, 'ekf_defines.h')

derivations = {
    'posne': (derivePosNEFusion, posnejson),
    'posd': (derivePosDFusion, posdjson),
    'velne': (deriveVelNEFusion, velnejson),
    'veld': (deriveVelDFusion, veldjson),
    'covpred': (deriveCovariancePrediction, predictionjson),
    'airspeed': (deriveAirspeedFusion, airspeedjson),
    'beta': (deriveBetaFusion, betajson),
    'mag': (deriveMagFusion, magjson),
    'yaw312': (deriveYaw312Fusion, yaw312json),
    'yaw321': (deriveYaw321Fusion, yaw321json),
    'flow': (deriveOptFlowFusion, flowjson),
    'declination': (deriveDeclinationFusion, declinationjson)
    }

assert set(args.derive).issubset(set(derivations.keys()+['all']))

if 'all' in args.derive:
    args.derive = derivations.keys()

workers = []
for k in args.derive:
    workers.append(Process(target=derivations[k][0], args=(derivations[k][1],)))
    workers[-1].start()

for p in workers:
    p.join()

if args.codegen:
    jsondict = {}
    for key, val in derivations.iteritems():
        jsondict[key] = val[1]

    generateCode(jsondict,c_header)

print 'End %s' % (datetime.datetime.now().time().isoformat(),)
