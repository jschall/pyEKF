from derivations import *
from codegen import *
import os
import argparse

derivations = ['prediction', 'airspeed', 'beta', 'mag', 'yaw', 'flow', 'declination']

parser = argparse.ArgumentParser()
parser.add_argument('--derive-all', dest='derive_all', action='store_true')
for x in derivations:
    parser.add_argument('--derive-%s'%(x,), dest='derive_%s'%(x,), action='store_true')
parser.add_argument('--codegen', dest='codegen', action='store_true')
args = parser.parse_args()

resultsdir = 'results'

if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)

predictionjson = os.path.join(resultsdir, 'covariancePrediction.json')
airspeedjson = os.path.join(resultsdir, 'airspeedFusion.json')
betajson = os.path.join(resultsdir, 'betaFusion.json')
magjson = os.path.join(resultsdir, 'magFusion.json')
yaw312json = os.path.join(resultsdir, 'yaw312Fusion.json')
yaw321json = os.path.join(resultsdir, 'yaw321Fusion.json')
flowjson = os.path.join(resultsdir, 'flowFusion.json')
declinationjson = os.path.join(resultsdir, 'declinationFusion.json')


if args.derive_all or args.derive_prediction:
    deriveCovariancePrediction(predictionjson)

if args.derive_all or args.derive_airspeed:
    deriveAirspeedFusion(airspeedjson)

if args.derive_all or args.derive_beta:
    deriveBetaFusion(betajson)

if args.derive_all or args.derive_mag:
    deriveMagFusion(magjson)

if args.derive_all or args.derive_yaw:
    deriveYaw312Fusion(yaw312json)
    deriveYaw321Fusion(yaw321json)

if args.derive_all or args.derive_flow:
    deriveOptFlowFusion(flowjson)

if args.derive_all or args.derive_declination:
    deriveDeclinationFusion(declinationjson)

if args.codegen:
    c_header = os.path.join(resultsdir, 'ekf_defines.h')
    generateCode(predictionjson, airspeedjson, c_header)
