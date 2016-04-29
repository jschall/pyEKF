from derivations import *
from codegen import *
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--derive-all', dest='derive_all', action='store_true')
parser.add_argument('--derive-prediction', dest='derive_prediction', action='store_true')
parser.add_argument('--derive-rotation', dest='derive_rotation', action='store_true')
parser.add_argument('--derive-airspeed', dest='derive_airspeed', action='store_true')
parser.add_argument('--derive-beta', dest='derive_beta', action='store_true')
parser.add_argument('--derive-mag', dest='derive_mag', action='store_true')
parser.add_argument('--codegen', dest='codegen', action='store_true')
args = parser.parse_args()

resultsdir = 'results'

if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)


if args.derive_all or args.derive_prediction:
    jsonfile = os.path.join(resultsdir, 'covariancePrediction.json')
    deriveCovariancePrediction(jsonfile)
    print('Covariance predicton derivation saved to %s' % (jsonfile,))

if args.derive_all or args.derive_rotation:
    jsonfile = os.path.join(resultsdir, 'covarianceRotation.json')
    deriveCovarianceRotation(jsonfile)
    print('Covariance rotation derivation saved to %s' % (jsonfile,))

if args.derive_all or args.derive_airspeed:
    jsonfile = os.path.join(resultsdir, 'airspeedFusion.json')
    deriveAirspeedFusion(jsonfile)
    print('Airspeed fusion derivation saved to %s' % (jsonfile,))

if args.derive_all or args.derive_beta:
    jsonfile = os.path.join(resultsdir, 'betaFusion.json')
    deriveBetaFusion(jsonfile)
    print('Beta fusion derivation saved to %s' % (jsonfile,))

if args.derive_all or args.derive_mag:
    jsonfile = os.path.join(resultsdir, 'magFusion.json')
    deriveMagFusion(jsonfile)
    print('Mag fusion derivation saved to %s' % (jsonfile,))

if args.codegen:
    try:
        jsonfile = os.path.join(resultsdir, 'covariancePrediction.json')
        cfile = os.path.join(resultsdir, 'covariancePrediction.c')
        generateCovariance(jsonfile, cfile)
        print('Covariance predicton c code saved to %s' % (cfile,))
    except:
        pass

    try:
        jsonfile = os.path.join(resultsdir, 'covarianceRotation.json')
        cfile = os.path.join(resultsdir, 'covarianceRotation.c')
        generateCovariance(jsonfile, cfile)
        print('Covariance rotation c code saved to %s' % (cfile,))
    except:
        pass
