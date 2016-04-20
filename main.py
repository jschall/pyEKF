from derivations import *
from codegen import *
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--no-derive', dest='derive', action='store_false')
parser.add_argument('--no-codegen', dest='codegen', action='store_false')
args = parser.parse_args()

resultsdir = 'results'

if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)




if args.derive:
    jsonfile = os.path.join(resultsdir, 'covariancePrediction.json')
    deriveCovariancePrediction(jsonfile)
    print('Covariance predicton derivation saved to %s' % (jsonfile,))

    jsonfile = os.path.join(resultsdir, 'airspeedFusion.json')
    deriveAirspeedFusion(jsonfile)
    print('Airspeed fusion derivation saved to %s' % (jsonfile,))

    jsonfile = os.path.join(resultsdir, 'betaFusion.json')
    deriveBetaFusion(jsonfile)
    print('Beta fusion derivation saved to %s' % (jsonfile,))

if args.codegen:
    jsonfile = os.path.join(resultsdir, 'covariancePrediction.json')
    cfile = os.path.join(resultsdir, 'covariancePrediction.c')
    generateCovariancePrediction(jsonfile, cfile)
    print('Covariance predicton c code saved to %s' % (cfile,))
