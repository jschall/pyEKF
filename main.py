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

jsonfile = os.path.join(resultsdir, 'covariancePrediction.json')
cfile = os.path.join(resultsdir, 'covariancePrediction.c')

if args.derive:
    deriveCovariancePrediction(jsonfile)
    print('Covariance predicton derivation saved to %s' % (jsonfile,))

if args.codegen:
    generateCovariancePrediction(jsonfile, cfile)
    print('Covariance predicton c code saved to %s' % (cfile,))
