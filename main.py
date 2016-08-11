#!/usr/bin/python
from derivations import *
import ccodegen
import pycodegen
import os
import sys
import signal
import argparse
from multiprocessing import Process
import datetime

workers = []
def sigint_handler(sig, frame):
    for p in workers:
        if p.is_alive():
            p.terminate()
    print 'Terminated by Ctrl-C'
    exit(0)

signal.signal(signal.SIGINT,sigint_handler)

parser = argparse.ArgumentParser()
parser.add_argument('--outputdir', nargs=1, dest='outdir', default=['output'])
parser.add_argument('--derive', nargs='*', dest='derive', default=[])
parser.add_argument('--codegen', dest='codegen', action='store_true')
args = parser.parse_args()

outdir = args.outdir[0]

if not os.path.exists(outdir):
    os.makedirs(outdir)

initializationjson = os.path.join(outdir, 'initialization.json')
predictionjson = os.path.join(outdir, 'prediction.json')
camerajson = os.path.join(outdir, 'cameraFusion.json')
velNEjson = os.path.join(outdir, 'velNEFusion.json')
velDjson = os.path.join(outdir, 'velDFusion.json')
heightjson = os.path.join(outdir, 'heightFusion.json')
conversionsjson = os.path.join(outdir, 'conversions.json')

c_header = os.path.join(outdir, 'ekf_defines.h')
pythonfile = os.path.join(outdir, 'ekf.py')

derivations = {
    'initialization': (deriveInitialization, initializationjson),
    'prediction': (derivePrediction, predictionjson),
    'camera': (deriveCameraFusion, camerajson),
    'velNE': (deriveVelNEFusion, velNEjson),
    'velD': (deriveVelDFusion, velDjson),
    'height': (deriveHeightFusion, heightjson),
    'conversions': (deriveConversions, conversionsjson),
    }

assert set(args.derive).issubset(set(derivations.keys()+['all']))

if 'all' in args.derive:
    args.derive = derivations.keys()

for k in args.derive:
    workers.append(Process(target=derivations[k][0], args=(derivations[k][1],)))
    workers[-1].start()

for p in workers:
    p.join()

if args.codegen:
    jsondict = {}
    for key, val in derivations.iteritems():
        jsondict[key] = val[1]

    ccodegen.generateCode(jsondict,c_header)
    pycodegen.generateCode(jsondict,pythonfile)
