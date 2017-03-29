#!/usr/bin/env python

"""
The main purpose of this script is to assess the distribution of evidence ratios and
upper limits produced by lalapps_pulsar_parameter_estimation_nested when running on the
Gaussian test likelihood. It will run with a uniform prior distribution bounded at 0,
with a variety of increasing maximum values, to see how the evidence calculation
performs. It can run with various proposal distributions.

This is all done by setting up a Condor DAG to run the required codes.

Copyright (C) 2017 Matthew Pitkin
"""

from __future__ import print_function

import os
import sys
import uuid
import argparse

parser = argparse.ArgumentParser( )
parser.add_argument("-i", "--niter", dest="niter", default=100, type=int, help="Set the number of runs for each maximum prior value [default: %(default)s]")
parser.add_argument("-N", "--Nlive", dest="nlive", default=1024, type=int, help="Set the number of live points to use [default: %(default)s]")
parser.add_argument("-r", "--rundir", dest="rundir", required=True, help="Set the run directory for outputs")
parser.add_argument("-c", "--run-code", dest="run", required=True, help="Set the test run script location")
parser.add_argument("-p", "--exec-path", dest="execpath", required=True, help="Set the path to the required executables")
parser.add_argument("-m", "--maxval", dest="maxval", action='append', default=None, help="A maximum limit for the prior (this can be specified multiple times to use more than one value).")
parser.add_argument("-u", "--uniformprop", dest="uniformprop", type=int, default=0, help="Set the amount of time to use the uniform proposal [default: %(default)s]")
parser.add_argument("-w", "--walkprop", dest="walkprop", type=int, default=1, help="Set the amount of time to use the ensemble walk proposal [default: %(default)s]")
parser.add_argument("-s", "--stretchprop", dest="stretchprop", type=int, default=0, help="Set the amount of time to use the ensemble stretch proposal [default: %(default)s]")
parser.add_argument("--gauss-mean", dest="mean", default="0.0", help="Set the test Gaussian likelihood mean value [default: %(default)s]")
parser.add_argument("--gauss-sigma", dest="sigma", default="1e-24", help="Set the test Gaussian likelihood standard deviation value [default: %(default)s]")

# parse input options
opts = parser.parse_args()

# the base directory
basedir = opts.rundir
if not os.path.isdir(basedir):
  print("Error... base directory '%s' does not exist." % basedir, file=sys.stderr)
  sys.exit(1)

if opts.maxval is None:
  # set the maximum values of the prior distribution
  maxvals = [1e-23, 1e-22, 1e-21, 1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13]
else:
  maxvals = [float(v) for v in opts.maxval]

# set the numbers of live points to run with
nlive = opts.nlive
# set the number of runs with each case
Nruns = opts.niter

# create log directory if it doesn't exist
logdir = os.path.join(basedir, 'log')
if not os.path.isdir(logdir):
  os.mkdir(logdir)

# setup sub file for extraction script
if not os.path.isfile(opts.run) or not os.access(opts.run, os.X_OK):
  print("Error... test run script '%s' does not exist, or is not executable." % opts.run, file=sys.stderr)
  sys.exit(1)

# check executable path is a directory
if not os.path.isdir(opts.execpath):
  print("Error... path for run executables does not exist.", file=sys.stderr)
  sys.exit(1)

# setup Condor sub file for runs
subfile = os.path.join(basedir, 'runtestgauss.sub')
fp = open(subfile, 'w')
subfiletxt = """universe = vanilla
executable = %s
arguments = " --outfile $(macrooutfile) --Nlive %d --gauss-mean %s --gauss-sigma %s --uniformprop %d --walkprop %d --stretchprop %d --exec-path %s --priorfile $(macroprior) "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = aluk.dev.o1.cw.targeted.bayesian
queue 1
""" % (opts.run, nlive, opts.mean, opts.sigma, opts.uniformprop, opts.walkprop, opts.stretchprop, opts.execpath, os.path.join(logdir, 'run-$(cluster).log'), \
       os.path.join(logdir,'run-$(cluster).err'), os.path.join(logdir,'run-$(cluster).out'))
fp.write(subfiletxt)
fp.close()

# create dag for all the jobs
dagfile = os.path.join(basedir, 'run.dag')
fp = open(dagfile, 'w')

for i, maxv in enumerate(maxvals):
  # create output directory
  maxvdir = os.path.join(basedir, '%03d' % i)
  if not os.path.isdir(maxvdir):
    os.mkdir(maxvdir)

  # create prior file
  priorfile = os.path.join(maxvdir, 'prior.txt')
  fpm = open(priorfile, 'w')
  fpm.write('H0 uniform 0.0 %.4e\n' % maxv)
  fpm.close()

  # loop over number of iterations
  for j in range(Nruns):
    # unique ID
    uippen = uuid.uuid4().hex

    # create output nested file
    outfilenest = os.path.join(maxvdir, 'nest_%04d.hdf' % j)
    
    # write out ppen job
    dagstr = 'JOB %s %s\nRETRY %s 0\nVARS %s macrooutfile=\"%s\" macroprior=\"%s\"\n' % (uippen, subfile, uippen, uippen, outfilenest, priorfile)
    fp.write(dagstr)

fp.close()
