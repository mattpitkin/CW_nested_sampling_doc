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

# import some pulsar utilities
from lalapps import pulsarpputils as pppu

parser = argparse.ArgumentParser( )
parser.add_argument("-i", "--niter", dest="niter", default=100, type=int, help="Set the number of runs for each maximum prior value [default: %(default)s]")
parser.add_argument("-N", "--Nlive", dest="nlive", default=1024, type=int, help="Set the number of live points to use [default: %(default)s]")
parser.add_argument("-r", "--rundir", dest="rundir", required=True, help="Set the run directory for outputs")
parser.add_argument("-c", "--extract-code", dest="extract", required=True, help="Set the odds and UL extraction script location")
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
  maxvals = [1e-22, 1e-21, 1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13]
else:
  maxvals = [float(v) for v in opts.maxval]

# get code path
cpath = opts.execpath
if not os.path.isdir(cpath):
  print("Error... path to executable files '%s' does not exist." % cpath, file=sys.stderr)
  sys.exit(1)

ppen = os.path.join(cpath, 'lalapps_pulsar_parameter_estimation_nested')
if not os.path.isfile(ppen):
  print("Error... executable file '%s' does not exist." % ppen, file=sys.stderr)
  sys.exit(1)

n2p = os.path.join(cpath, 'lalapps_nest2pos')
if not os.path.isfile(ppen):
  print("Error... executable file '%s' does not exist." % n2p, file=sys.stderr)
  sys.exit(1)

# set the numbers of live points to run with
nlive = opts.nlive
# set the number of runs with each case
Nruns = opts.niter

# create log directory if it doesn't exist
logdir = os.path.join(basedir, 'log')
if not os.path.isdir(logdir):
  os.mkdir(logdir)

# setup Condor sub file for runs
subfile = os.path.join(basedir, 'ppen.sub')
fp = open(subfile, 'w')
subfiletxt = """
universe = vanilla
executable = %s
arguments = " --prior-file $(macroprior) --Nmcmcinitial 0 --outfile $(macrooutfile) --Nlive %d --test-gaussian-likelihood --test-gaussian-mean %s --test-gaussian-sigma %s --uniformprop %d --ensembleWalk %d --ensembleStretch %d"
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.o1.cw.targeted.bayesian
queue 1
""" % (ppen, nlive, opts.mean, opts.sigma, opts.uniformprop, opts.walkprop, opts.stretchprop, os.path.join(logdir, 'ppen-$(cluster).log'), \
       os.path.join(logdir,'ppen-$(cluster).err'), os.path.join(logdir,'ppen-$(cluster).out'))
fp.write(subfiletxt)
fp.close()

# setup sub file for lalapps_nest2pos jobs
n2psubfile = os.path.join(basedir, 'n2p.sub')
fp = open(n2psubfile, 'w')
subfiletxt = """
universe = vanilla
executable = %s
arguments = "  -p $(macropost) $(macronest) "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.o1.cw.targeted.bayesian
queue 1
""" % (n2p, os.path.join(logdir, 'n2p-$(cluster).log'), \
       os.path.join(logdir,'n2p-$(cluster).err'), os.path.join(logdir,'n2p-$(cluster).out'))
fp.write(subfiletxt)
fp.close()

# setup sub file for extraction script
if not os.path.isfile(opts.extract):
  print("Error... extraction script '%s' does not exist." % opts.extract, file=sys.stderr)
  sys.exit(1)

extractsubfile = os.path.join(basedir, 'extract.sub')
fp = open(extractsubfile, 'w')
subfiletxt = """
universe = vanilla
executable = %s
arguments = " $(macropost) "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.o1.cw.targeted.bayesian
queue 1
""" % (opts.extract, os.path.join(logdir, 'extract-$(cluster).log'), \
       os.path.join(logdir,'extract-$(cluster).err'), os.path.join(logdir,'extract-$(cluster).out'))
fp.write(subfiletxt)
fp.close()

# sub file for removing files
rmsubfile = os.path.join(basedir, 'rm.sub')
fp = open(rmsubfile, 'w')
subfiletxt = """
universe = local
executable = /bin/rm
arguments = " $(macrofile) "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.o1.cw.targeted.bayesian
queue 1
""" % (os.path.join(logdir, 'rm-$(cluster).log'), os.path.join(logdir,'rm-$(cluster).err'), os.path.join(logdir,'rm-$(cluster).out'))
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

    # write out n2p job
    uin2p = uuid.uuid4().hex
    outfilepost = os.path.join(maxvdir, 'post_%04d.hdf' % i)
    dagstr = 'JOB %s %s\nRETRY %s 0\nVARS %s macronest=\"%s\" macropost=\"%s\"\n' % (uin2p, n2psubfile, uin2p, uin2p, outfilenest, outfilepost)
    fp.write(dagstr)

    # write out job for extracting values
    uiextract = uuid.uuid4().hex
    dagstr = 'JOB %s %s\nRETRY %s 0\nVARS %s macropost=\"%s\"\n' % (uiextract, extractsubfile, uiextract, uiextract, outfilepost)
    fp.write(dagstr)

    # write out job for removing nested samples files and posterior files
    uirm = uuid.uuid4().hex
    dagstr = 'JOB %s %s\nRETRY %s 0\nVARS %s macropost=\"%s %s\"\n' % (uirm, rmsubfile, uirm, uirm, outfilepost, outfilenest)

    # output parents and children
    dagstr = 'PARENT %s CHILD %s\n' % (uippen, uin2p)
    fp.write(dagstr)

    dagstr = 'PARENT %s CHILD %s\n' % (uin2p, uiextract)
    fp.write(dagstr)
    
    dagstr = 'PARENT %s CHILD %s\n' % (uiextract, uirm)
    fp.write(dagstr)

fp.close()
