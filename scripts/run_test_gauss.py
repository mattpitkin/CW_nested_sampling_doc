#!/usr/bin/env python

"""
Short script to extract the evidence ratio and a 95% H0 upper limit from a lalapps_pulsar_parameter_estimation_nested
produced posterior file. This takes one argument, which is the posterior HDF5 file.
"""

from __future__ import print_function, division

import os
import sys
import h5py
import numpy as np
import argparse
from scipy import stats
from scipy.special import erf
import subprocess as sp

from lalapps.pulsarpputils import pulsar_nest_to_posterior as pn2p
from lalapps.pulsarpputils import upper_limit_greedy as ulg

# define the CDF for a truncated Gaussian
def truncgauss_cdf(x, mu, sig, minv, maxv):
  Cv = 1./(np.sqrt(2.)*sig)
  areav = 0.5*(erf(Cv*(mu - minv)) - erf(Cv*(mu - maxv)))
  return 0.5*(erf(Cv*(mu-minv)) - erf(Cv*(mu-x)))/areav

parser = argparse.ArgumentParser( )
parser.add_argument("--outfile", dest="outfile", required=True, help="The full path to the output nested sampling HDF5 file") # the positional argument for the posterior file
parser.add_argument("--gauss-mean", dest="mean", type=float, required=True, help="Set the test Gaussian likelihood mean value")
parser.add_argument("--gauss-sigma", dest="sigma", type=float, required=True, help="Set the test Gaussian likelihood standard deviation value")
parser.add_argument("--min-val", dest="minv", type=float, default=None, help="Set the minumum of the prior range")
parser.add_argument("--max-val", dest="maxv", type=float, default=None, help="Set the maximum of the prior range")
parser.add_argument("-N", "--Nlive", dest="nlive", default=1024, type=int, help="Set the number of live points to use [default: %(default)s]")
parser.add_argument("-p", "--exec-path", dest="execpath", required=True, help="Set the path to the required executables")
parser.add_argument("-u", "--uniformprop", dest="uniformprop", type=int, default=0, help="Set the amount of time to use the uniform proposal [default: %(default)s]")
parser.add_argument("-w", "--walkprop", dest="walkprop", type=int, default=1, help="Set the amount of time to use the ensemble walk proposal [default: %(default)s]")
parser.add_argument("-s", "--stretchprop", dest="stretchprop", type=int, default=0, help="Set the amount of time to use the ensemble stretch proposal [default: %(default)s]")
parser.add_argument("-x", "--priorfile", dest="priorfile", default=None, help="Set a prior file to use, if not specifying max and min ranges")
parser.add_argument("-r", "--dont-remove", dest="drm", action='store_true', default=False, help="Set to NOT remove files (will be removed by default)")

# parse input options
opts = parser.parse_args()

outfile = opts.outfile

# get base directory from outfile name
basedir = os.path.dirname(outfile)

# check minimum and maximum prior ranges are valid, or if a prior file has been passed
priorfile = None
minv = 0.
maxv = np.inf
if opts.minv is not None and opts.maxv is not None:
  minv = opts.minv
  maxv = opts.maxv
elif opts.priorfile is not None:
  priorfile = opts.priorfile
  if not os.path.isfile(priorfile):
    print("Error... prior file '%s' does not exist." % priorfile, file=sys.stderr)
    sys.exit(1)
  try:
    fp = open(priorfile, 'r')
    l = fp.readline().split()
    fp.close()
    minv = float(l[-2].strip())
    maxv = float(l[-1].strip())
  except:
    print("Error... could not read minimum and maximum range from prior file.", file=sys.stderr)
    sys.exit(1)
else:
  print("Error... must either specify a max and min prior range, or give a prior file.", file=sys.stderr)
  sys.exit(1)

if minv > maxv:
  print("Error... minimum prior ranges is greater than maximum.", file=sys.stderr)
  sys.exit(1)

# get code path
cpath = opts.execpath
if not os.path.isdir(cpath):
  print("Error... path to executable files '%s' does not exist." % cpath, file=sys.stderr)
  sys.exit(1)

ppen = os.path.join(cpath, 'lalapps_pulsar_parameter_estimation_nested')
if not os.path.isfile(ppen) or not os.access(ppen, os.X_OK):
  print("Error... executable file '%s' does not exist or is not executable." % ppen, file=sys.stderr)
  sys.exit(1)

n2p = os.path.join(cpath, 'lalapps_nest2pos')
if not os.path.isfile(n2p) or not os.access(n2p, os.X_OK):
  print("Error... executable file '%s' does not exist or is not executable." % n2p, file=sys.stderr)
  sys.exit(1)

# set the prior file
if priorfile is None:
  priorfile = os.path.splitext(outfile)[0]+'.prior'
  fp = open(priorfile, 'w')
  fp.write("H0 uniform {} {}\n".format(minv, maxv))
  fp.close()

# check proposal ratios
uniformprop = 0
walkprop = 0
stretchprop = 0
if opts.uniformprop > 0: # only set values if greater than zero
  uniformprop = opts.uniformprop
if opts.walkprop > 0:
  walkprop = opts.walkprop
if opts.stretchprop > 0:
  stretchprop = opts.stretchprop

if not uniformprop and not walkprop and not stretchprop:
  print("Error... problem with given proposals.", file=sys.stderr)
  sys.exit(1)

if opts.nlive < 1:
  print("Error... problem with given number of live points.", file=sys.stderr)
  sys.exit(1)

# run lalapps_pulsar_parameter_estimation_nested
ppencodecall = "{} --prior-file {} --uniformprop {} --ensembleWalk {} --ensembleStretch {} --test-gaussian-likelihood --test-gaussian-mean {} --test-gaussian-sigma {} --outfile {} --Nlive {} --Nmcmcinitial 0"
p = sp.Popen(ppencodecall.format(ppen, priorfile, str(uniformprop), str(walkprop), str(stretchprop), str(opts.mean), str(opts.sigma), outfile, str(opts.nlive)), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
out, err = p.communicate()
print(err, file=sys.stderr)

# run lalapps_nest2pos
postfile = os.path.splitext(outfile)[0]+'_post.hdf'
n2pcodecall = "{} -p {} {}"
p = sp.Popen(n2pcodecall.format(n2p, postfile, outfile), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
out, err = p.communicate()
print(err, file=sys.stderr)

try:
  post, B, N = pn2p(postfile)
except:
  print("Error... could not extract samples from posterior file '%s'" % postfile, file=sys.stderr)
  sys.exit(1)

# get 95% upper limit for H0 samples
try:
  ul = ulg(post['H0'].samples)
except:
  print("Error... could not get 'H0' upper limit from posterior samples", file=sys.stderr)
  sys.exit(1)

# get the information and work out the expected error in the logZ value
hdf = h5py.File(postfile, 'r')
a = hdf['lalinference']['lalinference_nest']
info = a.attrs['information_nats']
nlive = a.attrs['number_live_points']

onesigerror = np.sqrt(info/nlive)

# get the two-sided KS statistic p-value
try:
  D, ksp = stats.kstest(post['H0'].samples[:,0], truncgauss_cdf, args=(opts.mean, opts.sigma, minv, maxv))
except:
  print("Warning... could not perform KS test.", file=sys.stderr)
  ksp = 0.

# output to text file named after input posterior
try:
  of = os.path.splitext(outfile)[0]+'_stats.txt'
  fp = open(of, 'w')
  fp.write("%.12e\t%.6e\t%.6e\t%.6e\n" % (B, ul, onesigerror, ksp))
  fp.close()
except:
  print("Error... could not output upper limit and odds to file.", file=sys.stderr)
  sys.exit(1)

# clean up posterior file, nested sampling file and prior file
if not opts.drm:
  try:
    if opts.priorfile is None:
      os.remove(priorfile) # remove prior file, unless is has been passed to the code
    os.remove(outfile)   # remove nested sample file
    os.remove(postfile)  # remove posterior sample file
  except:
    print("Warning... could not remove various files.", file=sys.stderr)

sys.exit(0)
