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

from lalapps.pulsarpputils import pulsar_nest_to_posterior as pn2p
from lalapps.pulsarpputils import upper_limit_greedy as ulg

# define the CDF for a truncated Gaussian
def truncgauss_cdf(x, mu, sig, minv, maxv):
  Cv = 1./(np.sqrt(2.)*sig)
  areav = 0.5*(erf(Cv*(mu - minv)) - erf(Cv*(mu - maxv)))
  return 0.5*(erf(Cv*(mu-minv)) - erf(Cv*(mu-x)))/areav

parser = argparse.ArgumentParser( )
parser.add_argument("infile", help="The input posterior HDF5 file") # the positional argument for the posterior file
parser.add_argument("--gauss-mean", dest="mean", type=float, required=True, help="Set the test Gaussian likelihood mean value")
parser.add_argument("--gauss-sigma", dest="sigma", type=float, required=True, help="Set the test Gaussian likelihood standard deviation value")
parser.add_argument("--min-val", dest="minv", type=float, required=True, help="Set the minumum of the prior range")
parser.add_argument("--max-val", dest="maxv", type=float, required=True, help="Set the maximum of the prior range")

# parse input options
opts = parser.parse_args()

infile = opts.infile

if not os.path.isfile(infile):
  print("Error... input file '%s' does not exist" % infile, file=sys.stderr)
  sys.exit(1)

if opts.minv > opts.maxv:
  print("Error... minumum prior ranges is greater than maximum.", file=sys.stderr)
  sys.exit(1)

try:
  post, B, N = pn2p(infile)
except:
  print("Error... could not extract samples from posterior file '%s'" % infile, file=sys.stderr)
  sys.exit(1)

# get 95% upper limit for H0 samples
try:
  ul = ulg(post['H0'].samples)
except:
  print("Error... could not get 'H0' upper limit from posterior samples", file=sys.stderr)
  sys.exit(1)

# get the information and work out the expected error in the logZ value
hdf = h5py.File(infile, 'r')
a = hdf['lalinference']['lalinference_nest']
info = a.attrs['information_nats']
nlive = a.attrs['number_live_points']

onesigerror = np.sqrt(info/nlive)

# get the two-sided KS statistic p-value
D, ksp = stats.kstest(post['H0'].samples[:,0], truncgauss_cdf, args=(opts.mean, opts.sigma, opts.minv, opts.maxv))

# output to text file named after input posterior
try:
  of = os.path.splitext(infile)[0]+'.txt'
  fp = open(of, 'w')
  fp.write("%.12e\t%.6e\t%.6e\t%.6f\n" % (B, ul, onesigerror, ksp))
  fp.close()
except:
  print("Error... could not output upper limit and odds to file.", file=sys.stderr)
  sys.exit(1)

sys.exit(0)