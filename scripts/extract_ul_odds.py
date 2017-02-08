#!/usr/bin/env python

"""
Short script to extract the evidence ratio and a 95% H0 upper limit from a lalapps_pulsar_parameter_estimation_nested
produced posterior file. This takes one argument, which is the posterior HDF5 file.
"""

from __future__ import print_function

import os
import sys

from lalapps.pulsarpputils import pulsar_nest_to_posterior as pn2p
from lalapps.pulsarpputils import upper_limit_greedy as ulg

infile = sys.argv[1]

if not os.path.isfile(infile):
  print("Error... input file '%s' does not exist" % infile, file=sys.stderr)
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

# create JSON output dictionary
outdic = {}
outdic['upper limit'] = ul
outdic['odds'] = B

# output to text file named after input posterior
try:
  of = os.path.splitext(infile)[0]+'.txt'
  fp = open(of, 'w')
  fp.write("%.12e\t%.6e\n" % (B, ul))
  fp.close()
except:
  print("Error... could not output upper limit and odds to file.", file=sys.stderr)
  sys.exit(1)

sys.exit(0)