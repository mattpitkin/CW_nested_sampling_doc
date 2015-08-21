#!/usr/bin/env python

"""
The main purpose of this script is to assess the distribution of evidence ratios produced by 
lalapps_pulsar_parameter_estimation_nested when running on the same piece of data for different numbers of 
live points. It will also check the consistency of the evidence ratio with that produced through a grid-based 
trapezium rule numerical integration for the same data. It will finally look at the value of the 95% upper 
limit on signal amplitude and how that varies with the number of live points used.

It will perform these tests both the Student's t and Gaussian likelihood functions.

Copyright (C) 2015 Matthew Pitkin
"""

import os
import subprocess
import numpy as np
import uuid
import json

# import some pulsar utilities
from lalapps import pulsarpputils as pppu

from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz

# the base directory
basedir = '/home/sismp2/projects/code_testing/evidence_ul_distribution'

# create a set of some fake noise
mu = 0. # mean
sigma = 1e-23 # standard deviation

fakedata = np.zeros((1440, 3))
fakedata[:,0] = np.arange(900000000., 900086400., 60) # GPS time stamps
fakedata[:,1:] = mu + sigma*np.random.randn(1440, 2)  # real and imaginary data

# output the data
datafile = os.path.join(basedir, 'data.txt.gz')
np.savetxt(datafile, fakedata, fmt='%.1f\t%.7e\t%.7e')

# include the standard deviation in the file
fakedatasigma = np.zeros((1440, 4))
fakedatasigma[:,:3] = fakedata
fakedatasigma[:,-1] = sigma*np.ones(1440)

datafilesigma = os.path.join(basedir, 'datasigma.txt.gz')
np.savetxt(datafilesigma, fakedatasigma, fmt='%.1f\t%.7e\t%.7e\t%.7e')

# create a prior file
priorfile = os.path.join(basedir, 'prior.txt')
fp = open(priorfile, 'w')
h0max = 1e-20 # maximum of h0 range
priortxt = """
H0 uniform 0 %.2e
COSIOTA uniform -1 1
PHI0 uniform 0 %.8f
PSI uniform 0 %.8f
""" % (h0max, np.pi, np.pi/2.)
fp.write(priortxt)
fp.close()

# create a par file
parfile = os.path.join(basedir, 'pulsar.par')
fp = open(parfile, 'w')
partxt = """
PSRJ J0000-0000
RAJ 00:00:00.0
DECJ 00:00:00.0
PEPOCH 55000.0
F0 100.0
"""
fp.write(partxt)
fp.close()

# set the numbers of live points to run with
nlives = [128, 256, 512, 1024, 2048, 4096]

# set the number of runs with each case
Nruns = 125

# create log directory if it doesn't exist
logdir = os.path.join(basedir, 'log')
if not os.path.isdir(logdir):
  os.mkdir(logdir)

# setup Condor sub file for runs
subfile = os.path.join(basedir, 'run.sub')
execu = os.path.join(os.environ['LSCSOFT_LOCATION'], 'bin/lalapps_pulsar_parameter_estimation_nested')
fp = open(subfile, 'w')
subfiletxt = """
universe = vanilla
executable = %s
arguments = " --prior-file %s --detectors H1 --par-file %s --Nmcmcinitial 200 --outfile $(macrooutfile) \
--Nlive $(macronlive) --gzip --input-files $(macroinput) $(macrolike) --oldChunks "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.s6.cw.targeted.bayesian
queue 1
""" % (execu, priorfile, parfile, os.path.join(logdir, '$(cluster).log'), \
       os.path.join(logdir,'$(cluster).err'), os.path.join(logdir,'$(cluster).out'))
fp.write(subfiletxt)
fp.close()

# setup sub file for lalapps_nest2pos jobs
n2psubfile = os.path.join(basedir, 'n2p.sub')
n2pexec = os.path.join(os.environ['LSCSOFT_LOCATION'], 'bin/lalapps_nest2pos')
fp = open(n2psubfile, 'w')
subfiletxt = """
universe = vanilla
executable = %s
arguments = " -N $(macronlive2) -p $(macropost) -H $(macroheader) -z $(macronest) "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.s6.cw.targeted.bayesian
queue 1
""" % (n2pexec, os.path.join(logdir, '$(cluster).log'), \
       os.path.join(logdir,'$(cluster).err'), os.path.join(logdir,'$(cluster).out'))
fp.write(subfiletxt)
fp.close()

# create dag for all the jobs
dagfile = os.path.join(basedir, 'run.dag')
fp = open(dagfile, 'w')

# perform analysis for both Student's t and Gaussian likelihoods
likelihoods = ['studentst', 'gaussian']

for n in nlives:
  # create output directory
  livedir = os.path.join(basedir, '%d' % n)
  if not os.path.isdir(livedir):
    os.mkdir(livedir)
  
  for l in likelihoods:
    likedir = os.path.join(livedir, l)
 
    if not os.path.isdir(likedir):
      os.mkdir(likedir)
 
    for i in range(Nruns):
      # create output file
      outfile = os.path.join(likedir, 'nest_%04d.txt' % i)
    
      lval = ''
      infile = datafile
      if l == 'gaussian':
        lval = '--gaussian-like'
        infile = datafilesigma
    
      # unique ID
      ui = uuid.uuid4().hex
      dagstr = 'JOB %s %s\nRETRY %s 0\nVARS %s macrooutfile=\"%s\" macronlive=\"%d\" macroinput=\"%s\" macrolike=\"%s\"\n' % (ui, subfile, ui, ui, outfile, n, infile, lval)
      fp.write(dagstr)

      # add nest2pos
      ui2 = uuid.uuid4().hex
      postfile = os.path.join(likedir, 'post_%04d.txt' % i)
      dagstr = 'JOB %s %s\nRETRY %s 0\nVARS %s macronlive2=\"%d\" macronest=\"%s\" macropost=\"%s\" macroheader=\"%s\"\n' % (ui2, n2psubfile, ui2, ui2, n, outfile+'.gz', postfile, outfile+'_params.txt')
      fp.write(dagstr)

      dagstr = 'PARENT %s CHILD %s\n' % (ui, ui2)
      fp.write(dagstr)

fp.close()

# run the grid-based posterior function to get UL and evidence ratio for comparison

# setup grid
paramranges = {}
paramranges['h0'] = (0., 2.*sigma, 100)
paramranges['psi'] = (0., np.pi/2., 50)
paramranges['phi0'] = (0., np.pi, 50)
paramranges['cosiota'] = (-1., 1., 50)
ra = 0.0
dec = 0.0
dets = 'H1'
ts = {}
ts[dets] = fakedata[:, 0]
data = {}
data[dets] = fakedata[:,1] + 1j*fakedata[:,2]

outdict = {}
outdict['Odds ratios'] = {}
outdict['Upper limits'] = {}

for l in likelihoods:
  sigmas = None
  if l == 'gaussian':
    sigmas = {}
    sigmas[dets] = fakedatasigma[:,3]


  L, h0pdf, phi0pdf, psipdf, cosiotapdf, grid, evrat = pppu.pulsar_posterior_grid(dets, ts, data, ra, dec, \
                                                                                  sigmas=sigmas, \
                                                                                  paramranges=paramranges)

  # scale evidence ratio
  evrat += np.log(2.*sigma/h0max)

  # get cumulative distribution of h0
  ct = cumtrapz(h0pdf, grid['h0'])

  # get 95% amplitude upper limits
  ctu, ui = np.unique(ct, return_index=True)
  intf = interp1d(ctu, grid['h0'][ui], kind='linear')
  h95ul = intf(0.95)
  
  outdict['Odds ratios'][l] = evrat
  outdict['Upper limits'][l] = float(h95ul)

infofile = os.path.join(basedir, 'gridoutput.txt')
fp = open(infofile, 'w')
json.dump(outdict, fp, indent=2)
fp.close()
