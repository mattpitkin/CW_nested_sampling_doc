#!/usr/bin/env python

"""
Run lalapps_pulsar_parameter_estimation_nested analysis over a range of SNRs for signals with parameters drawn randomly from
some set prior ranges for the four main GW parameters (h0, phi0, cosiota and psi) and f0, f1 and RA. Sources will be distributed
with the same declination, but with RA drawn from a Gaussian over a small range. f0 and f1 will also be drawn from Gaussian
distributions over a small range.

Everything will be hardcoded in this script.

These will be run as a Condor DAG.
"""

import os
import sys
import numpy as np
import uuid
import json

from ConfigParser import ConfigParser

import lalapps.pulsarpputils as pppu

# set flat prior ranges on other parameters
h0range = [0., 1e-20]
phi0range = [0., np.pi]
psirange = [0., np.pi/2.]
cirange = [-1., 1.]

# set range of injected h0 values
h0injrange = [0., 3.25e-22]

# run base directory (on RAVEN)
basedir = '/home/sismp2/projects/testing/pp_roq_uniform'
#basedir = '/home/matthew/testing/lalapps_knope_O2/outdir'

logdir = os.path.join(basedir, 'log')
if not os.path.isdir(logdir):
  os.makedirs(logdir)

# detectors to use
dets = ['H1', 'L1']
#dets = ['H1']

nsigs = 500 # total number of signals
#nsigs = 1

# generate signal parameters
snrs = np.zeros(nsigs) # randomly distributed SNRs
# generate RA and dec values uniformly over the sky
ramean = np.pi # central right ascension
rasigma = (10./(60.*60.*24))*2.*np.pi # 10 seconds
ras = ramean + np.random.randn(nsigs)*rasigma
decs = np.zeros(nsigs) # set declination at zero degrees

h, m, s = pppu.rad_to_hms(ramean)
rameanstr = pppu.coord_to_string(h, m, s)
d, m, s = pppu.rad_to_dms(decs[0])
decstr = pppu.coord_to_string(d, m, s)

# frequency parameters
f0mean = 100.
f0sigma = 5e-5
f1mean = -1e-9
f1sigma = 2e-10
f0s = f0mean + f0sigma*np.random.randn(nsigs)
f1s = f1mean + f1sigma*np.random.randn(nsigs)

# set injection parameters
h0s = h0injrange[0] + np.diff(h0injrange)[0]*np.random.rand(nsigs)
phi0s = phi0range[0] + np.diff(phi0range)[0]*np.random.rand(nsigs)  # phi0 value
psis = psirange[0] + np.diff(psirange)[0]*np.random.rand(nsigs)     # psi value
cis = cirange[0] + np.diff(cirange)[0]*np.random.rand(nsigs)        # cos(iota) value

nlive = 2048 # number of live points
#nlive = 256
nruns = 2    # number of parallel runs for each analysis
ncpus = 2    # number of cores required

# create sub file for running run_lppen.py
subfile = os.path.join(basedir, 'lppen.sub')
subdata = """universe = parallel
executable = /home/sismp2/repositories/CW_nested_sampling_doc/scripts/run_lppen.py
arguments = " -r $(macroinifile) "
getenv = True
log = %s
error = %s
output = %s
request_cpus = %d
notification = never
accounting_group = ligo.dev.o1.cw.targeted.bayesian
queue 1
""" % (os.path.join(logdir, 'run-$(cluster).log'), os.path.join(logdir,'run-$(cluster).err'), os.path.join(logdir,'run-$(cluster).out'), ncpus)

fp = open(subfile, 'w')
fp.write(subdata)
fp.close()

# create dag file
dagfile = os.path.join(basedir, 'lppen.dag')
fp = open(dagfile, 'w')

priorstr = """'H0 uniform {} {}\\nPHI0 uniform {} {}\\nPSI uniform {} {}\\nCOSIOTA uniform {} {}\\nF0 gaussian {} {}\\nF1 gaussian {} {}\\nRA gaussian {} {}'
""".format(h0range[0], h0range[1], phi0range[0], phi0range[1], psirange[0], psirange[1], cirange[0], cirange[1], f0mean, f0sigma, f1mean, f1sigma, ramean, rasigma)

lppenexec = '/home/sismp2/lscsoft/.virtualenvs/lalapps_knope_O2/bin/lalapps_pulsar_parameter_estimation_nested'
n2pexec = '/home/sismp2/lscsoft/.virtualenvs/lalapps_knope_O2/bin/lalapps_nest2pos'

datasigma = [1e-22] # data standard deviation

for i in range(nsigs):
  cprun = ConfigParser() # set config parser to output .ini configuration file

  # set executables configuration
  cprun.add_section('executables')
  cprun.set('executables', 'lppen', lppenexec)
  cprun.set('executables', 'n2p', n2pexec)

  # set run information
  outdir = os.path.join(basedir, '%04d' % i)
  if not os.path.isdir(outdir):
    os.makedirs(outdir)
  cprun.add_section('run')
  cprun.set('run', 'outdir', outdir)
  cprun.set('run', 'outname', 'nest_%04d' % i)
  cprun.set('run', 'detectors', json.dumps(dets))
  cprun.set('run', 'nruns', str(nruns))
  cprun.set('run', 'ncpus', str(ncpus))

  # set RA
  h, m, s = pppu.rad_to_hms(ras[i])
  rastr = pppu.coord_to_string(h, m, s)

  cprun.set('run', 'hetparams', json.dumps({'RAJ': rameanstr, 'DECJ': decstr, 'F0': f0mean, 'F1': f1mean}))

  # set simulated data parameters
  cprun.add_section('data')
  cprun.set('data', 'sigma', json.dumps(datasigma))
  cprun.set('data', 'start', json.dumps([900000000]))
  cprun.set('data', 'step', json.dumps([60]))
  cprun.set('data', 'length', json.dumps([1440]))

  # set nested sampling parameters
  cprun.add_section('nestedsampling')
  cprun.set('nestedsampling', 'nlive', str(nlive))
  cprun.set('nestedsampling', 'priordata', priorstr)

  # set use of ROQ
  cprun.set('nestedsampling', 'roq', 'True')
  cprun.set('nestedsampling', 'roqntraining', '2500')
  cprun.set('nestedsampling', 'roqtolerance', '5e-12')

  cprun.add_section('injection')
  cprun.set('injection', 'snr', str(snrs[i]))
  cprun.set('injection', 'injparams', json.dumps([{'RAJ': rastr, 'DECJ': decstr, 'H0': h0s[i], 'PHI0': phi0s[i], 'PSI': psis[i], 'COSIOTA': cis[i], 'F0': f0s[i], 'F1': f1s[i]}]))
  cprun.set('injection', 'coherent', 'True')

  # write out configuration file
  runconfig = os.path.join(outdir, 'config.ini')
  fc = open(runconfig, 'w')
  cprun.write(fc)
  fc.close()

  # write out to dag file
  uippen = uuid.uuid4().hex

  # write out ppen job
  dagstr = 'JOB {} {}\nRETRY {} 0\nVARS {} macroinifile=\"{}\"\n'.format(uippen, subfile, uippen, uippen, runconfig)
  fp.write(dagstr)

# close dag file
fp.close()
