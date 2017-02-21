#!/usr/bin/env python

"""
Run lalapps_pulsar_parameter_estimation_nested analysis over a range of SNRs for signals with parameters drawn randomly from
some set prior ranges for the four main GW parameters: h0, phi0, cosiota and psi. Sources will be distributed randomly on the sky.
NOTE: in reality prior if parameters are drawn from a flat prior, and the recovery prior is also flat, but much broader than
the drawing distributionm then p-p plot should still work (problems will occur is the prior is not flat).

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

# range of (multi-detector) SNRs to use
snrrange = [0., 20.]
# set flat prior ranges on other parameters
h0range = [0., 1e-20]
phi0range = [0., np.pi]
psirange = [0., np.pi/2.]
cirange = [-1., 1.]

# set fixed h0 value (that will be scale to SNR)
h0fixed = 1e-24

# run base directory (on RAVEN)
basedir = '/home/sismp2/projects/testing/pp_standard'
#basedir = '/home/matthew/testing/lalapps_knope_O2/outdir'

logdir = os.path.join(basedir, 'log')
if not os.path.isdir(logdir):
  os.makedirs(logdir)

# detectors to use
dets = ['H1', 'L1']

nsigs = 2000 # total number of signals

# generate signal parameters
snrs = np.random.rand(nsigs)*np.diff(snrrange)[0] # randomly distributed SNRs
# generate RA and dec values uniformly over the sky
ras = 2.*np.pi*np.random.rand(nsigs)
decs = -(np.pi/2.) + np.arccos(2.*np.random.rand(nsigs) - 1.)

# set injection parameters
h0s = np.ones(nsigs)*h0fixed
phi0s = phi0range[0] + np.diff(phi0range)[0]*np.random.rand(nsigs)  # phi0 value
psis = psirange[0] + np.diff(psirange)[0]*np.random.rand(nsigs)     # psi value
cis = cirange[0] + np.diff(cirange)[0]*np.random.rand(nsigs)        # cos(iota) value

nlive = 1024 # number of live points
nruns = 2    # number of parallel runs for each analysis

# create sub file for running run_lppen.py
subfile = os.path.join(basedir, 'lppen.sub')
subdata = """universe = vanilla
executable = /home/sismp2/repositories/CW_nested_sampling_doc/scripts/run_lppen.py
arguments = " $(macroinifile) "
getenv = True
log = %s
error = %s
output = %s
notification = never
accounting_group = ligo.dev.o1.cw.targeted.bayesian
queue 1
""" % (os.path.join(logdir, 'run-$(cluster).log'), os.path.join(logdir,'run-$(cluster).err'), os.path.join(logdir,'run-$(cluster).out'))

fp = open(subfile, 'w')
fp.write(subdata)
fp.close()

# create dag file
dagfile = os.path.join(basedir, 'lppen.dag')
fp = open(dagfile, 'w')

priorstr = """'H0 uniform {} {}\\nPHI0 uniform {} {}\\nPSI uniform {} {}\\nCOSIOTA uniform {} {}'
""".format(h0range[0], h0range[1], phi0range[0], phi0range[1], psirange[0], psirange[1], cirange[0], cirange[1])

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

  # set RA and DEC
  h, m, s = pppu.rad_to_hms(ras[i])
  rastr = pppu.coord_to_string(h, m, s)
  d, m, s = pppu.rad_to_dms(decs[i])
  decstr = pppu.coord_to_string(d, m, s)

  cprun.set('run', 'hetparams', json.dumps({'RAJ': rastr, 'DECJ': decstr}))

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

  cprun.add_section('injection')
  cprun.set('injection', 'snr', str(snrs[i]))
  cprun.set('injection', 'injparams', json.dumps([{'RAJ': rastr, 'DECJ': decstr, 'H0': h0s[i], 'PHI0': phi0s[i], 'PSI': psis[i], 'COSIOTA': cis[i]}]))
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
