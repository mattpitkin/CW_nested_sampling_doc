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

# set of SNRs to use
snrs = [0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.]

# run base directory (on RAVEN)
basedir = '/home/sismp2/projects/testing/pp_standard'
#basedir = '/home/matthew/testing/lalapps_knope_O2/outdir'

logdir = os.path.join(basedir, 'log')
if not os.path.isdir(logdir):
  os.makedirs(logdir)

# detectors to use
dets = ['H1', 'L1']

nsigs = 100 # number of signals per SNR

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

# set flat prior ranges
h0range = [0., 1e-20]
phi0range = [0., np.pi]
psirange = [0., np.pi/2.]
cirange = [-1., 1.]

priorstr = """'H0 uniform {} {}\\nPHI0 uniform {} {}\\nPSI uniform {} {}\\nCOSIOTA uniform {} {}'
""".format(h0range[0], h0range[1], phi0range[0], phi0range[1], psirange[0], psirange[1], cirange[0], cirange[1])

# set fixed h0 value (that will be scale to SNR)
h0fixed = 1e-24

lppenexec = '/home/sismp2/lscsoft/.virtualenvs/lalapps_knope_O2/bin/lalapps_pulsar_parameter_estimation_nested'
n2pexec = '/home/sismp2/lscsoft/.virtualenvs/lalapps_knope_O2/bin/lalapps_nest2pos'

datasigma = [1e-22, 1e-22] # data standard deviation

for snr in snrs:
  # generate RA and dec values uniformly over the sky
  ras = 2.*np.pi*np.random.rand(nsigs)
  decs = -(np.pi/2.) + np.arccos(2.*np.random.rand(nsigs) - 1.)

  snrdir = os.path.join(basedir, 'snr%.1f' % snr)
  if not os.path.isdir(snrdir):
    os.makedirs(snrdir)

  for i in range(nsigs):
    cprun = ConfigParser() # set config parser to output .ini configuration file

    # set executables configuration
    cprun.add_section('executables')
    cprun.set('executables', 'lppen', lppenexec)
    cprun.set('executables', 'n2p', n2pexec)

    # set run information
    outdir = os.path.join(snrdir, '%04d' % i)
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
    cprun.set('data', 'start', json.dumps([900000000, 900000000]))
    cprun.set('data', 'step', json.dumps([60, 60]))
    cprun.set('data', 'length', json.dumps([1440, 1440]))
    
    # set nested sampling parameters
    cprun.add_section('nestedsampling')
    cprun.set('nestedsampling', 'nlive', str(nlive))
    cprun.set('nestedsampling', 'priordata', priorstr)
    
    # set injection parameters
    if snr == 0.:
      thish0 = 0.
    else:
      thish0 = h0fixed
    thisphi0 = phi0range[0] + np.diff(phi0range)[0]*np.random.rand() # phi0 value
    thispsi = psirange[0] + np.diff(psirange)[0]*np.random.rand()    # psi value
    thisci = cirange[0] + np.diff(cirange)[0]*np.random.rand()       # cos(iota) value
    cprun.add_section('injection')
    cprun.set('injection', 'snr', str(snr))
    cprun.set('injection', 'injparams', json.dumps([{'RAJ': rastr, 'DECJ': decstr, 'H0': thish0, 'PHI0': thisphi0, 'PSI': thispsi, 'COSIOTA': thisci}]))
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
