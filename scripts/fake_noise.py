#!/usr/bin/env python

"""
Script to generate some fake noise (in the style of a fine heterodyned data file),
run lalapps_pulsar_parameter_estimation and lalapps_pulsar_parameter_estimation_nested
on it and compare the results.
"""

from __future__ import print_function, division

import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
import subprocess as sp
import os
import sys
import argparse
import json
from time import time
from lalapps.pulsarpputils import *

description = "Compare lalapps_pulsar_parameter_estimation_nested with a grid-based posterior"

parser = argparse.ArgumentParser( description = description )

parser.add_argument("-r", "--run-path", dest="runpath", default='.', help="Set the run directory [default: '%(default)s']")
parser.add_argument("-o", "--out-path", dest="outpath", default='output', help="Set the output directory (within the run directory) [default: '%(default)s']")
parser.add_argument("-d", "--detector", dest="detector", default='H1', help="Set the detector [default: '%(default)s']")
parser.add_argument("-p", "--data-path", dest="datapath", default='data', help="Set the data directory (within the run directory) [default: '%(default)s']")
parser.add_argument("-n", "--nlive", dest="nlive", type=int, default=2048, help="Set the number of live points [default: %(default)d]")

opts = parser.parse_args()

# set the run directory
rundir = opts.runpath
if not os.path.isdir(rundir): # make the directory
  os.makedirs(rundir)

detector = opts.detector

# set the output directory
outdir = os.path.join(rundir, opts.outpath)
if not os.path.isdir(outdir):
  os.makedirs(outdir)

# set the output json file
jsonout = os.path.join(outdir, 'data.json')
jsondic = {}
jsondic['h0uls'] = {}
jsondic['evrats'] = {}

# fake heterodyned data directory
datadir = os.path.join(rundir, opts.datapath)
if not os.path.isdir(datadir): # make the directory
  os.makedirs(datadir)

# create random sky positions
rarad = 2.*np.pi*np.random.rand()
decrad = np.arccos(-1.+2.*np.random.rand()) - np.pi/2.
rah, ram, ras = rad_to_hms(rarad)
decd, decm, decs = rad_to_dms(decrad)

if decd > 0:
  decds = '+%02d' % decd
else:
  decds = '%02d' % decd
psrname = 'J%02d%02d%s%02d' % (int(rah), int(ram), decds, int(decm))

datafile = os.path.join(datadir, 'data'+detector)
if not os.path.isdir(datafile): # make the directory
  os.makedirs(datafile)
datafile = os.path.join(datafile, 'finehet_'+psrname+'_'+detector)

# set data standard deviation
sigma = 1.0e-22

# create the fake data
gpstimes = np.arange(900000000, 900000000+864000, 60)
dlen = len(gpstimes)
data = sigma*np.random.randn(dlen, 2)

# get an estimate of the 95% credible upper limit to be expected
ulest = 10.8*np.sqrt(sigma**2/dlen)

# append times and data together
tad = np.vstack((gpstimes, data.T)).T

# output fake data
np.savetxt(datafile, tad, fmt='%.6f %.7e %.7e', delimiter='\t')

# create a fake pulsar parameter file
parfile = os.path.join(datadir, '%s.par' % psrname)
pardat = """PSRJ {}
RAJ {}
DECJ {}
F0 123.4567890
PEPOCH 56789.0
EPHEM DE405
"""
fp = open(parfile, 'w')
fp.write(pardat.format(psrname, coord_to_string(rah, ram, ras), coord_to_string(decd, decm, decs)))
fp.close()

# create a prior file (PHI0 in here is rotational phase, whereas for the older code it is GW phase for trixial emission l=m=2)
priorfile = os.path.join(datadir, '%s.prior' % psrname)
priordat = """H0 uniform 0 {}
PHI0 uniform 0 {}
PSI uniform {} {}
COSIOTA uniform -1 1
"""

fp = open(priorfile, 'w')
fp.write(priordat.format(ulest*10., np.pi, -np.pi/4., np.pi/4.))
fp.close()

### RUN lalapps_pulsar_parameter_estimation_nested
# set the executables (this assumes that you are using virtual environments with virtualenvwrapper.sh and
# have a WORKON_HOME environment variable set, but you can change the path as required)
try:
  virenv = 'local' # name of your virtual environment
  execpath = os.path.join(os.environ['WORKON_HOME'], virenv)
  execpath = os.path.join(execpath, 'bin')
  ppenexec = os.path.join(execpath, 'lalapps_pulsar_parameter_estimation_nested')
  n2pexec = os.path.join(execpath, 'lalapps_nest2pos') # script to convert nested samples to posterior samples
  ppeexec = os.path.join(execpath, 'lalapps_pulsar_parameter_estimation')
except: # assume execs are in path
  ppenexec = 'lalapps_pulsar_parameter_estimation_nested'
  n2pexec = 'lalapps_nest2pos'
  ppeexec = 'lalapps_pulsar_parameter_estimation'

nlive = str(opts.nlive)  # number of live points
codecall = ' '.join([ppenexec, '--detectors', detector,
                    '--par-file', parfile, '--prior-file', priorfile,
                    '--input-files', datafile, '--outhdf', os.path.join(outdir, 'fake_nest.hdf'),
                    '--gzip', '--Nlive', nlive, '--Nmcmcinitial', '0', '--oldChunks'])#, '--ensembleWalk 1 --uniformprop 0'])

print(codecall)
t0 = time()
p = sp.Popen(codecall, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
out, err = p.communicate()
t1 = time()
if p.returncode != 0:
  print("Error... call returned with error:", file=sys.stderr)
  print("\tstdout: %s" % out, file=sys.stderr)
  print("\tstderr: %s" % err, file=sys.stderr)
  sys.exit(1)

print("lalapps_pulsar_parameter_estimation_nested took %f s" % (t1-t0))

# run lalapps_nest2pos to convert nested samples to posterior samples
codecall = ' '.join([n2pexec, '-p', os.path.join(outdir, 'fake_post.hdf'), os.path.join(outdir, 'fake_nest.hdf')])

print(codecall)
p = sp.Popen(codecall, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
out, err = p.communicate()

post, evsig, evnoise = pulsar_nest_to_posterior(os.path.join(outdir, 'fake_post.hdf'))
h0ul = upper_limit_greedy(post['h0'].samples, upperlimit=0.95)

jsondic['h0uls']['nested'] = h0ul
jsondic['evrats']['nested'] = evsig - evnoise

print("Nested sampling 95%% credible upper limit = %.3e, evidence ratio = %.4e" % (h0ul, evsig-evnoise))

# RUN lalapps_pulsar_parameter_estimation in grid mode
# delete any previously created evidence file as things get appended to it
evfile = os.path.join(outdir, 'evidence_%s' % psrname)
if os.path.isfile(evfile):
  os.remove(evfile)

h0steps = '100' # number of grid points for each parameter
psisteps = '40'
phi0steps = '40'
cosiotasteps = '40'
h0max = '%.5e' % (10.*ulest) # maximum range of h0 values
h0ulc = '95'                 # % credible h0 upper limit to output

codecall = ' '.join([ppeexec, '--detectors', detector,
                     '--pulsar', psrname, '--par-file', parfile, '--input-dir', os.path.join(rundir, datadir),
                     '--output-dir', outdir, '--psi-bins', '1000', '--time-bins', '1440',
                     '--h0steps', h0steps, '--maxh0', h0max, '--phi0steps', phi0steps,
                     '--psisteps', psisteps, '--cisteps', cosiotasteps, '--dob-ul', h0ulc])

t0 = time()
p = sp.Popen(codecall, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
out, err = p.communicate()
t1 = time()
timegrid = (t1-t0)

print("lalapps_pulsar_parameter_estimation took %f s" % (t1-t0))

# get upper limit and evidence ratio
# evidence at end of first line, UL at end of second
fp = open(evfile, 'r')
evlines = fp.readlines()
h0ulgrid = float((evlines[1].split())[-1])
evratgrid = float((evlines[0].split())[-1])

# lalapps_pulsar_parameter_estimation does not apply the h0 prior or cos(iota) prior, so adjust evidence accordingly
# and also account for lalapps_pulsar_parameter_estimation using a 2pi phi0 range rather than pi
evratgrid = evratgrid - np.log(10.*ulest) - np.log(2.) + np.log(np.pi)

print("Grid-based 95%% credible upper limit = %.3e, evidence ratio = %.4e" % (h0ulgrid, evratgrid))

jsondic['h0uls']['grid'] = h0ulgrid
jsondic['evrats']['grid'] = evratgrid

### RUN python-ised grid-based code
#h0steps = 80 # number of grid points for each parameter
#psisteps = 25
#phi0steps = 25
#cosiotasteps = 25
#h0max = ulest*10.

#datacomp = {detector: data[:,0] + 1j*data[:,1]}
#tsdic = {detector: gpstimes}
#ra = rarad
#dec = decrad

#paramranges = {'h0': (0., h0max, h0steps),
#               'phi0': (0., np.pi, phi0steps),
#               'cosiota': (-1., 1., cosiotasteps),
#               'psi': (0., np.pi/2., psisteps)}

#t0 = time()
#L, h0pdf, phi0pdf, psipdf, cosiotapdf, lingrid, evrat = pulsar_posterior_grid(detector, tsdic, datacomp, ra, dec, paramranges=paramranges)
#t1 = time()

#print("Python grid-mode took %f s" % (t1-t0))
#print("Evidence ratio = %.12e" % evrat)

#ct = cumtrapz(h0pdf, lingrid['h0'])/np.trapz(h0pdf, lingrid['h0'])
#ctu, ui = np.unique(ct, return_index=True)
#intf = interp1d(ctu, lingrid['h0'][ui], kind='linear')

#jsondic['h0uls']['grid'] = float(intf(0.95))
#jsondic['evrats']['grid'] = evrat

fpjson = open(jsonout, 'w')
json.dump(jsondic, fpjson, indent=2)
fpjson.close()

# clean up nested samples posteriors and data files
#os.remove(os.path.join(outdir, 'fake_nest.hdf'))
#os.remove(os.path.join(outdir, 'fake_post.hdf'))
#os.remove(datafile)
