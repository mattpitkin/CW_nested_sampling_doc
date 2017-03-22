#!/usr/bin/env python

"""
Short script to run lalapps_pulsar_parameter_estimation_nested based on inputs given in a configuration file.
The outputs will be par
"""

from __future__ import print_function, division

import os
import sys
import h5py
import numpy as np
import argparse
import json
import uuid
import ast
import glob
from scipy import stats
import subprocess as sp
from ConfigParser import ConfigParser
from copy import copy
from multiprocessing.dummy import Pool

import lalapps.pulsarpputils as pppu
from lalapps.pulsarpputils import pulsar_nest_to_posterior as pn2p
from lalapps.pulsarpputils import upper_limit_greedy as ulg


def run_process(commands):
  p = sp.Popen(commands, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
  out, err = p.communicate()
  return 0


def get_with_default(config, section, name, default):
  if config.has_option(section, name):
    return config.get(section, name)
  else:
    return default


def find_exec_file(filename):
  """
  Search through the PATH environment for the first instance of the executable file "filename"
  """

  if os.path.isfile(filename) and os.access(filename, os.X_OK):
    return filename # already have the executable file

  # search through PATH environment variable
  for d in os.environ['PATH'].split(os.pathsep):
    filecheck = os.path.join(d, filename)
    if os.path.isfile(filecheck) and os.access(filecheck, os.X_OK):
      return filecheck

  return None


def credible_interval(samples, ci=95, paramval=None):
  # for injections, where we have parameter values (paramval) get the corresponding smallest credible
  # region that contains the parameter value)
  samples.sort()
  lowbound, highbound, _ = ci_loop(samples, ci)

  cival = None
  if paramval != None:
    cifound = False
    # loop over different credible intervals until finding the one that contains the injection
    for cit in range(1, 101):
      l, h, cr = ci_loop(samples, cit)
      if paramval >= l and paramval <= h:
        cifound = True
        cival = cit
        break
    if not cifound:
      cival = 100

  return cival, [lowbound, highbound]


def ci_loop(sortedsamples, ci):
  lowbound = sortedsamples[0]
  highbound = sortedsamples[-1]
  cregion = highbound - lowbound
  lsam = len(sortedsamples)
  cllen = int(lsam*float(ci)/100.)
  for j in range(lsam-cllen):
    if sortedsamples[j+cllen] - sortedsamples[j] < cregion:
      lowbound = sortedsamples[j]
      highbound = sortedsamples[j+cllen]
      cregion = highbound - lowbound

  return lowbound, highbound, cregion


def get_snr(pfile):
  # get SNR from file in pdir
  snr = 0.

  fp = open(pfile, 'r')
  lines = [line.strip() for line in fp.readlines()]

  if '# Recovered SNR' not in lines:
    print("Error... no recovered SNRs are given in the SNR file '%s'." % snrfile, file=sys.stderr)
    sys.exit(1)

  # just return the final (either coherent or single detector-single frequency value of the SNR)
  linevals = lines[-1].split()
  snr += float(linevals[-1])

  return snr


if __name__=='__main__':
  parser = argparse.ArgumentParser( )
  parser.add_argument("inifile", help="The configuration (.ini) file")
  parser.add_argument("-r", "--dont-remove", dest="drm", action='store_true', default=False, help="Set to NOT remove files (will be removed by default)")

  # parse input options
  opts = parser.parse_args()

  inifile = opts.inifile
  if not os.path.isfile(inifile):
    print("Error... configuration .ini file does not exist", file=sys.stderr)
    sys.exit(1)

  # open and parse config file
  cp = ConfigParser()
  try:
    cp.read(inifile)
  except:
    print("Error... cannot parse configuration file '%s'" % inifile, file=sys.stderr)
    sys.exit(1)

  # get information from the configuration file

  # get executables
  ppen = get_with_default(cp, 'executables', 'lppen', 'lalapps_pulsar_parameter_estimation_nested')
  ppen = find_exec_file(ppen)
  if ppen is None:
    print("Error... cannot find 'lalapps_pulsar_parameter_estimation_nested' executable", file=sys.stderr)
    sys.exit(1)

  n2p = get_with_default(cp, 'executables', 'n2p', 'lalapps_nest2pos')
  n2p = find_exec_file(n2p)
  if n2p is None:
    print("Error... cannot find 'lalapps_nest2pos' executable", file=sys.stderr)
    sys.exit(1)

  # get run information
  outdir = get_with_default(cp, 'run', 'outdir', None)
  if outdir is None:
    print("Error... an output directory is required", file=sys.stderr)
    sys.exit(1)

  # check if the directory exists
  if not os.path.isdir(outdir):
    try:
      os.makedirs(outdir)
    except:
      print("Error... could not make output directory '%s'" % outdir, file=sys.stderr)
      sys.exit(1)

  outname = get_with_default(cp, 'run', 'outname', 'nest'+uuid.uuid4().hex) # default to random string

  # get number of nested sample runs
  nruns = ast.literal_eval(get_with_default(cp, 'run', 'nruns', '1'))
  if nruns < 1:
    print("Warning... number of nested sample runs is less than 1, default to 1 run", file=sys.stderr)
    nruns = 1

  # get number of CPUs to run on
  ncpus = ast.literal_eval(get_with_default(cp, 'run', 'ncpus', '1'))

  # get the detectors to use
  dets = ast.literal_eval(get_with_default(cp, 'run', 'detectors', "['H1']"))
  idets = copy(dets)
  ndets = len(dets)
  if len(dets) > 1:
    dets.append(','.join(dets))

  # output files
  outnest = []
  outpost = []
  for j in range(len(dets)):
    outnestdet = []
    for i in range(nruns):
      outnestdet.append(os.path.join(outdir, outname+'_%s_%04d.hdf' % (dets[j].replace(',', ''), i))) # nested sample files
    outnest.append(outnestdet)
    outpost.append(os.path.join(outdir, outname+'_%s_post.hdf' % dets[j].replace(',', ''))) # posterior sample file

  # get the par file if exists, otherwise create it
  parfe = False
  parfile = get_with_default(cp, 'run', 'parfile', os.path.join(outdir, 'pulsar.par'))
  defaultpar = """'PSRJ J0000+0000\\n\
RAJ {RAJ}\\n\
DECJ {DECJ}\\n\
F0 {F0}\\n\
F1 {F1}\\n\
PEPOCH {PEPOCH}'"""

  if not os.path.isfile(parfile):
    parparams = "{'RAJ': '00:00:00.0', 'DECJ': '00:00:00.0', 'F0': 100, 'F1': -1e-9, 'PEPOCH': 54660}"
    hetpars = ast.literal_eval(get_with_default(cp, 'run', 'hetparams', parparams)) 
  
    if not isinstance(hetpars, dict):
      print("Error... heterodyne parameters must be a dictionary" % outdir, file=sys.stderr)
      sys.exit(1)
    else:
      pp = ast.literal_eval(parparams)
      for key in pp.keys():
        if key not in hetpars:
          hetpars[key] = pp[key]  

    pardata = ast.literal_eval(get_with_default(cp, 'run', 'pardata', defaultpar.format(**hetpars)))
    fp = open(parfile, 'w')
    fp.write(pardata)
    fp.close()
  else:
    parfe = True

  # read in parameters from par file for later
  hetpars = pppu.psr_par(parfile).__dict__

  # get simulated data parameters (for each detector specified), or pre-made data files
  datasigma = ast.literal_eval(get_with_default(cp, 'data', 'sigma', '[1e-22]'))     # noise standard deviation
  datastart = ast.literal_eval(get_with_default(cp, 'data', 'start', '[900000000]')) # start GPS time
  datatstep = ast.literal_eval(get_with_default(cp, 'data', 'step', '[60]'))         # data time step
  datalength = ast.literal_eval(get_with_default(cp, 'data', 'length', '[1440]'))      # number of data points
  datafiles = ast.literal_eval(get_with_default(cp, 'data', 'files', "['%s']" % os.path.join(outdir, 'data_{}.txt.gz'))) # pre-made data files, or names of files to be used

  # check if data files exist
  dfe = False
  if len(datafiles) == ndets:
    if np.array([os.path.isfile(datafiles[i]) for i in range(ndets)]).all():
      dfe = True # files exist

  if not dfe:
    if len(datafiles) == 1 and ndets > 1:
      df = datafiles[0]
      datafiles = [df for i in range(ndets)]

    if len(datasigma) == 1 and ndets > 1:
      ds = datasigma[0]
      datasigma = [ds for i in range(ndets)]

    if len(datastart) == 1 and ndets > 1:
      ds = datastart[0]
      datastart = [ds for i in range(ndets)]

    if len(datatstep) == 1 and ndets > 1:
      ds = datatstep[0]
      datatstep = [ds for i in range(ndets)]

    if len(datalength) == 1 and ndets > 1:
      ds = datalength[0]
      datalength = [ds for i in range(ndets)]

    if not np.all(np.array([len(datasigma), len(datasigma), len(datatstep), len(datalength)]) == ndets):
      print("Error... simulated data parameter numbers must match the number of detectors.", file=sys.stderr)
      sys.exit(1)

  # create the data (if files don't already exist)
  datafilesdict = {}
  injectoutputs = {}    # files output from each individual detector
  injjoint = []
  for j, d in enumerate(idets):
    datafilesdict[d] = datafiles[j].format(d)

    if not dfe:
      gpstimes = np.arange(datastart[j], datastart[j]+(datalength[j]-1.)*datatstep[j], datatstep[j])
      dlen = len(gpstimes)
      data = datasigma[j]*np.random.randn(dlen, 2)

      # append times and data together
      tad = np.vstack((gpstimes, data.T)).T
      np.savetxt(datafilesdict[d], tad, fmt='%.6f %.7e %.7e', delimiter='\t')
      injectoutputs[d] = os.path.join(outdir, 'injection_data.txt')
      injjoint.append(injectoutputs[d]+'_%s_2.0' % d) # append output format of actual files
    else:
      datafilesdict[d] = datafiles[j]

  if dfe:
    indfile = ','.join(datafilesdict.values())
  else:
    indfile = ','.join(injjoint)

  if len(dets) > 1:
    datafilesdict[dets[-1]] = indfile # add all data sets for joint analysis

  # get nested sampling parameters
  nlive = get_with_default(cp, 'nestedsampling', 'nlive', 1024)

  # get proposal ratios
  uniformprop = get_with_default(cp, 'nestedsampling', 'uniformprop', 1)
  walkprop = get_with_default(cp, 'nestedsampling', 'walkprop', 3)
  stretchprop = get_with_default(cp, 'nestedsampling', 'stretchprop', 0)

  # get prior file if exists, otherwise create prior
  priorfile = get_with_default(cp, 'nestedsampling', 'priorfile', os.path.join(outdir, 'prior.txt'))
  defaultprior = """'H0 uniform 0.0 1e-20\\n\
PHI0 uniform 0.0 {}\\n\
PSI uniform 0.0 {}\\n\
COSIOTA uniform -1.0 1.0\\n'"""
  priorfe = False
  if not os.path.isfile(priorfile):
    # check if specified prior file exists, and if not create it
    priordata = ast.literal_eval(get_with_default(cp, 'nestedsampling', 'priordata', defaultprior.format(np.pi, np.pi/2.)))
    fp = open(priorfile, 'w')
    fp.write(priordata)
    fp.close()
  else:
    priorfe = True

  # get the variable parameters from the prior file
  varpars = []
  fp = open(priorfile, 'r')
  for line in fp.readlines():
    if line[0] == '#' or len(line) == 0:
      continue
    else:
      varpars.append(line.split()[0].strip())

  # check if using ROQ
  roq = ast.literal_eval(get_with_default(cp, 'nestedsampling', 'roq', 'False'))
  roqtraining = ast.literal_eval(get_with_default(cp, 'nestedsampling', 'roqntraining', '2500'))
  roqtolerance = ast.literal_eval(get_with_default(cp, 'nestedsampling', 'roqtolerance', '5e-12'))

  # simulated injection information
  injsnr = ast.literal_eval(get_with_default(cp, 'injection', 'snr', 0.0)) # injection SNR (scale injection to this SNR)
  injpar = ast.literal_eval(get_with_default(cp, 'injection', 'parfile', "['%s']" % os.path.join(outdir, 'inj_{}_pulsar.par'))) # injection par file for each detector (can be different if wanting incoherent signals)
  defaultinjpar = """'PSRJ J0000+0000\\n\
RAJ {RAJ}\\n\
DECJ {DECJ}\\n\
F0 {F0}\\n\
F1 {F1}\\n\
PEPOCH {PEPOCH}\\n\
H0 {H0}\\n\
COSIOTA {COSIOTA}\\n\
PSI {PSI}\\n\
PHI0 {PHI0}\\n'"""
  defaultinjparams = "[{'RAJ': '00:00:00.0', 'DECJ': '00:00:00.0', 'F0': 100, 'F1': -1e-9, 'PEPOCH': 54660, 'H0': 1e-24, 'COSIOTA': 0.0, 'PSI': 0.5, 'PHI0': 1.4}]"
  injparfe = [False for i in range(ndets)]

  # get injection values
  if len(injpar) == 1: # use the same injection file for all detectors
    injpar = [injpar[0] for i in range(ndets)]

  ipe = False
  if len(injpar) == ndets:
    if np.array([os.path.isfile(datafiles[i]) for i in range(ndets)]).all():
      ipe = True # files exist

  if not ipe:
    dip = ast.literal_eval(defaultinjparams)[0]
    injpardet = ast.literal_eval(get_with_default(cp, 'injection', 'injparams', defaultinjparams)) 
  
    if len(injpardet) == 1:
      injpardet = [injpardet[0] for i in range(ndets)]

    for i in range(ndets):
      if not isinstance(injpardet[i], dict):
        print("Error... injection parameters must be a dictionary" % outdir, file=sys.stderr)
        sys.exit(1)
      else:
        for key in dip:
          if key not in injpardet[i]:
            injpardet[i][key] = dip[key]

    injpardata = ast.literal_eval(get_with_default(cp, 'injection', 'injdata', "["+defaultinjpar+"]"))
    if len(injpardata) == 1:
      injpardata = [injpardata[0] for i in range(ndets)]

    for i in range(ndets):
      fp = open(injpar[i].format(dets[i]), 'w')
      fp.write(injpardata[i].format(**injpardet[i]))
      fp.close()

  # calculate SNRs and scale amplitudes appropriately
  ampscale = 0.
  ampscales = []
  injsnrs = {}
  totsnr = 0.
  for i in range(ndets):
    pdict = {}
    injd = pppu.psr_par(injpar[i].format(dets[i])) # read in injection parameter file

    pdict['ra'] = injd['RAJ']
    pdict['dec'] = injd['DECJ']
    pdict['h0'] = injd['H0']
    pdict['psi'] = injd['PSI']
    pdict['cosiota'] = injd['COSIOTA']
    pdict['phi0'] = injd['PHI0']

    if not dfe:
      dtimes = np.arange(datastart[i], datastart[i]+datalength[i]*datatstep[i], datatstep[i])
      dsigma = datasigma[i]
    else: # get data times and sigma from actual data file
      try:
        data = np.loadtxt(datafilesdict[dets[i]], comments=['#', '%'])
      except:
        print("Error... could not load existing data file '%s'".format(datafilesdict[dets[i]]), file=sys.stderr)
        sys.exit(1)
      dtimes = data[:,0]
      dsigma = np.mean([np.std(data[:,1]), np.std(data[:,2])])

    ts, s = pppu.heterodyned_pulsar_signal(pdict, dets[i], datatimes=dtimes)
    thissnr = pppu.get_optimal_snr( s[0], dsigma )
    
    if injsnr != 0.:
      ampscales.append(thissnr**2)
      ampscale += thissnr**2
    else:
      ampscales.append(1.)
      ampscale = 1.

    if injsnr == 0.: # hold injected SNRs (calculated from the given amplitude)
      injsnrs[dets[i]] = thissnr
      totsnr += thissnr**2

  if injsnr == 0. and ndets > 1:
    injsnrs[dets[-1]] = np.sqrt(totsnr)

  ampscale = injsnr/np.sqrt(ampscale) # overall snr scale factor will be for the multi-detector signal
  for i in range(ndets):
    injd = pppu.psr_par(injpar[i].format(dets[i])).__dict__ # read in injection parameter file

    # scale any amplitude parameters
    if injsnr != 0.:
      for key in ['H0', 'C22', 'C21']:
        if key in injd:
          injd[key] *= ampscale

      fp = open(injpar[i].format(dets[i]), 'w')
      for key in injd:
        if '_ORIGINAL' in key:
          continue
        else:
          if key+'_ORIGINAL' in injd: # output 'ORIGINAL' (non-unit-converted) version if present
            fp.write('{} {}\n'.format(key, injd[key+'_ORIGINAL']))
          else:
            fp.write('{} {}\n'.format(key, injd[key]))
      fp.close()

      if ndets > 1:
        asn = (injsnr/ampscale)**2
        injsnrs[dets[i]] = np.sqrt((injsnr**2/asn)*(asn-sum([x for j,x in enumerate(ampscales) if j!=i])))
      else:
        injsnrs[dets[i]] = injsnr

  if injsnr != 0. and ndets > 1:
    injsnrs[dets[-1]] = injsnr

  injcoh = ast.literal_eval(get_with_default(cp, 'injection', 'coherent', 'True')) # check whether injection is coherent between detectors

  ppencodecall = "{} --detectors {} --input-files {} --prior-file {} --par-file {} --uniformprop {} --ensembleWalk {} --ensembleStretch {} --outfile {} --Nlive {} --Nmcmcinitial 0 --tolerance 0.1 {}"
  n2pcodecall = "{} -p {} {}"

  Zs = {}        # signal evidences
  Ns = {}        # noise evidences
  uls = {}       # 95% h0 upper limits
  Zerrs = {}     # estimated errors on the evidence
  recSNRs = {}   # recovered SNR
  cis = {}
  cisp = {}
  ainjpars = {}
  for j in range(len(dets)):
    ci = {}
    cip = {}

    extraargs = ""
    if j < ndets and not dfe: # doing single detector analyses
      extraargs += "--inject-file {} --inject-output {}".format(injpar[j].format(dets[j]), injectoutputs[dets[j]])

    if roq:
      extraargs += " --roq --ntraining {} --roq-tolerance {}".format(roqtraining, roqtolerance)

    snrs = []

    # run on multiple CPUs
    npools = ncpus if nruns > ncpus else nruns
    pool = Pool(npools) # ncpu concurrent commands at a time
    commands = []
    for i in range(nruns):
      commands.append(ppencodecall.format(ppen, dets[j], datafilesdict[dets[j]], priorfile, parfile, uniformprop, walkprop, stretchprop, outnest[j][i], nlive, extraargs))

    # run pool for lalapps_pulsar_parameter_estimation_nested
    pool.map(run_process, commands)
    pool.close()

    for i in range(nruns):    
      # get SNRs
      snrfile = os.path.splitext(outnest[j][i])[0]+'_SNR'
      snr = get_snr(snrfile)
      snrs.append(snr)

    recSNRs[dets[j]] = np.mean(snrs)

    # run lalapps_nest2pos
    p = sp.Popen(n2pcodecall.format(n2p, outpost[j], ' '.join(outnest[j])), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    out, err = p.communicate()

    try:
      post, Z, N = pn2p(outpost[j])
      Zs[dets[j]] = Z
      Ns[dets[j]] = N
    except:
      print("Error... could not extract samples from posterior file '%s'" % outpost[j], file=sys.stderr)
      sys.exit(1)

    # get 95% upper limit for H0 samples
    try:
      uls[dets[j]] = ulg(post['H0'].samples)
    except:
      print("Error... could not get 'H0' upper limit from posterior samples", file=sys.stderr)
      sys.exit(1)

    # get error on the Z value
    hdf = h5py.File(outpost[j], 'r')
    a = hdf['lalinference']['lalinference_nest']
    info = a.attrs['information_nats']
    nlive = a.attrs['number_live_points']
    Zerrs[dets[j]] = np.sqrt(info/nlive)

    # injection parameters
    if j < ndets:
      ainjpars[dets[j]] = pppu.psr_par(injpar[j].format(dets[j])).__dict__

    # get smallest credible interval containing the injection (for p-p plots) - only do this for coherent injections
    if injcoh or ndets == 1:
      # get heterodyne parameters from par file and injection parameters from par file
      thetpars = []
      tinjpars = []
      for vp in varpars:
        vpa = vp
        if vp == 'RAJ' or vp == 'DECJ' or vp == 'RA' or vp == 'DEC':
          vpa = vp.strip('J') + '_RAD'

        tinjpars.append(ainjpars[dets[0]][vpa]) # get from first detector as all values should be the same for coherent injections

        # credible interval
        ci[vp], cip[vp] = credible_interval(post[vp].samples[:,0], paramval=ainjpars[dets[0]][vpa])

      cis[dets[j]] = ci   # interval in which injection is found
      cisp[dets[j]] = cip # 95% credible bounds

  # output dictionary
  outdic = {}
  outdic['Detectors'] = dets
  outdic['Injected SNR'] = injsnrs
  outdic['Recovered SNR'] = recSNRs
  outdic['Signal evidence'] = Zs
  outdic['Noise evidence'] = Ns
  outdic['Injection credible intervals'] = cis
  outdic[r'95% amplitude upper limits'] = uls
  outdic[r'Parameter 95% credible intervals'] = cisp
  outdic['Heterodyne parameters'] = hetpars

  # output injection parameters
  if len(ainjpars) > 0:
    outdic['Injection parameters'] = ainjpars

  outfile = get_with_default(cp, 'output', 'outfile', os.path.join(outdir, 'stats.json'))
  fp = open(outfile, 'w')
  json.dump(outdic, fp, indent=2)
  fp.close()

  # clean up files if required
  if not opts.drm:
    try:
      # remove nested sample files
      for f in outnest:
        for ff in f:
          os.remove(ff)

      # remove posterior samples
      for f in outpost:
        os.remove(f)

      # remove SNR and Znoise files
      fs = glob.glob(os.path.join(outdir, '*_SNR'))
      for f in fs:
        os.remove(f)

      fs = glob.glob(os.path.join(outdir, '*_Znoise'))
      for f in fs:
        os.remove(f)

      # remove data files (if not previously existing files)
      if not dfe:
        for d in idets:
          os.remove(datafilesdict[d])
          fs = glob.glob(injectoutputs[d]+'*')
          for f in fs:
            os.remove(f)

      # remove par file
      if not parfe:
        os.remove(parfile)

      # remove prior file
      if not priorfe:
        os.remove(priorfile)

      # remove injection files
      if not ipe:
        for i, f in enumerate(injpar):
          os.remove(f.format(dets[i]))
    except:
      print("Warning... could not remove various files.", file=sys.stderr)

  sys.exit(0)
