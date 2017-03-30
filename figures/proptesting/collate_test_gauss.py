#!/usr/bin/env python

"""
Script to collate the results from running the analysis on the test Gaussian runs
"""

from __future__ import print_function, division

import sys
import os
import numpy as np
import argparse
import subprocess as sp
from scipy.special import erf
import tarfile

import matplotlib as mpl
from matplotlib import pyplot as pl
import matplotlib.patches as mpatches

mplparams = { \
      'backend': 'Agg',
      'text.usetex': True, # use LaTeX for all text
      'axes.linewidth': 0.5, # set axes linewidths to 0.5
      'axes.grid': True, # add a grid
      'grid.linewidth': 0.5,
      'font.family': 'sans-serif',
      'font.sans-serif': 'Avant Garde, Helvetica, Computer Modern Sans serif',
      'font.size': 15 }

mpl.rcParams.update(mplparams)

parser = argparse.ArgumentParser( )
parser.add_argument("-i", "--input-file", dest="infile", required=True, help="The input tarball file")
#parser.add_argument("-i", "--ndirs", dest="ndirs", required=True, type=int, help="Set the number directories used in the runs.")
#parser.add_argument("-b", "--base-dir", dest="basedir", help="Set the base directory for the runs")
parser.add_argument("-o", "--output-file", dest="outfile", help="Set the output figure file without extension")
#parser.add_argument("-N", "--Nlive", dest="nlives", action='append', default=None, help="Set number of nested samples used (each value passed will be assumed as the base directory name for that data")

# parse input options
opts = parser.parse_args()

# the base directory
#basedir = opts.basedir
#if not os.path.isdir(basedir):
#  print("Error... base directory '%s' does not exist." % basedir, file=sys.stderr)
#  sys.exit(1)

#if opts.ndirs < 1:
#  print("Error... there must be a positive number of directories.", file=sys.stderr)
#  sys.exit(1)

tar = tarfile.open(opts.infile, "r:gz")

maxranges = {}
trueuls = {}
trueevs = {}
truekls = {}
evidences = {}
upperlimits = {}
evidenceerrs = {}
kspvalues = {}
timings = {}

oformat = '*_stats.txt'
pfile = 'prior.txt'
truefile = 'test_gauss.txt'

# calculate the true KL-divergence (required as the calculation used for the value in test_gauss.txt was wrong in the code version that was run for these tests)
def kltrue(maxv, sigmah, lnZ):
  prior = 1./maxv
  p_Z = np.exp(np.log(prior)- lnZ)

  L = lnZ + 0.5*np.log(2.*np.pi*sigmah**2)

  D = -(1.+2.*L)*erf(-maxv/(np.sqrt(2.)*sigmah))
  G = -(1./(np.sqrt(2.*np.pi)*sigmah))*(maxv*np.exp(-0.5*maxv**2/sigmah**2))

  return -0.25*p_Z*(D + 2.*G)

#if opts.nlives is None:
#  nlives = ['']
#else:
#  nlives = []
#  for nlive in opts.nlives:
#    try:
#      nlives.append(str(int(nlive)))
#    except:
#      print("Error... could not convert '%s' number of live points to integer" % nlive, file=sys.stderr)
#      sys.exit(1)

# hard code in run parameters
nlives = ['512', '1024', '2048', '4096', '8192']
ndirs = ['%03d' % i for i in range(11)]

for j in range(len(nlives)):
  #livedir = os.path.join(basedir, nlives[j])
  #if len(nlives[j]) == 0:
  #  lname = 'Unknown'
  #else:
  #  lname = nlives[j]
  lname = nlives[j]

  # initialise dictionaries for different numbers of live points
  maxranges[lname] = []
  trueuls[lname] = []
  trueevs[lname] = []
  truekls[lname] = []
  evidences[lname] = []
  upperlimits[lname] = []
  evidenceerrs[lname] = []
  kspvalues[lname] = []
  timings[lname] = []

  #if not os.path.isdir(livedir):
  #  print("Error... '%s' directory does not exist." % livedir, file=sys.stderr)
  #  sys.exit(1)
  
  #for i in range(opts.ndirs):
  for fdir in ndirs:
    #fdir = os.path.join(livedir, "%03d" % i)
    #if not os.path.isdir(fdir):
    #  print("Error... '%s' directory does not exist." % fdir, file=sys.stderr)
    #  continue
    
    # get directories for the given number of live points and fdir value
    lpdirs = [(tar.getmember(name), name) for name in tar.getnames() if nlives[j]+'/'+fdir  in name]
  
    #l = os.listdir(fdir)
    #nf = 0
    #for f in l: # count number of files
    #  if '_stats.txt' in f:
    #    nf += 1

    a = []
    for lpv in lpdirs:
      if lpv[0].isfile():
        if '_stats.txt' in lpv[1]:
          fp = tar.extractfile(lpv[0])
          a.append([float(v.strip()) for v in fp.readline().split()])
        elif pfile in lpv[1]:
          # get limits from prior file
          fp = tar.extractfile(lpv[0])
          l = fp.readline().split()
          maxv = float(l[-1].strip())
          maxranges[lname].append(maxv)
    
    for lpv in lpdirs: # run again so that maxv is alway defined
      if lpv[0].isfile():
        if truefile in lpv[1]:
          # get true values of evidence and upper limits
          fp = tar.extractfile(lpv[0])
          l = fp.readline().split()
          trueevs[lname].append(float(l[0].strip()))
          trueuls[lname].append(float(l[1].strip())) 
          #truekls[lname].append(float(l[2].strip())) # this value was wrongly calculated
          truekls[lname].append(kltrue(maxv, 1e-24, trueevs[lname][-1]))

    a = np.array(a)

    # concatenate output of 'stats' files and parse it
    #p = (sp.check_output("cat "+os.path.join(fdir, oformat), shell=True)).split('\n')
    #a = np.array([[float(v.strip()) for v in l.split()] for l in p if len(l.split()) == 4])

    #if a.shape[0] != nf:
    #  print("Warning... number of files ('%d') and number of values read in ('%d') is not consistent." % (nf, a.shape[0]), file=sys.stderr)

    evidences[lname].append(a[:,0])
    upperlimits[lname].append(a[:,1])
    evidenceerrs[lname].append(a[:,2])
    kspvalues[lname].append(a[:,3])

    if a.shape[1] == 5:
      timings[lname].append(a[:,4])
    else:
      timings = None

    # get limits
    #if not os.path.isfile(os.path.join(fdir, pfile)):
    #  print("Error... no prior file given in '%s'" % fdir, file=sys.stderr)
    #  sys.exit(1)

    #fp = open(os.path.join(fdir, pfile), 'r')
    #l = fp.readline().split()
    #fp.close()
    #maxv = float(l[-1].strip())
    #maxranges[lname].append(maxv)

    # get true values of evidence and upper limits
    #if not os.path.isfile(os.path.join(fdir, truefile)):
    #  print("Error... no 'true value' file given in '%s'" % fdir, file=sys.stderr)
    #  sys.exit(1)

    #fp = open(os.path.join(fdir, truefile), 'r')
    #l = fp.readline().split()
    #fp.close()
    #trueevs[lname].append(float(l[0].strip()))
    #trueuls[lname].append(float(l[1].strip()))
    #truekls[lname].append(float(l[2].strip())) # this value was wrongly calculated
    #truekls[lname].append(kltrue(maxv, 1e-24, trueevs[lname][-1]))

    #print("Mean Z = %.7e +/- %.7e, true Z = %.7e" % (np.mean(evidences[lname][-1]), np.mean(evidenceerrs[lname][-1]), trueevs[lname][-1]), file=sys.stdout)

tar.close()

# create figure for evidences
fige = pl.figure(figsize=(12,7))
axe = fige.add_subplot(111)

# create figure for upper limits
figul = pl.figure(figsize=(12,7))
axul = figul.add_subplot(111)

# create figure for the K-S test
figks = pl.figure(figsize=(12,7))
axks = figks.add_subplot(111)

# create figure for timings
if timings is not None:
  figtim = pl.figure(figsize=(12,7))
  axtim = figtim.add_subplot(111)

colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow'] # some colours

# proxy patches for legend
handles = []

# linear fits to KL divergence vs mean evidence ratio
lfits = []

# loop over live points in order
for j, nlive in enumerate(nlives):
  nlivei = int(nlive)

  handles.append(mpatches.Patch(color=colors[j], alpha=0.2, lw=2, label=nlive))

  # got linear fit to mean evidence offset vs information
  p = np.polyfit(truekls[nlive], [np.mean(ev-trueevs[nlive][i]) for i, ev in enumerate(evidences[nlive])], deg=1)
  lfits.append(p)
  print("N live: %s, linear fit ln(Z/Z_true) = %.2f + %.2f(KL-div)" % (nlive, p[1], p[0]), file=sys.stdout)

  # violin plot for evidence ratios
  logpos = np.log10(maxranges[nlive])

  # offset the error bars so they don't overlap
  if len(nlives) > 1:
    dloffset = 0.2*(logpos[1]-logpos[0])
    dlp = 0.4*(logpos[1]-logpos[0])/(len(nlives)-1.)
  else:
    dloffset = 0.
    dlp = 0.

  logposoff = logpos-dloffset+j*dlp

  vd = axe.violinplot([ev-trueevs[nlive][i] for i, ev in enumerate(evidences[nlive])], logposoff, showextrema=False, showmedians=True, widths=0.09)
  # set colors
  for ps in vd['bodies']:
    ps.set_facecolor(colors[j])
    ps.set_edgecolor(colors[j])
    ps.set_alpha(0.2)
    ps.set_lw(1)
    ps = vd['cmedians']
    ps.set_color(colors[j])
  yerr = [np.mean(ee) for ee in evidenceerrs[nlive]]

  axe.errorbar(logposoff, np.zeros(len(logpos)), yerr=[yerr, yerr], fmt='o', capsize=3, capthick=1, color=colors[j], markersize=3)

  # add new log-style axis labels
  if j == 0:
    axe.set_xlim([logpos[0]-1, logpos[-1]+1])

    axenew = axe.twiny()

    # turn-off current (non-log) axis lables
    axe.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    axe.xaxis.grid(False)

    xlim = axe.get_xlim()
    axenew.set_xlim(10**np.array(xlim))
    axenew.set_xscale('log', nonposx='clip')
    axenew.set_xticks(np.logspace(logpos[0], logpos[-1], len(logpos)))
    axenew.xaxis.set_ticks_position('bottom')
    axenew.xaxis.grid(b=True, which='minor', color='k', linestyle='--', linewidth=0.5, alpha=0.3)

    axenew.set_xlabel('Prior range')
    axenew.xaxis.set_label_position('bottom')
    axe.set_ylabel(r'$\ln{(Z/Z_{\rm true})}$')

    # add the KL divergence (information gain)
    axenew2 = axenew.twiny()
    axenew2.grid(b=False)
    axenew2.set_xticks(logpos)
    axenew2.set_xbound(axe.get_xbound())
    axenew2.set_xticklabels(['$%.1f$' % yv for yv in truekls[nlive]])
    axenew2.set_xlabel('Information Gain (nats)')

  # show equivalent percentage evidence (not log evidence) offsets
  if j == len(nlives)-1: # only do for last set
    ayenew = axe.twinx()
    DZmin, DZmax = axe.get_ybound()
    pDZvals = np.arange(np.around(100.*(np.exp(DZmin)-1.), decimals=-1), np.around(100.*(np.exp(DZmax)-1.), decimals=-1), 50)
    equivlnZ = np.log((pDZvals/100.)+1.)

    ayenew.grid(b=False)
    ayenew.set_yticks(equivlnZ)
    ayenew.set_ybound(axe.get_ybound())

    ayenew.get_yaxis().set_tick_params(which='both', direction='out')
    ayenew.set_yticklabels(['$%d$' % yv for yv in pDZvals])
    ayenew.set_ylabel(r'$(Z-Z_{\rm true})/Z_{\rm true} \%$')

  # produce violin plot of h0 upper limits
  vd = axul.violinplot(upperlimits[nlive], logposoff, showextrema=False, showmedians=True, widths=0.09)
  for ps in vd['bodies']:
    ps.set_facecolor(colors[j])
    ps.set_alpha(0.2)
    ps.set_edgecolor(colors[j])
    ps.set_lw(1)
    ps = vd['cmedians']
    ps.set_color(colors[j])

  if j == 0:
    axul.set_xlim([logpos[0]-1, logpos[-1]+1])

    axulnew = axul.twiny()

    # turn-off current (non-log) axis lables
    axul.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    axul.xaxis.grid(False)

    # plot true analytic upper limit
    xlim = axul.get_xlim()
    axul.plot(xlim, [trueuls[nlive][0], trueuls[nlive][0]], 'r--')

    axulnew.set_xlim(10**np.array(xlim))
    axulnew.set_xscale('log', nonposx='clip')
    axulnew.set_xticks(np.logspace(logpos[0], logpos[-1], len(logpos)))
    axulnew.xaxis.set_ticks_position('bottom')
    axulnew.xaxis.grid(b=True, which='minor', color='k', linestyle='--', linewidth=0.5, alpha=0.3)

    axulnew.set_xlabel('Prior range')
    axulnew.xaxis.set_label_position('bottom')
    axul.set_ylabel(r'95\% credible upper limit')

  # convert into a fractional upper limit difference
  if j == len(nlives)-1:
    ayulnew = axul.twinx()
    ulmin, ulmax = axul.get_ybound()
    pdulvals = np.flipud(-np.arange(0., -np.floor(100.*(ulmin-trueuls[nlive][0])/trueuls[nlive][0]), 2))
    pdulvals = np.concatenate((pdulvals, np.arange(2., np.floor(100.*(ulmax-trueuls[nlive][0])/trueuls[nlive][0]), 2)))
    equivuls = ((pdulvals/100.)*trueuls[nlive][0]) + trueuls[nlive][0]

    ayulnew.grid(b=False)
    ayulnew.set_yticks(equivuls)
    ayulnew.set_ybound(axul.get_ybound())

    ayulnew.get_yaxis().set_tick_params(which='both', direction='out')
    ayulnew.set_yticklabels(['$%d$' % yv for yv in pdulvals])
    ayulnew.set_ylabel(r'$({\rm UL}-{\rm UL}_{\rm true})/{\rm UL}_{\rm true} \%$')

  # produce violin plots of K-S test p-values
  vd = axks.violinplot([np.log10(v[v > 0.0]) for v in kspvalues[nlive]], logposoff, showextrema=False, showmedians=True, widths=0.09)
  for ps in vd['bodies']:
    ps.set_facecolor(colors[j])
    ps.set_alpha(0.2)
    ps.set_edgecolor(colors[j])
    ps.set_lw(1)
    ps = vd['cmedians']
    ps.set_color(colors[j])

  if j == 0:
    axks.set_xlim([logpos[0]-1, logpos[-1]+1])

    axksnew = axks.twiny()

    # turn-off current (non-log) axis lables
    axks.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    axks.xaxis.grid(False)

    xlim = axks.get_xlim()
    axksnew.set_xlim(10**np.array(xlim))
    axksnew.set_xscale('log', nonposx='clip')
    axksnew.set_xticks(np.logspace(logpos[0], logpos[-1], len(logpos)))
    axksnew.xaxis.set_ticks_position('bottom')
    axksnew.xaxis.grid(b=True, which='minor', color='k', linestyle='--', linewidth=0.5, alpha=0.3)

    axksnew.set_xlabel('Prior range')
    axksnew.xaxis.set_label_position('bottom')

    ayksnew = axks.twinx()
    axks.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ylim = axks.get_ylim()
    ayksnew.set_ylim(10**np.array(ylim))
    ayksnew.set_yscale('log', nonposx='clip')
    ayksnew.yaxis.set_ticks_position('left')
    ayksnew.set_ylabel(r'KS test $p$-values')
    ayksnew.yaxis.set_label_position('left')

  logposoff = logpos-dloffset+j*dlp

  if timings is not None:
    # offset the error bars so they don't overlap
    #if j == 0:
    #  if len(nlives) > 1:
    #    dkloffset = 0.2*(truekls[nlive][1]-truekls[nlive][0])
    #    dklp = 0.4*(truekls[nlive][1]-truekls[nlive][0])/(len(nlives)-1.)
    #  else:
    #    dkloffset = 0.
    #    dklp = 0.
    #kloff = truekls[nlive]-dkloffset+j*dklp
    
    timescaling = 5.9e-6
    #print(kloff)
    #print(np.array(timings[nlive]).shape)
    vd = axtim.violinplot(np.array(timings[nlive]).T/timescaling, truekls[nlive], showextrema=False, showmedians=True)
    for ps in vd['bodies']:
      ps.set_facecolor(colors[j])
      ps.set_alpha(0.2)
      ps.set_edgecolor(colors[j])
      ps.set_lw(1)
      ps = vd['cmedians']
      ps.set_color(colors[j])
      
    # get linear fit to mean evidence offset vs information
    p = np.polyfit(truekls[nlive], [np.median(tims)/timescaling for tims in timings[nlive]], deg=1)
    print("N live: %s, linear fit T = %.2f + %.2f(KL-div)" % (nlive, p[1], p[0]), file=sys.stdout)

# add legend
axe.legend(handles=handles, loc='upper left')
axul.legend(handles=handles, loc='best')
axks.legend(handles=handles, loc='lower left')

fige.tight_layout()
figul.tight_layout()
figks.tight_layout()

fige.savefig(opts.outfile+'_evidences.png', dpi=300)
fige.savefig(opts.outfile+'_evidences.pdf')
p = sp.Popen('pdftops -eps %s' % opts.outfile+'_evidences.pdf', shell=True)
p.communicate()

figul.savefig(opts.outfile+'_uls.png', dpi=300)
figul.savefig(opts.outfile+'_uls.pdf')
p = sp.Popen('pdftops -eps %s' % opts.outfile+'_uls.pdf', shell=True)
p.communicate()

figks.savefig(opts.outfile+'_ks.png', dpi=300)
figks.savefig(opts.outfile+'_ks.pdf')
p = sp.Popen('pdftops -eps %s' % opts.outfile+'_ks.pdf', shell=True)
p.communicate()

if timings is not None:
  axtim.set_yscale('log')
  axtim.legend(handles=handles, loc='lower right')
  axtim.set_xlabel('Information Gain (nats)')
  axtim.set_ylabel(r'median run time ($\mathcal{T}_{L}^{-1}$)')
  figtim.tight_layout()
  figtim.savefig(opts.outfile+'_timings.png', dpi=300)
  figtim.savefig(opts.outfile+'_timings.pdf')
  p = sp.Popen('pdftops -eps %s' % opts.outfile+'_timings.pdf', shell=True)
  p.communicate()