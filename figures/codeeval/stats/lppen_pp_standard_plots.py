#!/usr/bin/env python

"""
Script to generate p-p plots for standard parameters
"""

import os
import sys
import json
import numpy as np
import tarfile
import subprocess as sp
import itertools

from scipy.stats import kstest

import matplotlib as mpl
from matplotlib import pyplot as pl
import matplotlib.gridspec as gridspec

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

# tarball containing (coherent) data
tarball = 'pp_standard.tar.gz'

# tarball containing (incoherent) data
itarball = 'pp_incoherent.tar.gz'

# tarball containing (coherent) higher dimensional runs
htarball = 'pp_roq.tar.gz'

tar = tarfile.open(tarball, "r:gz")
itar = tarfile.open(itarball, "r:gz")
htar = tarfile.open(htarball, "r:gz")

# number of sub-directories
ndirs = 2000

statsfile = 'stats.json'

pars = ['H0', 'PHI0', 'COSIOTA', 'PSI'] # parameters to pp-plot
labels = ['$h_0$', '$\Phi_{22}^C$', '$\cos{\iota}$', '$\psi$']
hpars = ['H0', 'PHI0', 'COSIOTA', 'PSI', 'F0', 'F1', 'RA'] # parameters to pp-plot
hlabels = ['$h_0$', '$\Phi_{22}^C$', '$\cos{\iota}$', '$\psi$', '$f$', '$\dot{f}$', r'$\alpha$']

detectors = ['H1', 'L1', 'H1,L1']
dcolors = ['red', 'green', 'black']

# snrs to use
snrlim = 10 # use signals with SNRs greater than this

injcreds = []

injsnrs = []
recsnrs = []

sigev = []
noiseev = []

# COHERENT: loop through directories and grab required data (for SNRs above limit)
for m in tar.getmembers():
  if m.isfile():
    if statsfile in m.name:
      fp = tar.extractfile(m)
      d = json.load(fp)

      # get injection credible intervals for joint analysis
      if d['Injected SNR'][detectors[-1]] > snrlim:
        injcreds.append(d['Injection credible intervals'][detectors[-1]])
    
      # injected and recovered SNRs
      injsnrs.append(d['Injected SNR'])
      recsnrs.append(d['Recovered SNR'])

      # noise evidences and signal evidences
      sigev.append(d['Signal evidence'])
      noiseev.append(d['Noise evidence'])

tar.close() # close tarball

# INCOHERENT: loop through directories and grab required data (for SNRs above limit)
isigev = []
inoiseev = []
iinjsnrs = []
irecsnrs = []

for m in itar.getmembers():
  if m.isfile():
    if statsfile in m.name:
      fp = itar.extractfile(m)
      d = json.load(fp)

      # injected and recovered SNRs
      iinjsnrs.append(d['Injected SNR'])
      irecsnrs.append(d['Recovered SNR'])

      # noise evidences and signal evidences
      isigev.append(d['Signal evidence'])
      inoiseev.append(d['Noise evidence'])

itar.close() # close tarball

# COHERENT HIGHER DIMENSIONAL RUNS
hsigev = []
hnoiseev = []
hinjsnrs = []
hrecsnrs = []
hinjcreds = []

for m in htar.getmembers():
  if m.isfile():
    if statsfile in m.name:
      fp = htar.extractfile(m)
      d = json.load(fp)

      # get injection credible intervals for joint analysis
      if d['Injected SNR'][detectors[-1]] > snrlim:
        hinjcreds.append(d['Injection credible intervals'][detectors[-1]])

      # injected and recovered SNRs
      hinjsnrs.append(d['Injected SNR'])
      hrecsnrs.append(d['Recovered SNR'])

      # noise evidences and signal evidences
      hsigev.append(d['Signal evidence'])
      hnoiseev.append(d['Noise evidence'])

htar.close() # close tarball

fig, axs = pl.subplots(2,2, figsize=(9,9))

nlightning = 100

# loop over parameters and produce p-p plots
for i, ax in enumerate(axs.flatten()):
  parcreds = np.array([inj[pars[i]] for inj in injcreds])/100.
  parcreds.sort() # sort into order

  # add "lightning"
  for j in range(nlightning):
    #x = np.random.randint(1, 101, len(parcreds))
    x = np.random.rand(len(parcreds))
    x.sort()
    nb, binedges = np.histogram(x, bins=len(x), normed=True)
    cn = np.cumsum(nb)/len(x)
    ax.step(binedges[:-1], cn, color='g', alpha=0.02)

  # perform K-S test (conparing to uniform distribution)
  D, pv = kstest(parcreds, lambda y: y)

  # create p-p plot histogram
  nb, binedges = np.histogram(parcreds, bins=len(parcreds), normed=True)
  cn = np.cumsum(nb)/len(parcreds)
  ax.step(binedges[:-1], cn, color='k', label='KS $p$-value: %.3f' % pv)
  ax.plot([0., 1.], [0., 1.], color='darkred', ls='--', lw=0.5)
  ax.set_xlim([0., 1.])
  ax.set_ylim([0., 1.])
  ax.set_xlabel('Credible interval (%s)' % labels[i])
  ax.set_ylabel('Fraction of true values within CI')

  ax.legend()
  
fig.tight_layout()

# output P-P plots
ppdir = os.path.join(os.getcwd(), 'pp_standard')
if not os.path.isdir(ppdir):
  os.makedirs(ppdir)
fig.savefig(os.path.join(ppdir, 'pp_plots_standard.png'), dpi=300)
fig.savefig(os.path.join(ppdir, 'pp_plots_standard.pdf'))

p = sp.Popen('pdftops -eps %s' % os.path.join(ppdir, 'pp_plots_standard.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:pp_standard}
The cumulative fraction of true parameter values (from simulated signals) found within
a particular credible interval.
"""

fp = open(os.path.join(ppdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()

# produce larger p-p plots
pl.clf()

fig, axs = pl.subplots(3,3, figsize=(13,13))

# loop over parameters and produce p-p plots
counter = 0
axf = axs.flatten()
axi = 1
axf[0].axis('off')
axf[2].axis('off')
for i in range(len(hpars)):
  parcreds = np.array([inj[hpars[i]] for inj in hinjcreds])/100.
  parcreds.sort() # sort into order

  # add "lightning"
  for j in range(nlightning):
    #x = np.random.randint(1, 101, len(parcreds))
    x = np.random.rand(len(parcreds))
    x.sort()
    nb, binedges = np.histogram(x, bins=len(x), normed=True)
    cn = np.cumsum(nb)/len(x)
    axf[axi].step(binedges[:-1], cn, color='g', alpha=0.02)

  # perform K-S test (conparing to uniform distribution)
  D, pv = kstest(parcreds, lambda y: y)

  # create p-p plot histogram
  nb, binedges = np.histogram(parcreds, bins=len(parcreds), normed=True)
  cn = np.cumsum(nb)/len(parcreds)
  axf[axi].step(binedges[:-1], cn, color='k', label='KS $p$-value: %.2f' % pv)
  axf[axi].plot([0., 1.], [0., 1.], color='darkred', ls='--', lw=0.5)
  axf[axi].set_xlim([0., 1.])
  axf[axi].set_ylim([0., 1.])
  axf[axi].text(0.85, 0.1, hlabels[i])

  if axi > 5:
    axf[axi].set_xlabel('Credible interval')
  else:
    axf[axi].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
  if axi == 3 or axi == 6:
    axf[axi].set_ylabel('Fraction of true values within CI')
  else:
    axf[axi].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')

  axf[axi].legend(loc='upper left')

  if i==0:
    axi += 2
  else:
    axi += 1
  
fig.tight_layout()


# output P-P plots
ppdir = os.path.join(os.getcwd(), 'pp_extra')
if not os.path.isdir(ppdir):
  os.makedirs(ppdir)
fig.savefig(os.path.join(ppdir, 'pp_plots_extra.png'), dpi=300)
fig.savefig(os.path.join(ppdir, 'pp_plots_extra.pdf'))

p = sp.Popen('pdftops -eps %s' % os.path.join(ppdir, 'pp_plots_extra.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:pp_extra}
The cumulative fraction of true parameter values (from simulated signals) found within
a particular credible interval.
"""

fp = open(os.path.join(ppdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()


#################################################
# Create a plot of injected SNR versus recovered SNR for each detector
pl.clf()
fig = pl.figure(figsize=(6,5))
ax = pl.gca()
msizes = [4, 4, 2]
labels = [detectors[0], detectors[1], 'Coherent']
for i, det in enumerate(detectors):
  ax.plot([isv[det] for isv in injsnrs], [rsv[det] for rsv in recsnrs], color=dcolors[i], ls='none', marker='o', ms=msizes[i], alpha=0.2, label=labels[i])

# make sure axes start at zero
xlims = ax.get_xlim()
ylims = ax.get_ylim()
ax.set_xlim([0., 20.])
ax.set_ylim([0., ylims[-1]])

ax.set_xlabel('Injected SNR')
ax.set_ylabel('Recovered SNR')

ax.legend(loc='upper left')

fig.tight_layout()

# output snr plots
ppdir = os.path.join(os.getcwd(), 'snrs')
if not os.path.isdir(ppdir):
  os.makedirs(ppdir)
fig.savefig(os.path.join(ppdir, 'snr_plot.png'), dpi=300)
fig.savefig(os.path.join(ppdir, 'snr_plot.pdf'))

p = sp.Popen('pdftops -eps %s' % os.path.join(ppdir, 'snr_plot.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:snrplot}
The injected signal signal-to-noise ratios plotted against the recovered signal-to-noise ratios (Equation~\ref{eq:snr}),
where the recovered values where calculated from the maximum likelihood signal template.
"""

fp = open(os.path.join(ppdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()


#################################################
# Create plot of coherent Bayes factors versus incoherent Bayes factors (with SNR as the colour scale)

# COHERENT SIGNALS: get the incoherent odds (with equal priors for all hypotheses)
nifos = len(detectors[:-1])
bcins = []
bcis = []
cohO = []  # coherent odds ratio
for se, ne in zip(sigev, noiseev):
  ifossn = []
  idd = 0.
  for j in range(nifos):
    ifossn.append({'s': se[detectors[j]], 'n': ne[detectors[j]]})
    idd += se[detectors[j]]
  cohO.append(se[detectors[-1]] - ne[detectors[-1]])

  combs = [list(i) for i in itertools.product(['s', 'n'], repeat=nifos)]
  incoherentcombs = -np.inf
  for comb in combs:
    combsum = 0.
    for i, cval in enumerate(comb):
      combsum += ifossn[i][cval]
    incoherentcombs = np.logaddexp(incoherentcombs, combsum)
  bcins.append(se[detectors[-1]] - incoherentcombs)
  bcis.append(se[detectors[-1]] - idd)

# INCOHERENT SIGNALS: get the incoherent odds (with equal priors for all hypotheses)
ibcins = []
icohO = []  # coherent odds ratio
for se, ne in zip(isigev, inoiseev):
  ifossn = []
  for j in range(nifos):
    ifossn.append({'s': se[detectors[j]], 'n': ne[detectors[j]]})
  icohO.append(se[detectors[-1]] - ne[detectors[-1]])

  combs = [list(i) for i in itertools.product(['s', 'n'], repeat=nifos)]
  incoherentcombs = -np.inf
  for comb in combs:
    combsum = 0.
    for i, cval in enumerate(comb):
      combsum += ifossn[i][cval]
    incoherentcombs = np.logaddexp(incoherentcombs, combsum)
  ibcins.append(se[detectors[-1]] - incoherentcombs)

pl.clf()

fig = pl.figure(figsize=(10,8))

gs = gridspec.GridSpec(2, 2, width_ratios=[25,1], height_ratios=[1,1], wspace=0.03, hspace=0.03)

ax1 = pl.subplot(gs[0])

#fig, axs = pl.subplots(2,1, figsize=(10,8))
#ax = pl.gca()

# on left plot just plot coherent results and incoherent ones with log(O_S/I) greater than range of plot. 
cohO = np.array(cohO)/np.log(10.)
bcins = np.array(bcins)/np.log(10.)
bcis = np.array(bcis)/np.log(10.)

icohO = np.array(icohO)/np.log(10.)
ibcins = np.array(ibcins)/np.log(10.)

cohsnrs = np.array([isv[detectors[-1]] for isv in injsnrs])
#hb = ax.scatter(cohO, bcins, c=cohsnrs, s=cohsnrs, alpha=0.5)
hb = ax1.hexbin(cohO, bcins, gridsize=75, C=[isv[detectors[-1]] for isv in injsnrs], alpha=0.9)

# set color bar
axcb = pl.subplot(gs[1])
cb = pl.colorbar(hb, cax=axcb)
cb.set_label('signal-to-noise ratio')

xlims = ax1.get_xlim()
ylims = ax1.get_ylim()

ax1.plot(icohO, ibcins, 'ko', ms=5, alpha=0.3)
ax1.set_ylim(ylims)

axty = ax1.twiny()
axtx = ax1.twinx()

ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
ax1.xaxis.grid(False)
ax1.yaxis.grid(False)

#axty.set_xlim(10**np.array(xlims))
#axty.set_xscale('log', nonposx='clip')
#logpos = np.arange(np.around(xlims[0], decimals=-1)+10, np.around(xlims[1], decimals=-1)+10., 10)
#axty.set_xticks(np.logspace(logpos[0], logpos[-1], len(logpos)))
axtx.set_ylim(10**np.array(ylims))
axtx.set_yscale('log', nonposx='clip')

axty.xaxis.set_ticks_position('bottom')
axty.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='off')
axtx.yaxis.set_ticks_position('left')

#axty.set_xlabel(r'$\mathcal{O}_{{\rm S}/{\rm N}}$')
axtx.set_ylabel(r'$\mathcal{O}_{{\rm S}/{\rm I}}$')
axty.xaxis.set_label_position('bottom')
axtx.yaxis.set_label_position('left')

# on right plot show both coherent and incoherent over full ranges
ax2 = pl.subplot(gs[2]) #, sharex=axty)
ax2.loglog(10**cohO, 10**bcins, 'ro', ms=5, alpha=0.8, label='Coherent')

# overplot incoherent signal values
ax2.loglog(10**icohO, 10**ibcins, 'ko', ms=5, alpha=0.3, label='Incoherent')

ax2.set_xlabel(r'$\mathcal{O}_{{\rm S}/{\rm N}}$')
ax2.set_ylabel(r'$\mathcal{O}_{{\rm S}/{\rm I}}$')

ax2.legend(loc='best')

# re-set top plot x-axis to match bottom plot
xlims = ax2.get_xlim()
axty.set_xlim(xlims)
axty.set_xscale('log', nonposx='clip')
logpos = np.arange(np.around(np.log10(xlims[0]), decimals=-1)+10, np.around(np.log10(xlims[1]), decimals=-1), 10)
axty.set_xticks(np.logspace(logpos[0], logpos[-1], len(logpos)))
ax2.set_xticks(np.logspace(logpos[0], logpos[-1], len(logpos)))
ax1.set_xlim(np.log10(xlims))

#fig.tight_layout()

# output odds plots
ppdir = os.path.join(os.getcwd(), 'odds')
if not os.path.isdir(ppdir):
  os.makedirs(ppdir)
fig.savefig(os.path.join(ppdir, 'odds_plot.png'), dpi=300)
fig.savefig(os.path.join(ppdir, 'odds_plot.pdf'))

p = sp.Popen('pdftops -eps %s' % os.path.join(ppdir, 'odds_plot.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:oddsplot}
The coherent signal versus incoherent signal {\it or} noise odds ($\mathcal{O}_{{\rm S/I}_{\rm simple}}$ given
in Equation~\ref{eq:cohvincoh2}) plotted against the coherent signal versus noise odds ($\mathcal{O}_{\rm S/N}$ given
in Equation~\ref{eq:sigvsnoise}), for simulated signals injected as if into the two LIGO Hanford and Livingston
detectors. Half the simulations were created with coherent signals between the two detectors (hexagonal-binned histogram
points in the upper plots, and red points in the lower plot), whilst the other half where incoherent between the
two detectors (black points in both plots). The upper plot is a zoomed in version of the lower plot, with hexagonal
histogram bins coloured to show the coherent signal SNR.
"""

fp = open(os.path.join(ppdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()

# plot SNR vs bcins
ppdir = os.path.join(os.getcwd(), 'snr_vs_odds')
pl.clf()
fig = pl.figure(figsize=(6,5))
ax = pl.gca()

ax.plot(cohsnrs, bcis, 'o', ms=3, color='orange', alpha=0.1, label=r'$\mathcal{O}_{{\rm S/I}_{\rm simple}}$')
ax.plot(cohsnrs, bcins, '*', ms=3, color='blue', alpha=0.1, label=r'$\mathcal{O}_{\rm S/I}$')
ax.plot(cohsnrs, cohO, 'o', ms=3, color='green', alpha=0.1, label=r'$\mathcal{O}_{\rm S/N}$')
ax.set_ylim([-5., 12.])

# get fits to the log(O^S_N) distribution as a function of SNR^2
f = np.polyfit(cohsnrs, cohO, 2)
csnrs = np.linspace(0.01, 20., 50)
ax.plot(csnrs, f[2] + f[1]*csnrs + f[0]*csnrs**2, 'g-', lw=2)

# get fits to the max of O^S_Isimple (note: not the log) as a function of SNR
dsnr = 0.1 # steps in SNR
nsnrs = np.arange(1.5, 20, dsnr)
maxbcis = np.zeros(len(nsnrs))
snrvals = np.zeros(len(nsnrs))
for i, nsnr in enumerate(nsnrs):
  isnr = np.argmax(bcis[(cohsnrs >= nsnr) & (cohsnrs < nsnr+dsnr)])
  maxbcis[i] = bcis[(cohsnrs >= nsnr) & (cohsnrs < nsnr+dsnr)][isnr]
  snrvals[i] = cohsnrs[(cohsnrs >= nsnr) & (cohsnrs < nsnr+dsnr)][isnr]

f = np.polyfit(np.log10(snrvals), maxbcis, 1)
ax.plot(csnrs[csnrs>1.5], f[1] + np.log10(csnrs[csnrs>1.5])*f[0], '-', color='orange', lw=2)

ax.set_xlabel(r'$\rho_{\rm coh}$')
ax.set_ylabel(r'$\mathcal{O}$')
ax.legend(loc='best')
ax.set_xlim([0., 20.])

fig.canvas.draw()
# change y labels to be of the form 10^{N}
ylabels = []
for ticklabel in ax.get_yticklabels():
  xl = ticklabel.get_text()
  xl = '{' + xl.replace("$", "") + '}'
  ylabels.append(r'$10^{%s}$' % xl)
ax.set_yticklabels(ylabels)

fig.tight_layout()

# output snr plots
if not os.path.isdir(ppdir):
  os.makedirs(ppdir)
fig.savefig(os.path.join(ppdir, 'snr_v_odds_plot.png'), dpi=300)
fig.savefig(os.path.join(ppdir, 'snr_v_odds_plot.pdf'))

p = sp.Popen('pdftops -eps %s' % os.path.join(ppdir, 'snr_v_odds_plot.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:snrvsodds}
Various odds values plotted as a function of the injected coherent signal-to-noise ratio ($\rho_{\text coh}$).
Shown are the (simple) coherent signal versus incoherent signal odds ($\mathcal{O}_{{\rm S/I}_{\rm simple}}$ given
in Equation~\ref{eq:cohvincoh1}), more complete coherent signal versus incoherent signal {\it or} noise odds
($\mathcal{O}_{\rm S/I}$ given in Equation~\ref{eq:cohvincoh2}), and coherent signal versus Gaussian noise odds
($\mathcal{O}_{\rm S/N}$ given in Equation~\ref{eq:sigvsnoise}). Also shown are lines with a quadratic fit
to the $\mathcal{O}_{\rm S/N}$ values, to show the $\log{\mathcal{O}_{\rm S/N}} \propto \rho_{\text{coh}}^2$ relation,
and a fit to show the $\mathcal{O}_{{\rm S/I}_{\rm simple}} \propto \rho_{\text{coh}}$ relation.
"""

fp = open(os.path.join(ppdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()
