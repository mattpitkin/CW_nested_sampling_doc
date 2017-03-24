#!/usr/bin/env python

"""
Script to generate plots of upper limit comparisons, num. posterior sample vs live points, evidence comparisons etc
"""

from __future__ import division

import os
import sys
import json
import numpy as np
import tarfile
import subprocess as sp
import collections

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

# tarball containing data
tarball = 'upper_limits.tar.gz'
tar = tarfile.open(tarball, "r:gz")

Nlives = ['256', '512', '1024', '2048', '4096', '8192', '16384']

allh0s_nest = collections.OrderedDict()
allh0s_grid = collections.OrderedDict()
allh0s_diff = collections.OrderedDict()

allevs_nest = collections.OrderedDict()
allevs_grid = collections.OrderedDict()
allevs_diff = collections.OrderedDict()
allevs_err = collections.OrderedDict()

allnposts = collections.OrderedDict()
allinfos = collections.OrderedDict()

nli = []

for nl in Nlives:
    nli.append(int(nl))
    allh0s_nest[nl] = []
    allh0s_grid[nl] = []
    allh0s_diff[nl] = []
    
    allevs_nest[nl] = []
    allevs_grid[nl] = []
    allevs_diff[nl] = []
    allevs_err[nl] = []
    
    allnposts[nl] = []
    allinfos[nl] = []

for m in tar.getmembers():
  if m.isfile():
    if 'data.json' in m.name:
      fp = tar.extractfile(m)
      d = json.load(fp)

      nlive = str(d['information']['nlive'])

      # injected and recovered SNRs
      h0g = d['h0uls']['grid']
      h0n = d['h0uls']['nested']
      allh0s_grid[nlive].append(h0g)
      allh0s_nest[nlive].append(h0n)
      allh0s_diff[nlive].append(100.*(h0n-h0g)/h0g)
      
      evg = d['evrats']['grid']
      evn = d['evrats']['nested']
      allevs_grid[nlive].append(evg)
      allevs_nest[nlive].append(evn)
      allevs_diff[nlive].append(100.*(np.exp(evn-evg)-1.))
      allevs_err[nlive].append(d['information']['everr'])
      allinfos[nlive].append(d['information']['information'])

      allnposts[nlive].append(d['posterior']['num_samples'])

tar.close() # close tarball

# plot number of posterior samples as a function of number of live points
outdir = 'numposts'

fig = pl.figure(figsize=(7,5))
ax = pl.gca()

npmeds = []
npmeans = [] 
for npvals in allnposts.values():
    npmeds.append(np.median(npvals))
    npmeans.append(np.mean(npvals))

f = np.polyfit(np.log2(nli), np.log2(npmeans), 1)

ax.violinplot(allnposts.values(), np.log2(nli), showextrema=False, showmeans=True)
ax.text(8, 70000, r'$\langle N_{\rm post} \rangle \approx %.1f N_{\rm live}$' % (2**f[1])) 
ax.set_xlabel(r'No. live points $(N_{\rm live})$')
ax.set_ylabel(r'No. posterior samples $(N_{\rm post})$')

ax.set_xticks(np.log2(nli))
ax.set_xticklabels([r'$2^{%d}$' % int(nl) for nl in np.log2(nli)])

fig.tight_layout()

fig.savefig(os.path.join(outdir, 'numposts.png'), dpi=300)
fig.savefig(os.path.join(outdir, 'numposts.pdf'))
p = sp.Popen('pdftops -eps %s' % os.path.join(outdir, 'numposts.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:numposts}
The number of posterior samples produced as a function of the number of live points used in the nested sampling algorithm.
"""

fp = open(os.path.join(outdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()

# clear figure
pl.clf()

# plot the distribution of signal versus noise odds as a function of number of live points used
fig = pl.figure(figsize=(7,5))
ax = pl.gca()

outdir = 'nest_evs'

evmeans = []
everrs = []
evstds = []
infos = []
for evs, evstrue, errs, info in zip(allevs_nest.values(), allevs_grid.values(), allevs_err.values(), allinfos.values()):
    evmeans.append(np.mean(np.array(evs)-np.array(evstrue))/np.log(10.))
    everrs.append(np.mean(errs)/np.log(10.))
    evstds.append(np.std(np.array(evs)-np.array(evstrue))/np.log(10.))
    infos.append(np.mean(info))

print infos

ax.violinplot((np.array(allevs_nest.values())-np.array(allevs_grid.values())).T/np.log(10.), np.log2(nli), showextrema=False, showmeans=True)
ax.errorbar(np.log2(nli), evmeans, yerr=everrs, fmt='o', capsize=3, capthick=1, color='r', markersize=3, label=r'$\sigma = \sqrt{H/N_{\rm live}}$')
a, c, b = ax.errorbar(np.log2(nli), evmeans, yerr=evstds, fmt='o', capsize=3, capthick=1, color='g', markersize=3, elinewidth=1, label=r'$\sigma$ distribution')
b[0].set_linestyle('--')
ax.legend(loc='upper right')
ax.set_ylabel(r'$\log{}_{10}\left(\mathcal{Z}_{\rm nest}/\mathcal{Z}_{\rm grid}\right)$')
ax.set_xlabel(r'No. live points $(N_{\rm live})$')

ax.set_xticks(np.log2(nli))
ax.set_xticklabels([r'$2^{%d}$' % int(nl) for nl in np.log2(nli)])

# show equivalent percentage evidence (not log evidence) offsets
axnew = ax.twinx()
DZmin, DZmax = ax.get_ybound()
pDZvals = np.arange(np.around(100.*(10**(DZmin)-1.), decimals=-1), np.around(100.*(10**(DZmax)-1.), decimals=-1), 10)
equivlnZ = np.log10((pDZvals/100.)+1.)

axnew.grid(b=False)
axnew.set_yticks(equivlnZ)
axnew.set_ybound(ax.get_ybound())

axnew.get_yaxis().set_tick_params(which='both', direction='out')
axnew.set_yticklabels(['$%d$' % yv for yv in pDZvals])
axnew.set_ylabel(r'$(\mathcal{Z}_{\rm nest}-\mathcal{Z}_{\rm grid})/\mathcal{Z}_{\rm grid} \%$')

fig.tight_layout()

fig.savefig(os.path.join(outdir, 'nest_evs.png'), dpi=300)
fig.savefig(os.path.join(outdir, 'nest_evs.pdf'))
p = sp.Popen('pdftops -eps %s' % os.path.join(outdir, 'nest_evs.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:nest_evs}
The distributions of ratios of evidence values calculated as a function of the number of live points used, based on searches over $h_0$, $\cos{\iota}$,
$\psi$ and $\phi0$ on 500 simulated Gaussian noise realisations for each $N_{\text{live}}$ value. The ratio compares the evidence produced by the
nested sampling algorithm in \lppen to that produced from a grid over the parameter space with \lppe (which we assume is a representation of the
true value). A comparison of the estimated error on the evidence value using the information gain, $H$ (which is $\sim 2.4$ here), with the measured
standard deviation of the distributions is also plotted and agree rather well.
"""

fp = open(os.path.join(outdir, 'caption.tex'), 'w')
fp.write(caption)
fp.close()

# clear figure
pl.clf()

# plot distributions of 95% h0 upper upper_limit
fig = pl.figure(figsize=(7,5))
ax = pl.gca()

outdir = 'uls'
    
r = ax.violinplot(allh0s_diff.values(), np.log2(nli), showextrema=True, showmeans=True)
ax.set_xlabel(r'No. live points $(N_{\rm live})$')
ax.set_ylabel(r'$\left(\{h_0^{95\%}\}_{\rm nest} - \{h_0^{95\%}\}_{\rm grid}\right)/\{h_0^{95\%}\}_{\rm grid} \%$')

yv = np.max([np.max(np.abs(hd)) for hd in allh0s_diff.values()])
ax.set_ylim([-1.1*yv, 1.1*yv])

r['cmaxes'].set_linewidth(0.5)
r['cmaxes'].set_color('darkmagenta')
r['cmins'].set_linewidth(0.5)
r['cmins'].set_color('darkmagenta')
r['cbars'].set_linewidth(1.0)
r['cbars'].set_color('darkmagenta')

fig.tight_layout()

fig.savefig(os.path.join(outdir, 'uls.png'), dpi=300)
fig.savefig(os.path.join(outdir, 'uls.pdf'))
p = sp.Popen('pdftops -eps %s' % os.path.join(outdir, 'uls.pdf'), shell=True)
p.communicate()

caption = r"""\label{fig:uls}
The distributions of ratios 95\% upper limits on $h_0$ as a function of the number of live points used, based on searches over $h_0$, $\cos{\iota}$,
$\psi$ and $\phi0$ on 500 simulated Gaussian noise realisations for each $N_{\text{live}}$ value. The ratio compares the upper limits produced from the
output of the nested sampling algorithm in \lppen to that produced from a grid over the parameter space with \lppe (which we assume is a representation of the
true value).
"""
