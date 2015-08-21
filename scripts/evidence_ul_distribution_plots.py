#!/usr/bin/env python

"""
Script to produce plots from the output of evidence_ul_distribution.py
"""

from matplotlib import pyplot as pl
from matplotlib.patches import Polygon

import numpy as np
import json
import os
import subprocess as sp
import gzip

from lalapps import pulsarpputils as pppu

# function for filling in boxplots
def fillboxes(axis, boxplot, data):
  numBoxes = len(data)
  medians = range(numBoxes)
  for i in range(numBoxes):
    box = boxplot['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
      boxX.append(box.get_xdata()[j])
      boxY.append(box.get_ydata()[j])
    boxCoords = zip(boxX,boxY)

    boxPolygon = Polygon(boxCoords, facecolor='grey')
    axis.add_patch(boxPolygon)

    med = boxplot['medians'][i]
    # Finally, overplot the sample averages, with horizontal alignment
    # in the center of each box
    axis.plot([np.average(med.get_xdata())], [np.average(data[i])], color='w', marker='*', markeredgecolor='k', ms=8)

# set some plot defaults
pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.rc('font', size=14)

basedir = '/home/sismp2/projects/code_testing/evidence_ul_distribution'

nlives = [128, 256, 512, 1024, 2048, 4096]

# get values from grid-based evaluation
griddata = os.path.join(basedir, 'gridoutput.txt')
fp = open(griddata, 'r')
info = json.load(fp)
fp.close()

# produce a boxplot of the number of live points versus the mean and standard deviation of the distribution 
# of evidence ratio (also plot the grid-based value of the evidence ratio)
datast = []
datagauss = []
for n in nlives:
  # get the evidence ratios from the '_B' files
  p = (sp.check_output("cat "+os.path.join(basedir, "%d/gaussian/post_*.txt_B.txt" % n), shell=True)).split('\n')
  TAodds = np.array([float(line.split()[0]) for line in p if len(line) > 0])
  datagauss.append(TAodds)
  
  p = (sp.check_output("cat "+os.path.join(basedir, "%d/studentst/post_*.txt_B.txt" % n), shell=True)).split('\n')
  TAodds = np.array([float(line.split()[0]) for line in p if len(line) > 0])
  datast.append(TAodds)

fig, ax = pl.subplots(figsize=(7,5))
bp = pl.boxplot(datast, whis=[5, 95], notch=0, sym='x', positions=np.log2(nlives)) 
pl.setp(bp['fliers'], color='blue')
pl.setp(bp['boxes'], color='black')
pl.setp(bp['whiskers'], color='black')

# fill in the boxes
fillboxes(ax, bp, datast)

pl.xticks(np.log2(nlives).tolist(), nlives)
ax.set_xlim((np.log2(nlives[0]/2.), np.log2(nlives[-1]*2.)))
pl.axhline(y=info['Odds ratios']['studentst'], color='k', linewidth=2, ls='dashed')
pl.xlabel('Number of live points')
pl.ylabel('log(Odds Ratio)')
outfig = os.path.join(basedir, 'oddsratio_distr_studentst.pdf')
fig.savefig(outfig)
fig.clf()
pl.close(fig)

fig, ax = pl.subplots(figsize=(7,5))
bp = pl.boxplot(datagauss, whis=[5, 95], notch=0, sym='x', positions=np.log2(nlives)) 
pl.setp(bp['fliers'], color='blue')
pl.setp(bp['boxes'], color='black')
pl.setp(bp['whiskers'], color='black')

# fill in the boxes
fillboxes(ax, bp, datagauss)

pl.xticks(np.log2(nlives).tolist(), nlives)
ax.set_xlim((np.log2(nlives[0]/2.), np.log2(nlives[-1]*2.)))
pl.axhline(y=info['Odds ratios']['gaussian'], color='k', linewidth=2, ls='dashed')
pl.xlabel('Number of live points')
pl.ylabel('log(Odds Ratio)')
outfig = os.path.join(basedir, 'oddsratio_distr_gaussian.pdf')
fig.savefig(outfig)
fig.clf()
pl.close(fig)

# produce plots of the distribution of h0 upper limits as a function of number of live points
h95sst = []
h95sgauss = []
nsampsst = []
nsampsgauss = []

for n in nlives:
  h95ssttmp = []
  h95sgausstmp = []
  nsampssttmp = []
  nsampsgausstmp = []  

  for i in range(125):
    postfile = os.path.join(basedir, "%d/studentst/post_%04d.txt.gz" % (n, i))

    # check position of h0 value
    if i == 0:
      # read in header
      fp = gzip.open(postfile, 'r')
      header = (fp.readline()).split()
      header = [h.upper() for h in header]
      fp.close()

      # get required parameter columns
      usecols = (header.index('H0'),)
  
    # read in h0 values
    pos = np.loadtxt(postfile, skiprows=1, usecols=usecols) # skip header line

    # get the length of the posterior samples
    nsampssttmp.append(pos.shape[0])

    # get 95% upper limit
    #h95ssttmp.append(pppu.upper_limit(pos, upperlimit=0.95, parambounds=[0., np.inf], nbins=50))
    h95ssttmp.append(pppu.upper_limit_greedy(pos, upperlimit=0.95, nbins=200))

    # now for the Gaussian likelihood
    postfile = os.path.join(basedir, "%d/gaussian/post_%04d.txt.gz" % (n, i))

    pos = np.loadtxt(postfile, skiprows=1, usecols=usecols)
    nsampsgausstmp.append(pos.shape[0])
    #h95sgausstmp.append(pppu.upper_limit(pos, upperlimit=0.95, parambounds=[0., np.inf], nbins=50))
    h95sgausstmp.append(pppu.upper_limit_greedy(pos, upperlimit=0.95, nbins=200))  

  h95sst.append(np.array(h95ssttmp))
  h95sgauss.append(np.array(h95sgausstmp))
  nsampsst.append(np.array(nsampssttmp))
  nsampsgauss.append(np.array(nsampsgausstmp))

fig, ax = pl.subplots(figsize=(7,5))
ax2 = ax.twiny()
bp = ax.boxplot(h95sst, whis=[5, 95], notch=0, sym='x', positions=np.log2(nlives))
pl.setp(bp['fliers'], color='blue')
pl.setp(bp['boxes'], color='black')
pl.setp(bp['whiskers'], color='black')

# fill in the boxes
fillboxes(ax, bp, h95sst)

ax.set_xticks(np.log2(nlives))
ax.set_xticklabels(["%d" % n for n in nlives])
ax.set_xlim((np.log2(nlives[0]/2.), np.log2(nlives[-1]*2.)))
ax.axhline(y=info['Upper limits']['studentst'], color='k', linewidth=2, ls='dashed')
ax.set_xlabel('Number of live points')
ax.set_ylabel('upper limit $h_0^{95\%}$')

# add another x-axis label for the average number of posterior samples (rounded to the nearest 10)
ax2.set_xticks(np.log2(nlives))
ax2.set_xticklabels(["$%d \pm %d$" % (int(np.around(np.mean(na), decimals=-1)), int(np.around(np.std(na), decimals=-1))) for na in nsampsst], rotation=45)
ax2.set_xlabel("Number of posterior samples")

outfig = os.path.join(basedir, 'h95_distr_studentst.pdf')
fig.savefig(outfig)
fig.clf()
pl.close(fig)

fig, ax = pl.subplots(figsize=(7,5))
ax2 = ax.twiny()
bp = ax.boxplot(h95sgauss, whis=[5, 95], notch=0, sym='x', positions=np.log2(nlives))
pl.setp(bp['fliers'], color='blue')
pl.setp(bp['boxes'], color='black')
pl.setp(bp['whiskers'], color='black')

# fill in the boxes
fillboxes(ax, bp, h95sgauss)

ax.set_xticks(np.log2(nlives))
ax.set_xticklabels(["%d" % n for n in nlives])
ax.set_xlim((np.log2(nlives[0]/2.), np.log2(nlives[-1]*2.)))
ax.axhline(y=info['Upper limits']['gaussian'], color='k', linewidth=2, ls='dashed')
ax.set_xlabel('Number of live points')
ax.set_ylabel('upper limit $h_0^{95\%}$')

# add another x-axis label for the average number of posterior samples (rounded to the nearest 10)
ax2.set_xticks(np.log2(nlives))
ax2.set_xticklabels(["$%d \pm %d$" % (int(np.around(np.mean(na), decimals=-1)), int(np.around(np.std(na), decimals=-1))) for na in nsampsgauss], rotation=45)
ax2.set_xlabel("Number of posterior samples")

outfig = os.path.join(basedir, 'h95_distr_gaussian.pdf')
fig.savefig(outfig)
fig.clf()
pl.close(fig)

