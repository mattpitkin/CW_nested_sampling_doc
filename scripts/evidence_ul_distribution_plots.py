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

fig, ax = pl.subplots(figsize=(8,5))
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

fig, ax = pl.subplots(figsize=(8,5))
bp = pl.boxplot(datagauss, whis=[5, 95], notch=0, sym='x', positions=np.log2(nlives)) 
pl.setp(bp['fliers'], color='blue')
pl.setp(bp['boxes'], color='black')
pl.setp(bp['whiskers'], color='black')

# fill in the boxes
fillboxes(ax, bp, datagauss)

pl.xticks(np.log2(nlives).tolist(), nlives)
ax.set_xlim((np.log2(nlives[0]/2.), np.log2(nlives[-1]*2.)))
pl.axhline(y=info['Odds ratios']['gaussian'], color='k', linewidth=2)
pl.xlabel('Number of live points')
pl.ylabel('log(Odds Ratio)')
outfig = os.path.join(basedir, 'oddsratio_distr_gaussian.pdf')
fig.savefig(outfig)
fig.clf()
pl.close(fig)

# produce plots of the distribution of h0 upper limits as a function of number of posterior samples (which
# should be proportional to the number of live points)
