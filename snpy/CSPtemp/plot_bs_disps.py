#!/usr/bin/env python

import sys,os,string
from numpy import *
from scipy.interpolate import *
from myplotlib import PanelPlot
from matplotlib import pyplot
import pickle

tck_file0 = 'tck.pickle'
tck_file1 = 'bs_tck.pickle'

f = open(tck_file0)
all_tck0 = pickle.load(f)
f.close()

f = open(tck_file1)
all_tcks = pickle.load(f)
f.close()

dm15_low = 0.6
dm15_high = 2.0
t_low = -10.0
t_high = 70

dm15s = [0.9, 1.1, 1.3, 1.5, 1.7, 1.9]
ts = arange(81)/80.0*(t_high - t_low) + t_low

#mp = PanelPlot(3,3)
mp2 = PanelPlot(3,3)

bands = ['u','B','V','g','r','i','Y','J','H']
for j in range(len(bands)):
   band = bands[j]
   tck0 = all_tck0[band]
   tcks = all_tcks[band]
   #mp.axes[j].text(0.5, 0.9, band, transform=mp.axes[j].transAxes,
   #      verticalalignment='top', horizontalalignment='center')
   mp2.axes[j].text(0.5, 0.9, band, transform=mp2.axes[j].transAxes,
         verticalalignment='top', horizontalalignment='center')
   for k in range(len(dm15s)):
      dzs = []
      z0 = bisplev(ts, dm15s[k], tck0)[:,0]
      for i in range(len(tcks)):
         z1 = bisplev(ts, dm15s[k], tcks[i])[:,0]
         if any(absolute(z1 - z0) > 0.15):
            continue
         dzs.append(z1 - z0)
         if k == 1:
            mp2.axes[j].plot(ts, dzs[-1], color='0.65')
      dzs = array(dzs)
      rms = sqrt(mean(power(dzs, 2), axis=0))
      mads = median(absolute(dzs), axis=0)
      #mp.axes[j].plot(ts, 1.49*mads, '-', label='%.1f' % (dm15s[k]))
      if k == 1:
         mp2.axes[j].plot(ts, rms, '-', color='red', linewidth=2)


#mp.axes[0].legend(prop={'size':8})
#mp.xlabel('$t - t_{max}(B)$ (days)')
mp2.xlabel('$t - t_{max}(B)$ (days)')
#mp.ylabel('Median Absolute Deviation')
mp2.ylabel('Master - Bootstrap')
#mp.set_limits(all_equal=1)
mp2.set_limits(all_equal=1)
#mp.draw()
mp2.draw()
pyplot.show()
#mp.close()
mp2.close()
