'''
Module for SNooPy to download/parse data from the Open Supernova Catalog.
'''
import urllib.request, urllib.parse, urllib.error
from snpy import lc
from numpy import array

SWIFT_URL = '''http://people.physics.tamu.edu/pbrown/SwiftSN/{}_uvotB15.1.dat'''

def load_lcs(obj):

   try:
      u = urllib.request.urlopen(SWIFT_URL.format(obj.name))
   except:
      raise ValueError("Object not found in SWIFT database")
   
   lines = u.readlines()
   lines = [line.strip().split() for line in lines if line[0] != '#']
   lines = [line for line in lines if len(line) == 12]

   data = {}
   for line in lines:
      if line[0] not in data:  data[line[0]] = []
      if line[2] == 'NULL' or line[3] == 'NULL': continue
      data[line[0]].append([float(line[1]), float(line[2]), float(line[3])])
   for band in data:
      data[band] = array(data[band])

   for band in ['UVW1','UVW2','UVM2']:
      if band not in data: continue
      obj.data[band] = lc(obj, band, data[band][:,0], data[band][:,1],
            data[band][:,2])
