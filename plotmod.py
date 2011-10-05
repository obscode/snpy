import os

# Check to see if the user has asked for a particular plotting module
if "PLOTMOD" in os.environ:
   if os.environ['PLOTMOD'] == 'PGPLOT':
      try:
         from plot_sne_pg import *
      except:
         raise ImportError, "Sorry, PGPLOT wrappers failed to load"
   elif os.environ['PLOTMOD'] == "MPL":
      #try:
      from plot_sne_mpl import *
      #except:
      #   raise ImportError, "Sorry, matplotlib wrappers failed to load"
   else:
      raise ImportError, "Sorry, unrecognized PLOTMOD '%s'" % \
            os.environ['PLOTMOD']
else:
   # Okay, try MPL first and if not, PYGPLOT
   try:
      from plot_sne_mpl import *
      version='mpl'
   except:
      try:
         from plot_sne_pg import *
         version='mpl'
      except:
         st = "\nSorry, you don't seem to have matplotlib or pygplot wrappers.\n"+\
               "You can get pygplot from http://astro.swarthmore.edu/~burns/pygplot/\n"+\
               "Or matplotlib from http://matplotlib.sourceforge.net/"
         raise ImportError, st

