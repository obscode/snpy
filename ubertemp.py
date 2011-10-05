'''This is a container class for all the light-curve generators
we have.  At present, there are 3:

   - Jose-Louis Prieto's template generator for bands Bs,Vs,Rs,Is
   - Burns 2008 generator for u,B,V,g,r,i,Y,J,H
   - Kevin Krisciunas's polynomial fits for J,H,K

selecting one of these filters as the restband in SNOOPY will
automatically select the appropriate generator.'''
import CSPtemp as dm15temp2
import dm15temp
try:
   import dm15temp_dest
   have_dest = 1
except ImportError:
   have_dest = 0

dest_filts = ['d%d'%i for i in range(15)]
dest_filts = dest_filts + ['dw%d'%i for i in range(10)]

class template:

   def __init__(self):
      self.Pt = dm15temp.template()
      self.Ct = dm15temp2.template()
      if have_dest:  self.Dt = dm15temp_dest.template()


      for band in ['Bs','Vs','Rs','Is']:
         self.__dict__[band] = self.Pt.__dict__[band[0]]
         self.__dict__['e'+band] = self.Pt.__dict__['e'+band[0]]
      for band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         self.__dict__[band] = self.Ct.__dict__[band]
         self.__dict__['e'+band] = self.Ct.__dict__['e'+band]
      if have_dest:
         for i in range(15):
            self.__dict__['d%d'%i] = self.Dt.__dict__['d%d'%i]
            self.__dict__['ed%d'%i] = self.Dt.__dict__['ed%d'%i]
         for i in range(10):
            self.__dict__['dw%d'%i] = self.Dt.__dict__['dw%d'%i]
            self.__dict__['edw%d'%i] = self.Dt.__dict__['edw%d'%i]


   def __getattr__(self, name):
      if name in ['Bs','Vs','Rs','Is','eBs','eVs','eRs','eIs']:
         # Use Prieto values
         return self.Pt.__dict__[name[:-1]]
      elif name in ['u','B','V','g','r','i','Y','J','H','H_K','J_K',
                    'eu','eB','eV','eg', 'er','ei','eY','eJ','eH','eJ_K','eH_K']:
         return self.Ct.__dict__[name]
      elif name == 't':
         return self.Ct.t
      elif name in dest_filts and have_dest:
         return self.Dt.__dict__[name]
      raise AttributeError, "Error:  attribute %s not defined" % (name)


   def __setattr__(self, name, value):
      if name == 'Rv':
         self.Ct.Rv = value
      else:
         self.__dict__[name] = value

   def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
      args1 = {'dm15':dm15, 'method':1, 'colors':'none','generate':generate}
      args2 = {'dm15':dm15, 'dm15_int':dm15_int, 'dm15_colors':'int',
               'generate':generate}
      self.Pt.mktemplate(**args1)
      self.Ct.mktemplate(**args2)
      if have_dest:  self.Dt.mktemplate(**args2)

   def eval(self, band, times, z=0):
      if band in ['Bs','Vs','Rs','Is']:
         return(self.Pt.eval(band[0], times, z))
      elif band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         return(self.Ct.eval(band, times, z))
      elif band in dest_filts and have_dest:
         return(self.Dt.eval(band, times, z))
      else:
         raise AttributeError,"Sorry, band %s is not supported" % band

   def MMax(self, band):
      if band in ['Bs','Vs','Rs','Is']:
         return(self.Pt.MMax(band[0]))
      elif band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         return(self.Ct.MMax(band))
      elif band in dest_filts and have_dest:
         return(self.Dt.MMax(band))
      else:
         raise AttributeError,"Sorry, band %s is not supported" % band

