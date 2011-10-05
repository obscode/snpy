'''This is a container class for all the light-curve generators
we have.  At present, there are 3:

   - Jose-Louis Prieto's template generator for bands Bs,Vs,Rs,Is
   - Burns 2008 generator for u,B,V,g,r,i,Y,J,H
   - Kevin Krisciunas's polynomial fits for J,H,K
   - Burns 2011 2nd generation generator for uBVgriYJH, based on stretch

selecting one of these filters as the restband in SNOOPY will
automatically select the appropriate generator.  We now have a couple new
keywords:  [gen], which corresponds to the generation of the lightcurve,
and [param] which can refer to 'dm15' or 'st'   

Got rid of the MMax function, as this is now part of the model class.'''
import CSPtemp as dm15temp2
import dm15temp

class template:

   def __init__(self):
      self.Pt = dm15temp.template()
      self.Ct = dm15temp2.dm15_template()

      for band in ['Bs','Vs','Rs','Is']:
         self.__dict__[band] = self.Pt.__dict__[band[0]]
         self.__dict__['e'+band] = self.Pt.__dict__['e'+band[0]]
      for band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         self.__dict__[band] = self.Ct.__dict__[band]
         self.__dict__['e'+band] = self.Ct.__dict__['e'+band]

   def __getattr__(self, name):
      if name in ['Bs','Vs','Rs','Is','eBs','eVs','eRs','eIs']:
         # Use Prieto values
         return self.Pt.__dict__[name[:-1]]
      elif name in ['u','B','V','g','r','i','Y','J','H','H_K','J_K',
                    'eu','eB','eV','eg', 'er','ei','eY','eJ','eH','eJ_K','eH_K']:
         return self.Ct.__dict__[name]
      elif name == 't':
         return self.Ct.t
      raise AttributeError, "Error:  attribute %s not defined" % (name)


   def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
      args1 = {'dm15':dm15, 'method':1, 'colors':'none','generate':generate}
      args2 = {'dm15':dm15, 'dm15_int':dm15_int, 'dm15_colors':'int',
               'generate':generate}
      self.Pt.mktemplate(**args1)
      self.Ct.mktemplate(**args2)

   def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True):
      if band in ['Bs','Vs','Rs','Is']:
         return(self.Pt.eval(band[0], times, z))
      elif band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         return(self.Ct.eval(band, times, z, mag, sextrap, gen, toff))
      else:
         raise AttributeError,"Sorry, band %s is not supported" % band

class stemplate:

   def __init__(self):
      self.Ct = dm15temp2.st_template()

      for band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         self.__dict__[band] = self.Ct.__dict__[band]
         self.__dict__['e'+band] = self.Ct.__dict__['e'+band]

   def __getattr__(self, name):
      if name in ['u','B','V','g','r','i','Y','J','H','H_K','J_K',
                    'eu','eB','eV','eg', 'er','ei','eY','eJ','eH','eJ_K','eH_K']:
         return self.Ct.__dict__[name]
      elif name == 't':
         return self.Ct.t
      raise AttributeError, "Error:  attribute %s not defined" % (name)


   def mktemplate(self, st, dm15_int=None, dm15_colors='int', generate=0):
      args = {'st':st, 'dm15_int':dm15_int, 'dm15_colors':'int',
               'generate':generate}
      self.Ct.mktemplate(**args)

   def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True):
      if band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         return(self.Ct.eval(band, times, z, mag, sextrap, gen, toff))
      else:
         raise AttributeError,"Sorry, band %s is not supported" % band

