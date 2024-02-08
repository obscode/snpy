from __future__ import print_function
import os,sys
import numpy
import getpass

try:
   import pymysql as sql
   have_sql = 1
except:
   have_sql = 0

class sqlbase:
   '''A base class for an sql connection object.  Custom SQL query objects
   should inherit from this base class and override the member variables and
   read_sql() and update_sql() functions if necessary.'''

   host = "hostname"       # Host name of SQL server
   user = "username"       # user name on SQL server
   passwd = None           # password for username
   port = 3306             # SQL server port on hostname
   readonly = 0            # Can we read only?

   PHOTO_DB = "database"    # Database containing photometry
   SPEC_DB = None          # Database containing spectroscopy

   SN_TABLE = "SNe"        # SQL Table of supernova attributes
   SN_ID = "name"          # Unique SN ID field in SN_TABLE
   SN_ID2 = None           # Possible seoncary name (ugh!)
   SN_COND = ""            # Any extra SN query conditions
   PHOTO_TABLE = "Photo"   # SQL Table of photometry
   PHOTO_ID = "name"       # Unique SN ID field in PHOTO_TABLE
   PHOTO_FILT = "filter"   # SQL field for filter name in PHOTO_TABLE
   SPEC_TABLE = None       # SQL Table of spectroscopy
   SPEC_INFO = None        # SQL Table of spectal info
   SPEC_JD = "JD"          # SQL field in SPEC_INFO that gives JD
   SPEC_LAMB = "LAMBDA"    # SQL fieldin SPEC_TABLE that gives wavelength
   SPEC_FLUX = "FLUX"      # SQL field in SPEC_TABLE that gives flux
   SPEC_INDEX = "FILE"     # SQL index common to SPEC_INFO and SPEC_TABLE for JOINS
   SPEC_NAME = "SN"  


   # dictonary of SQL-sn object attributes ATTR_KEYS[attr] refers to the column
   # name in the SQL table SN_TABLE that corresponds to attribute attr in the
   # sn object.  Use None to indicate no information available
   ATTR_KEYS = {}
   PHOTO_ATTR_KEYS = {}

   FILTER_KEYS = {}        # Dictionary that maps SQL filters to SNpy filters
   PHOTO_JD = "JD"        # SQL field in PHOTO_TABLE that gives time (JD)
   JD_OFFSET = -2400000.5   #  An offset added to the JD before returning
   PHOTO_MAG = "m"       # SQL field in PHOTO_TABLE that gives magnitude
   PHOTO_EMAG = "e_m"    # SQL field in PHOTO_TABLE that gives error
   PHOTO_SNR = None      # SQL field that gives the signal-to-noise ratio
   PHOTO_K = "K"           # SQL field in PHOTO_TABLE that gives K-corrections
   PHOTO_FILT = "filter"   # SQL field in PHOTO_TABLE that gives filter name
   PHOTO_COND = ""         # Any extra WHERE conditions
   
   def __init__(self):
      self.connected=0

      # Allow for environment variable override
      for att in ['user','passwd','photo_db','spec_db','host','port']:
         key = 'SQL_' + att.upper()
         if key in os.environ:
            self.__dict__[att] = os.environ[key]

   def __getstate__(self):
      # Need to do this because the sql connection instances
      # themselves are not pickleable.
      dict = self.__dict__.copy()
      if 'con' in dict:  del dict['con']
      if 'c' in dict:  del dict['c']
      if 'c2' in dict:  del dict['c2']
      if 'con2' in dict:  del dict['con2']
      dict['connected'] = 0
      return dict

   def filter_key(self, name, epoch=None):
      '''Function to map an SQL filter string to a filter name.'''
      if name in self.FILTER_KEYS:
         return self.FILTER_KEYS[name]
      else:
         return name[0]

   def connect(self, name):
      '''Connect to the database.'''
      if self.passwd is None:
         passwd = getpass.getpass(prompt='SQL passwd for %s@%s:\n' % (self.user,
            self.host))
      else:
         passwd = self.passwd
      if not self.connected:
         try:
            self.con = sql.connect(host=self.host,user=self.user, passwd=passwd,
                  db=self.PHOTO_DB, port=self.port)
            if self.SPEC_DB is not None:
               if self.SPEC_DB != self.PHOTO_DB:
                  self.con2 = sql.connect(host=self.host, user=self.user, 
                        passwd=passwd, db=self.SPEC_DB)
               else:
                  self.con2 = self.con
            else:
               self.con2 = None
         except:
            raise
         else:
            # Cache the passwd for future use if it worked
            if self.passwd is None:  self.passwd = passwd
            
         self.c = self.con.cursor()
         if self.con2 is not None:
            self.c2 = self.con2.cursor()
         else:
            self.c2 = self.c
         self.connected = 1
         # Collect info about the SN_TABLE and PHOTO_TABLE:
         self.SN_table_info = self.get_table_info(self.SN_TABLE)
         self.PHOTO_table_info = self.get_table_info(self.PHOTO_TABLE)

      self.name = name

      # See if the SN object exists
      if self.SN_ID2 is not None:
         slct = '''SELECT %s,%s from %s where %s = %%s %s''' % \
             (self.SN_ID,self.SN_ID2,self.SN_TABLE, self.SN_ID, self.SN_COND)
      else:
         slct = '''SELECT %s from %s where %s = %%s %s''' % \
             (self.SN_ID,self.SN_TABLE, self.SN_ID, self.SN_COND)
      N = self.c.execute(slct, (self.name))
      if N == 0:
         return 0
      if N > 1:
         print("Warning!  %s is not unique in the database, taking first" %\
               (self.name))
         return N
      if self.SN_ID2 is not None:
         l = self.c.fetchall()[0]
         self.name2 = l[1]
      else:
         self.name2 = None
      return 1

   def get_table_info(self, table):
      '''Gets info from the table and returns it as a dictionary of dictionaries
      indexed by field name:
         {'field':{'type': field type,
                   'null': boolean,
                   'key': (PRI,MUL,UNI),
                   'default': default value,
                   'extra': extra info}, ..., } '''
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")
      if type(table) is not type([]):
         tables = [table]
      else:
         tables = table

      data = {}
      for table in tables: 
         data[table] = {}
         self.c.execute('''DESCRIBE %s''' % table)
         l = self.c.fetchall()
         for item in l:
            data[table][item[0]] = {'type':item[1],
                             'null':item[2],
                             'key':item[3],
                             'default':item[4],
                             'extra':item[5]}
      if len(tables) == 1:
         data = data[tables[0]]
      return(data)

   def close(self):
      '''Close connection to the database.'''
      if not self.connected: return
      self.con.close()
      if self.con2 is not None and self.con2 != self.con:
         self.con2.close()
      self.connected = 0

   def sql_field(self, attr):
      '''Use self.ATTR_KEYS to convert an SN object attribute into an
      SQL field.'''
      if attr in self.ATTR_KEYS:
         return self.ATTR_KEYS[attr]
      else:
         return attr

   def create_SN(self, ra, decl, z):
      '''Create a new SN with values for ra, decl, and z'''
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")
      insrt = '''INSERT INTO %s (%s,%s,%s,%s) VALUES (%%s,%%s,%%s,%%s)''' % \
            (self.SN_TABLE, self.SN_ID, self.sql_field('ra'), self.sql_field('decl'), 
                  self.sql_field('z'))
      self.c.execute(insrt, (self.name, ra, decl, z))

   def create_SN_photometry(self, data):
      '''Insert the photomtery data into the database.  See get_SN_photometry
      for an explanation of data'''
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")
      filters = list(data.keys())
      if type(PHOTO_TABLE) is type([]):
         raise TypeError("function not supported for table joins")
      insrt = '''INSERT INTO %s (%s,%s,%s,%s,%s,%s) VALUES (%%s,%%s,%%s,%%s,%%s,%%s)''' %\
            (self.PHOTO_TABLE,self.PHOTO_ID,self.PHOTO_FILT,self.PHOTO_JD,
             self.PHOTO_MAG, self.PHOTO_EMAG, self.PHOTO_K)
      for filter in data:
         d = data[filter]
         for i in range(len(d[0])):
            if d[3] is None:
               self.c.execute(insrt, (self.name,filter,d[0][i],d[1][i],
                  d[2][i],0.0))
            else:
               self.c.execute(insrt, (self.name,filter,d[0][i],d[1][i],
                  d[2][i],d[3][i]))


   def delete_SN(self, name):
      '''Delete SN from SN_TABLE.'''
      delete = '''DELETE from %s where %s = %%s''' %\
            (self.SN_TABLE, self.SN_ID)
      self.c.execute(delete, (self.name))

   def get_SN_parameter(self, attr, verbose=0):
      '''Gets a single attribute of a SN in the database.'''
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")

      if attr not in self.SN_table_info and attr not in self.ATTR_KEYS:
         raise AttributeError("SN table does not have this attribute: %s" % \
               attr)
      field = self.sql_field(attr)

      slct = '''SELECT %s from %s where %s = %%s %s''' % (field, self.SN_TABLE, 
            self.SN_ID, self.SN_COND)
      if verbose:
         print("executing...  ",slct % (self.name))
      N = self.c.execute(slct, (self.name))
      if N == 0:
         raise AttributeError("Attribute %s(%s) not found in SQL db" % \
               (attr,field))
      return(self.c.fetchall()[0][0])

   def set_SN_parameter(self, attr, value):
      '''Sets a single attribute of a SN in the datase.'''
      if self.readonly:
         raise ValueError("Database is read-only")
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")
      if attr not in list(self.SN_table_info.keys()) and attr not in self.ATTR_KEYS:
         raise AttributeError("SN table does not have this attribute: %s" % \
               attr)
      field = self.sql_field(attr)
      if field not in self.SN_table_info:
         raise AttributeError("SN table does not have this attribute: %s" % \
               field)

      updt = '''UPDATE %s SET %s = %%s where %s = %%s''' % \
            (self.SN_TABLE, field, self.SN_ID)
      self.c.execute(updt, (value, self.name))

   def get_SN_photometry(self, verbose=0):
      '''Get the photometry form the SQL database.  Returns a dictionary
      indexed by filter name.  Each element is a 5-tuple of arrays:
      (MJD, mag, e_mag, K).'''
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")
      if self.name2 is not None:
         name_where = '(%s="%s" or %s="%s")' % (self.PHOTO_ID,self.name,
                                                self.PHOTO_ID,self.name2)
      else:
         name_where = '%s="%s"' % (self.PHOTO_ID,self.name)

      if type(self.PHOTO_TABLE) is type([]):
         phot_table = ','.join(self.PHOTO_TABLE)
      else:
         phot_table = self.PHOTO_TABLE

      slct = '''SELECT %s,%s,%s,%s''' % \
            (self.PHOTO_FILT, self.PHOTO_JD, self.PHOTO_MAG, self.PHOTO_EMAG)

      if self.PHOTO_K is not None:
         slct += ",%s" % (self.PHOTO_K)
      if self.PHOTO_SNR is not None:
         slct += ",%s" % (self.PHOTO_SNR)

      slct += " from %s where %s %s ORDER by %s" % \
            (phot_table, name_where, self.PHOTO_COND, self.PHOTO_JD)
      if verbose:  print("executing... ",slct)
      N = self.c.execute(slct)
      if N == 0:
         raise ValueError("No photometry for %s found" % \
                  (self.name))
      list = self.c.fetchall()

      data = {}
      for l in list:
         # We now make this a function to handle more complicated cases
         filter = self.filter_key(l[0], l[1] + self.JD_OFFSET)
         if filter in data:
            data[filter]['t'].append(l[1] + self.JD_OFFSET)
            data[filter]['m'].append(l[2])
            data[filter]['em'].append(l[3])
            ii = 0
            if self.PHOTO_K is not None:
               data[filter]['K'].append(l[4])
               ii += 1
            if self.PHOTO_SNR is not None:
               data[filter]['SNR'].append(l[4+ii])
         else:
            data[filter] = {'t':[l[1] + self.JD_OFFSET],
                            'm':[l[2]], 'em':[l[3]]}
            ii = 0
            if self.PHOTO_K is not None:
               data[filter]['K'] = [l[4]]
            if self.PHOTO_SNR is not None:
               data[filter]['SNR'] = [l[4+ii]]
      for key in data:
         for key2 in data[key]:
            data[key][key2] = numpy.array(data[key][key2])
      
      return(data,'CSP')

   def get_SN_spectra(self):
      '''Get the spectra form the SQL database.  Returns a tuple:
      (JDs, waves, fluxes)
      JDs is 1D array of Julian Dates
      waves is a list of 1D arrays of wavelengths
      fluxes is a list of 1D arrays of fluxes.'''
      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")
      slct1 = '''SELECT %s,%s from %s where %s=%%s''' % \
            (self.SPEC_JD,self.SPEC_INDEX,self.SPEC_INFO,self.SPEC_NAME)
      slct2 = '''SELECT %s,%s from %s where %s=%%s and %s=%%s order by %s''' % \
            (self.SPEC_LAMB,self.SPEC_FLUX, self.SPEC_TABLE,self.SPEC_NAME,
             self.SPEC_INDEX, self.SPEC_LAMB)
      N = self.c2.execute(slct1, (self.name))
      if N == 0:
         raise ValueError("No spectroscopy for %s found" % \
                  (self.name))
      list = self.c2.fetchall()

      JDs = []
      waves = []
      fluxes = []
      for jd,file in list:
         JDs.append(jd)
         self.c2.execute(slct2, (self.name, file))
         data = self.c2.fetchall()
         waves.append(numpy.array([datum[0] for datum in data]))
         fluxes.append(numpy.array([datum[1] for datum in data]))
      return(numpy.array(JDs),waves,fluxes)

   def update_photometry(self, filter, times, attr, values, tol=1e-6):
      '''Update an entries in the SN photometry table.  The affected rows are
      those for which filter matches and the times parameter is less than
      tol from the values in the table.  attr for row i is changed to
      values[i]'''
      if type(self.PHOTO_TABLE) is type([]):
         raise TypeError("This function is not available on table joins")

      if not self.connected:
         raise RuntimeError("Not connected to SQL database.")

      if not numpy.shape(times) == numpy.shape(values):
         raise ValueError("shapes of times and values must agree")
      if len(numpy.shape(times)) == 0:
         times = numpy.array([times])
         values = numpy.array([values])

      if attr in self.PHOTO_ATTR_KEYS:
         field = self.PHOTO_ATTR_KEYS[attr]
      else:
         field = attr
      
      if field not in self.PHOTO_table_info:
         raise AttributeError("Photometry table has no field %s")

      # First, find the primary key for the PHOTO_TABLE:
      PRI = None
      for f in self.PHOTO_table_info:
         if self.PHOTO_table_info[f]['key'] == "PRI":
            PRI = f
            break
      if PRI is not None:
         slct = '''SELECT %s,%s from %s where %s=%%s and %s=%%s %s''' % \
               (self.PHOTO_JD, PRI, self.PHOTO_TABLE, self.PHOTO_ID,
               self.PHOTO_FILT, self.PHOTO_COND)
         self.c.execute(slct, (self.name, filter))
         l = self.c.fetchall()
         ids = numpy.array([item[1] for item in l])
         db_time = numpy.array([item[0] + self.JD_OFFSET for item in l])
         con = numpy.less(numpy.absolute(db_time[numpy.newaxis,:] - times[:,numpy.newaxis]),
               tol)
         num_match = numpy.sum(con, axis=1)
         # Check to see there's at least one match per requested time:
         if not numpy.all(num_match):
            bids = numpy.nonzero(numpy.equal(num_match, 0))[0]
            btimes = ""
            for i in range(len(bids)):  btimes += " %.1f " % (times[bids[i]])
            raise AttributeError("No photometry entry for %s" % btimes)

         # check if multiple matches:
         if numpy.any(numpy.greater(num_match,1)):
            bids = numpy.nonzero(numpy.greater(num_match,1))[0]
            btimes = ""
            for i in range(len(bids)):  btimes += " %.1f " % (times[bids[i]])
            raise AttributeError("Multiple time matches for %s" % btimes)

         # Should be all clear now
         match_ids = (ids[numpy.newaxis,:]*numpy.ones(numpy.shape(times))[:,numpy.newaxis])[con]
         
         # Do the update
         updt = '''UPDATE %s set %s = %%s where %s = %%s''' % (self.PHOTO_TABLE,
                field, PRI)
         for i in range(len(times)):
            self.c.execute(updt, (values[i], match_ids[i]))


      else:
         # Gotta do it the easier, but slower way
         updt = '''UPDATE %s set %s = %%s where %s = %%s and %s = %%s and abs(%s - %%s) < %f''' %\
               (self.PHOTO_TABLE, field, self.PHOTO_ID, self.PHOTO_FILT,
                self.PHOTO_JD, tol)
         for i in range(len(times)):
            self.c.execute(updt, (values[i], self.name, filter, times[i]-self.JD_OFFSET))


class sql_highz(sqlbase):
   host = "kepler.obs.carnegiescience.edu"
   user = "CSP"
   passwd = None
   port = 3306

   PHOTO_DB = "SN"
   SPEC_TABLE = "SPECTRA"
   SPEC_INFO = "SP_INFO"

   #FILTER_KEYS = {'Bs':'B',
   #               'Vs':'V'}


class sql_lowz(sqlbase):
   host = "csp2.lco.cl"
   user = "CSP"
   passwd = None
   PHOTO_DB = "Phot"
   port = 3306
   readonly = 1

   #FILTER_KEYS = {'B':'Bs',    # Over-ride Swope B and V
   #               'V':'Vs'}

   SN_TABLE = "SNList"
   SN_ID = "SN"
   ATTR_KEYS = {'z':'zc/300000.0',
                'ra':'ra*15.0',
                'decl':'de'}
   PHOTO_TABLE = "SNPHOT"
   PHOTO_ID = "SN"
   PHOTO_JD = "JD"
   PHOTO_MAG = "MNAT"
   PHOTO_EMAG = "NATERR"
   PHOTO_FILT = "FILT"
   PHOTO_K = None      # No K-corrections in the DB
   PHOTO_COND = "and MNAT is not NULL and MNAT > 0"
   JD_OFFSET = 52999.5    # database is JD - 2453000, so +2999.5 gives MJD

class sql_oldlocal(sqlbase):
   host = "kepler.obs.carnegiescience.edu"
   user = "CSP"
   passwd = None
   PHOTO_DB = "CSP"
   port = 3306
   readonly = 1

   #FILTER_KEYS = {'B':'Bs',    # Over-ride Swope B and V
   #               'V':'Vs'}

   SN_TABLE = "SNList"
   SN_ID = "SN"
   ATTR_KEYS = {'z':'zc/300000.0',
                'ra':'ra*15.0',
                'decl':'de'}
   PHOTO_TABLE = "SNPHOT"
   PHOTO_ID = "SN"
   PHOTO_JD = "JD"
   PHOTO_MAG = "MNAT"
   PHOTO_EMAG = "NATERR"
   PHOTO_FILT = "FILT"
   PHOTO_K = None      # No K-corrections in the DB
   PHOTO_COND = "and MNAT is not NULL and MNAT > 0"
   JD_OFFSET = 52999.5    # database is JD - 2453000, so +2999.5 gives MJD

class sql_csp2(sqlbase):
   host = "csp2.lco.cl"
   user = "cburns"
   passwd = None
   PHOTO_DB = "Phot"
   port = 3306
   readonly = 1

   SN_TABLE = "SNList"
   SN_ID = "SN"
   SN_ID2 = "NAME_CSP"
   ATTR_KEYS = {'z':'zc/300000.0',
                'ra':'ra*15.0',
                'decl':'de'}
   PHOTO_TABLE = ["MAGSN","MAGINS"]
   PHOTO_ID = "MAGSN.field"
   PHOTO_JD = "MAGSN.JD"
   PHOTO_MAG = "MAGSN.mag"
   PHOTO_EMAG = "sqrt(MAGSN.err*MAGSN.err + MAGSN.fiterr*MAGSN.fiterr)"
   PHOTO_FILT = "CONCAT(tel,ins,MAGSN.filt)"
   #PHOTO_FILT = "CASE WHEN (ins='RC' and MAGSN.filt='J' and MAGSN.JD < 2454846.0) THEN 'J' WHEN (ins='WI' and MAGSN.filt='Y') THEN 'Ydw' WHEN (ins='RC' and MAGSN.filt='J' and MAGSN.JD >= 2454846.0) THEN 'Jrc2' ELSE MAGSN.filt END"
   PHOTO_K = None      # No K-corrections in the DB
   PHOTO_SNR = "MAGINS.flux/sqrt(MAGINS.flux + 200*MAGINS.sky)"
   PHOTO_COND = "and MAGSN.obj=-1 and MAGSN.mag > 0 and MAGSN.fits=MAGINS.fits and MAGSN.obj=MAGINS.obj"
   JD_OFFSET = -2400000.5    # database is JD, JD-2400000.5 gives MJD

   FILTER_KEYS = {
         'BAAFSJ':'JFS',
         'BAAFSH':'HFS',
         'BAAFSY':'J1FS',
         'DUPRCY':'Yd',
         'DUPRCJ':'Jd',
         'DUPRCH':'Hd',
         'DUPWIY':'Ydw',
         'DUPWIJ':'J',
         'DUPWIH':'H',
         'DUPWIK':'Kdw',
         'DUPDCB':'B',
         'DUPDCV':'V',
         'DUPDCu':'u',
         'DUPDCg':'g',
         'DUPDCr':'r',
         'DUPDCi':'i',
         'SWODCB':'B',
         'SWODCu':'u',
         'SWODCg':'g',
         'SWODCr':'r',
         'SWODCi':'i',
         'SWONCB':'B2',
         'SWONCV':'V2',
         'SWONCu':'u2',
         'SWONCg':'g2',
         'SWONCr':'r2',
         'SWONCi':'i2',
         'SWORCY':'Y',
         'SWORCH':'H'}
   
   def filter_key(self, name, epoch):
      if name in self.FILTER_KEYS:
         return self.FILTER_KEYS[name]

      if name == 'SWORCJ':
         # J was changed
         if epoch < 54845.0:
            return 'J'
         else:
            return 'Jrc2'
      if name == 'SWODCV':
         if epoch < 53748.0:
            return 'V0'
         elif 53748.0 < epoch < 53759.0:
            return 'V1'
         else:
            return 'V'
      print("Warning: don't know filter {}".format(name))
      return name[-1]

class sql_SBS_csp2(sql_csp2):
   host = 'obsns09.obs.carnegiescience.edu'
   user = 'CSP'

   PHOTO_DB = "CSP"
   SPEC_DB = "CSP"
   SPEC_TABLE = "SPECTRA"
   SPEC_INFO = "SP_INFO"
   SPEC_JD = "JD"
   SPEC_LAMB = "LAMBDA"
   SPEC_FLUX = "FLUX"
   SPEC_INDEX = "FILE"
   SPEC_NAME = "SN"

class sql_csp2_pub(sql_csp2):
   PHOTO_DB = "PubPhot"

class sql_SBS_csp2_pub(sql_csp2):
   host = 'kepler.obs.carnegiescience.edu'
   user = 'CSP'
   PHOTO_DB = "PubPhot"


databases = \
   {'default':(sql_csp2, "Working CSP2 database at LCO"),
    'SBS':(sql_SBS_csp2, "Working CSP2 database at SBS"),
    'LCOpub':(sql_csp2_pub, "Published CP2 database at LCO"),
    'SBSpub':(sql_SBS_csp2_pub, "Published CSP2 database at SBS"),
    'highz':(sql_highz, "Highz database at SBS")}

default_sql = databases['default'][0]()

def setSQL(name):
   global default_sql
   if name not in databases:
      raise ValueError("Error, unknown database %s" % name)
   default_sql = databases[name][0]()
