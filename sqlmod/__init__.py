import os,sys
import numpy
import getpass
#from snpy import filters

try:
   import MySQLdb as sql
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
   db = "database"         # database containing SN info
   readonly = 0            # Can we read only?

   SN_TABLE = "SNe"        # SQL Table of supernova attributes
   SN_ID = "name"          # Unique SN ID field in SN_TABLE
   SN_COND = ""            # Any extra SN query conditions
   PHOTO_TABLE = "Photo"   # SQL Table of photometry
   PHOTO_ID = "name"       # Unique SN ID field in PHOTO_TABLE
   PHOTO_FILT = "filter"   # SQL field for filter name in PHOTO_TABLE

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
   PHOTO_K = "K"           # SQL field in PHOTO_TABLE that gives K-corrections
   PHOTO_FILT = "filter"   # SQL field in PHOTO_TABLE that gives filter name
   PHOTO_COND = ""         # Any extra WHERE conditions
   
   def __init__(self):
      self.connected=0

      # Allow for environment variable override
      for att in ['user','passwd','db','host','port']:
         key = 'SQL_' + att.upper()
         if key in os.environ:
            self.__dict__[att] = os.environ[key]

   def __getstate__(self):
      # Need to do this because the sql connection instances
      # themselves are not pickleable.
      dict = self.__dict__.copy()
      if 'con' in dict:  del dict['con']
      if 'c' in dict:  del dict['c']
      dict['connected'] = 0
      return dict


   def connect(self, name):
      '''Connect to the database.'''
      if self.passwd is None:
         passwd = getpass.getpass(prompt='SQL passwd for %s@%s:\n' % (self.user,
            self.host))
      else:
         passwd = self.passwd
      if not self.connected:
         try:
            self.con = sql.connect(host=self.host, user=self.user, passwd=passwd,
                  db=self.db, port=self.port)
         except:
            raise
         else:
            # Cache the passwd for future use if it worked
            if self.passwd is None:  self.passwd = passwd
            
         self.c = self.con.cursor()
         self.connected = 1
         self.name = name
         # Collect info about the SN_TABLE and PHOTO_TABLE:
         self.SN_table_info = self.get_table_info(self.SN_TABLE)
         self.PHOTO_table_info = self.get_table_info(self.PHOTO_TABLE)

         # See if the SN object exists
         slct = '''SELECT * from %s where %s = %%s %s''' % (self.SN_TABLE, self.SN_ID,
               self.SN_COND)
         N = self.c.execute(slct, (self.name))
         if N == 0:
            return 0
         if N > 1:
            print "Warning!  %s is not unique in the database, taking first" %\
                  (self.name)
            return N
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
         raise RuntimeError, "Not connected to SQL database."
      self.c.execute('''DESCRIBE %s''' % table)
      l = self.c.fetchall()
      data = {}
      for item in l:
         data[item[0]] = {'type':item[1],
                          'null':item[2],
                          'key':item[3],
                          'default':item[4],
                          'extra':item[5]}
      return(data)

   def close(self):
      '''Close connection to the database.'''
      self.con.close()
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
         raise RuntimeError, "Not connected to SQL database."
      insrt = '''INSERT INTO %s (%s,%s,%s,%s) VALUES (%%s,%%s,%%s,%%s)''' % \
            (self.SN_TABLE, self.SN_ID, self.sql_field('name'),
             self.sql_field('ra'), self.sql_field('decl'), self.sql_field('z'))
      self.c.execute(insrt, (self.name, ra, decl, z))

   def create_SN_photometry(self, data):
      '''Insert the photomtery data into the database.  See get_SN_photometry
      for an explanation of data'''
      if not self.connected:
         raise RuntimeError, "Not connected to SQL database."
      filters = data.keys()
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

   def get_SN_parameter(self, attr):
      '''Gets a single attribute of a SN in the database.'''
      if not self.connected:
         raise RuntimeError, "Not connected to SQL database."

      if attr not in self.SN_table_info and attr not in self.ATTR_KEYS:
         raise AttributeError, "SN table does not have this attribute: %s" % \
               attr
      field = self.sql_field(attr)

      slct = '''SELECT %s from %s where %s = %%s %s''' % (field, self.SN_TABLE, 
            self.SN_ID, self.SN_COND)
      N = self.c.execute(slct, (self.name))
      if N == 0:
         raise AttributeError, "Attribute %s(%s) not found in SQL db" % \
               (attr,field)
      return(self.c.fetchall()[0][0])

   def set_SN_parameter(self, attr, value):
      '''Sets a single attribute of a SN in the datase.'''
      if self.readonly:
         raise ValueError, "Database is read-only"
      if not self.connected:
         raise RuntimeError, "Not connected to SQL database."
      if attr not in self.SN_table_info.keys() and attr not in self.ATTR_KEYS:
         raise AttributeError, "SN table does not have this attribute: %s" % \
               attr
      field = self.sql_field(attr)
      if field not in self.SN_table_info:
         raise AttributeError, "SN table does not have this attribute: %s" % \
               field

      updt = '''UPDATE %s SET %s = %%s where %s = %%s''' % \
            (self.SN_TABLE, field, self.SN_ID)
      self.c.execute(updt, (value, self.name))

   def get_SN_photometry(self):
      '''Get the photometry form the SQL database.  Returns a dictionary
      indexed by filter name.  Each element is a 4-tuple of arrays:
      (MJD, mag, e_mag, K).'''
      if not self.connected:
         raise RuntimeError, "Not connected to SQL database."
      if self.PHOTO_K is not None:
         slct = '''SELECT %s,%s,%s,%s,%s from %s where %s=%%s %s ORDER by %s''' % \
               (self.PHOTO_FILT, self.PHOTO_JD, self.PHOTO_MAG, self.PHOTO_EMAG,
                self.PHOTO_K, self.PHOTO_TABLE, self.PHOTO_ID, self.PHOTO_COND,
                self.PHOTO_JD)
      else:
         slct = '''SELECT %s,%s,%s,%s from %s where %s=%%s %s ORDER by %s''' % \
               (self.PHOTO_FILT, self.PHOTO_JD, self.PHOTO_MAG, self.PHOTO_EMAG,
                self.PHOTO_TABLE, self.PHOTO_ID, self.PHOTO_COND, self.PHOTO_JD)
      N = self.c.execute(slct, (self.name))
      if N == 0:
         raise ValueError, "No photometry for %s found" % \
                  (self.name)
      list = self.c.fetchall()

      data = {}
      for l in list:
         filter = l[0]
         if l[0] not in self.FILTER_KEYS:  
            filter = l[0]
         else:
            filter = self.FILTER_KEYS[l[0]]
         if filter in data:
            data[filter][0].append(l[1] + self.JD_OFFSET)
            data[filter][1].append(l[2])
            data[filter][2].append(l[3])
            if self.PHOTO_K is not None:
               data[filter][3].append(l[4])
         else:
            if self.PHOTO_K is not None:
               data[filter] = [[l[1] + self.JD_OFFSET],
                               [l[2]], [l[3]], [l[4]]]
            else:
               data[filter] = [[l[1] + self.JD_OFFSET],
                               [l[2]], [l[3]], None]
      for key in data:
         for i in range(3):
            data[key][i] = numpy.array(data[key][i])
         if data[key][-1] is not None:
            data[key][-1] = numpy.array(data[key][-1])
      return(data)

   def update_photometry(self, filter, times, attr, values, tol=1e-6):
      '''Update an entries in the SN photometry table.  The affected rows are
      those for which filter matches and the times parameter is less than
      tol from the values in the table.  attr for row i is changed to
      values[i]'''
      if not self.connected:
         raise RuntimeError, "Not connected to SQL database."

      if not numpy.shape(times) == numpy.shape(values):
         raise ValueError, "shapes of times and values must agree"
      if len(numpy.shape(times)) == 0:
         times = numpy.array([times])
         values = numpy.array([values])

      if attr in self.PHOTO_ATTR_KEYS:
         field = self.PHOTO_ATTR_KEYS[attr]
      else:
         field = attr
      
      if field not in self.PHOTO_table_info:
         raise AttributeError, "Photometry table has no field %s"

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
         if not numpy.alltrue(num_match):
            print type(num_match)
            bids = numpy.nonzero(numpy.equal(num_match, 0))[0]
            print type(bids)
            btimes = ""
            for i in range(len(bids)):  btimes += " %.1f " % (times[bids[i]])
            raise AttributeError, "No photometry entry for %s" % btimes

         # check if multiple matches:
         if numpy.sometrue(numpy.greater(num_match,1)):
            bids = numpy.nonzero(numpy.greater(num_match,1))[0]
            btimes = ""
            for i in range(len(bids)):  btimes += " %.1f " % (times[bids[i]])
            raise AttributeError, "Multiple time matches for %s" % btimes

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
   host = "newton.obs.carnegiescience.edu"
   user = "CSP"
   passwd = None
   db = "SN"
   port = 3306

   FILTER_KEYS = {'Bs':'B',
                  'Vs':'V'}


class sql_lowz(sqlbase):
   host = "csp2.lco.cl"
   user = "cburns"
   passwd = None
   db = "Phot"
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
   JD_OFFSET = 2999.5    # database is JD - 2453000, so +2999.5 gives MJD

default_sql = None
