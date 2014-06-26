# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2014
#  Fredrik Strömberg <fredrik314@gmail.com>,
#  Stephan Ehlen <stephan.j.ehlen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPLv2)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

r"""
Programs for inserting modular forms data into a mongo database. By default we are also synchronizing with a file-system based datbase.

"""


import sage.all
import os,re
import pymongo
import gridfs
from compmf.filesdb import FilenamesMFDBLoading
from compmf.compute import ComputeMFData
from compmf.character_conversions import conrey_from_sage_character_number,sage_character_galois_orbit_rep_from_number
from sage.all import prime_pi,parallel,loads,dimension_new_cusp_forms,RR,ceil,load,dumps,save,euler_phi


from compmf import clogger

class MongoMF(object):
    def __init__(self,host='localhost',port=37010,verbose=0,db='modularforms2',**kwds):
        r"""
        Mongo database for modular forms.
        
        INPUT:

        - host    -- string:  host where mongodb is located
        - port    -- integer: port where we connect ot the mongodb
        - db      -- string:  name of the mongo database to use.
        - verbose -- integer: set verbosity

        KEYWORDS:
        
        - compute -- Boolean (default True) set to False if you do not want to compute anything.
        """
        self._mongodb = pymongo.Connection('{0}:{1}'.format(host,port))[db]
        self._mongo_conn = pymongo.Connection('{0}:{1}'.format(host,port))

        ## Our databases
      
        self._modular_symbols_collection = 'Modular_symbols'
        self._newform_factors_collection = 'Newform_factors'
        self._aps_collection = 'ap'
        self._atkin_lehner_collection = 'Atkin_Lehner'
        self._modular_symbols = self._mongodb["{0}.files".format(self._modular_symbols_collection)]
        self._newform_factors = self._mongodb["{0}.files".format(self._newform_factors_collection)]
        self._aps = self._mongodb["{0}.files".format(self._aps_collection)]
        self._atkin_lehner = self._mongodb["{0}.files".format(self._atkin_lehner_collection)]

        self._sage_version = sage.version.version

    def __repr__(self):
        r"""
        String representation of self.
        """
        s="Modular forms database at mongodb: {1}".format(self._db._data,self._mongo_conn)
        return s


    def show_existing_mongo(self,db='fr'):
        r"""
        Display an overview of the existing records in the mongo db.
        """
        files = self._modular_symbols
        factors = self._newform_factors 
        levels = files.distinct('N')
        weights = files.distinct('k')
        print "Modular symbols: \n"
        if levels<>[]:
            print "{0} records with levels in range {1} -- {2}".format(files.count(),min(levels),max(levels))
        print "Newforms: \n"
        levels = factors.distinct('N')
        weights = factors.distinct('k')
        if levels<>[]:
            print "{0} records with levels in range {1} -- {2}".format(factors.count(),min(levels),max(levels))        

          ### Routines for accessing the objects stored in the mongo database.

    def load_from_mongo(self,col,fid):
        r"""
        Load a file from gridfs collection col and id nr. fid.
        
        """
        col_files = '{0}.files'.format(col)
        if not col_files in self._mongodb.collection_names():
            return None
        if self._mongodb[col_files].find({'_id':fid}).count()==0:
            return None
        fs = gridfs.GridFS(self._mongodb,col)
        return loads(fs.get(fid).read())
       
    def get_ambient(self,N,k,i,**kwds):
        r"""
        Return the ambient space M(N,k,i)

        """
        ambient_id = kwds.get('ambient_id',None)
        if ambient_id is None:
            if kwds.get('compute',True):
                ambient_id = self.compute_ambient(N,k,i,**kwds)
            else:
                ids = self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(i)}).distinct('_id')
                if ids==[]:
                    return None
                ambient_id = ids[0]
        return self.load_from_mongo('Modular_symbols',ambient_id)
       
    def number_of_factors(self,N,k,i,**kwds):
        r"""
        Return the number of factors. 
        """
        return self._newform_factors.find({'N':int(N),'k':int(k),'chi':int(i)}).count()

    def get_factors(self,N,k,i,d=None):
        r"""
        Get factor nr. d of the space M(N,k,i)
        """
        s = {'N':int(N),'k':int(k),'chi':int(i)}
        if not d is None:
            s['newform']=int(d)
        res = {}
        for r in self._newform_factors.find(s):
            d = r['newform']
            fid = r['_id']
            f = self.load_from_mongo('Newform_factors',fid)
            t = (int(N),int(k),int(i),int(d))
            res[t] = f
        return res
    
    def get_aps(self,N,k,i,d=None,character_naming='sage'):
        r"""
        Get the lists of Fourier coefficients for the space M(n,k,i) and orbit nr. d
        
        """
        if character_naming=='sage':
            s = {'N':int(N),'k':int(k),'chi':int(i)}
        else:
            ## Fetch the record in the database corresponding to the Galois orbit of
            ## chi_N,i  (in Conrey's character naming scheme)
            s = {'N':int(N),'k':int(k),'character_galois_orbit':{"$all":[int(i)]}}
        if not d is None:
            s['newform']=int(d)
        res = {}
        for r in self._aps.find(s):
            fid=r['_id']
            d = r['newform']
            meta = {'cputime':r.get('cputime'),'version':r.get('sage_version')}
            E,v = self.load_from_mongo('ap',fid)
            t = (int(N),int(k),int(i),int(d))
            if not res.has_key(t):
                res[t]=[]
            res[t].append((E,v,meta))
        return res
        
class CompMF(MongoMF):
    r"""
    Class for computing modular forms and inserting records in a mongodb as well as a file system database.
    """
    def __init__(self,datadir='',host='localhost',port=37010,verbose=0,db='modularforms2',**kwds):
        r"""

        INPUT:

        - datadir -- string:  root directory of the file system database
        - host    -- string:  host where mongodb is located
        - port    -- integer: port where we connect ot the mongodb
        - db      -- string:  name of the mongo database to use.
        - verbose -- integer: set verbosity

        KEYWORDS:
        
        - compute -- Boolean (default True) set to False if you do not want to compute anything.
        - save_to_file -- Boolean (default True) set to False if you do not want to save to file


        """
        super(CompMF,self).__init__(host,port,verbose,db,**kwds)
        self._db = FilenamesMFDBLoading(datadir)
        self._computedb = ComputeMFData(datadir)
        assert str(db).isalnum()
        
        self._do_computations = kwds.get('compute',True)
        self._save_to_file = kwds.get('save_to_file',True)
        
    def __repr__(self):
        r"""
        String representation of self.
        """
        s="Connecting files at {0} and mongodb: {1}".format(self._db._data,self._mongo_conn)
        return s

    
    

    def show_existing_files(self):
        r"""
        Display an overview of the existing records in the files database.
        """
        s=""
        print "Directory:{0}".format(self._db._data)
        s+=self._db.show_existing_files_db()
        return s
    
    ## The main routines for inserting records into the mongo database.
    ## We insert records from the files database and if these records are incomplete
    ## (or not computed at all) we first compute them.
    
    def convert_to_mongo_all(self,par=0,**kwds):
        r"""
        Convert all records in the files to MongoDB
        """
        nrange = self._db.known_levels()
        nrange.sort()
        clogger.debug("nrange={0}".format(nrange))
        return self.convert_to_mongo(nrange)
           
    def convert_to_mongo(self,N=None,k=None,ncpus=1,**kwds):
        r"""
        Converting records which exists (possibly partly) in the files for a given level and insert into the mongodb.
        """
        k0 = kwds.get('k',None)
        clogger.debug("Converting N={0} and k={1}".format(N,k))
        s = ""
        if N<>None:
            s = "N={0}".format(N)
        if k<>None:
            s+= " k={0}".format(k)            
        
        args = []
        if kwds.get('trivial'):
            for (N,k,i,newforms,nap) in self._db.known(s):
                if i==0:
                    args.append((N,k,i))
        else:
            for (N,k,i,newforms,nap) in self._db.known(s):
                args.append((N,k,i))
        self.get_or_compute_spaces(args,**kwds)
        
    def get_or_compute_spaces(self,args,**kwds):
        r"""
        Get or compute a list of spaces in parallel.
        
        """
        ncpus = kwds.get('ncpus',0)
        if ncpus==8:
            return list(self._compute_and_insert_one_space8(args,**kwds))
        elif ncpus==16:
            return list(self._compute_and_insert_one_space16(args,**kwds))            
        elif ncpus==32:
            return list(self._compute_and_insert_one_space32(args,**kwds))                        
        else:
            return [self.compute_and_insert_one_space(x[0],x[1],x[2],**kwds) for x in args]

    ## Different levels of parallelization
    @parallel(ncpus=8)
    def _compute_and_insert_one_space8(self,N,k,i,**kwds):
        return self.compute_and_insert_one_space(N,k,i,**kwds)        
    @parallel(ncpus=16)
    def _compute_and_insert_one_space16(self,N,k,i,**kwds):
        return self.compute_and_insert_one_space(N,k,i,**kwds)
    @parallel(ncpus=32)
    def _compute_and_insert_one_space32(self,N,k,i,**kwds):
        return self.compute_and_insert_one_space(N,k,i,**kwds)        
        
    def compute_and_insert_one_space(self,N,k,i,**kwds):
        r"""
        If data for the given space M(N,k,i) exists in either mongo or files database we fetch this data (and possible complete if e.g. more coefficients are needed) and then insert the result into both databases unless explicitly told not to.

        """
        clogger.debug("converting {0}".format(N,k,i))
        clogger.debug("Computing ambient modular symbols")
        if kwds.get('Nmax',0)<>0 and kwds.get('Nmax')>N:
            return 
        ambient_fid = self.compute_ambient(N,k,i,**kwds)
        # Insert ambient modular symbols
        kwds['ambient_id']=ambient_fid
        clogger.debug("Getting factors!")
        factor_ids = self.compute_factors(N,k,i,**kwds)
        kwds['factor_ids']=factor_ids 
        # Necessary for L-function computations.
        
        pprec = precision_needed_for_L(N,k,pprec=100)
        clogger.debug("Computing aps to prec {0}".format(pprec))
        aps = self.compute_aps(N,k,i,pprec,**kwds)
        clogger.debug("Computing atkin lehner to prec {0} for i={0}".format(pprec,i))
        try:
            atkin_lehner = self.compute_atkin_lehner(N,k,i,**kwds)
        except:
            pass 
        return True    

    def compute_atkin_lehner(self,N,k,i,**kwds):
        r"""
        Compute the Atkin-Lehner eigenvalues of space (N,k,i).
        """
        verbose = kwds.get('verbose',0)
        c = conrey_from_sage_character_number(N,i)
        ci = c.number()
        if not c.is_trivial(): #or c.multiplicative_order()==2):
            return []
        al_in_mongo = self._atkin_lehner.find({'N':int(N),'k':int(k),'chi':int(i),'cchi':int(ci)}).distinct('_id')
        fs = gridfs.GridFS(self._mongodb, 'Atkin_Lehner')
        if len(al_in_mongo)==0:
            ambient = self.get_ambient(N,k,i,**kwds)
            number_of_factors = self.number_of_factors(N,k,i)
            self._computedb.compute_atkin_lehner(N,k,i,M=ambient,m=number_of_factors,verbose=verbose)
            m = self._computedb._db.number_of_known_factors(N,k,i)
            for d in range(m):

                atkin_lehner_file = self._computedb._db.factor_atkin_lehner(N,k,i,d,True)
                #atkin_lehner = load(atkin_lehner_file)
                try:
                    atkin_lehner = open(atkin_lehner_file,'r').read()
                except IOError:
                    raise ArithmeticError,"Error with opening the Atkin-Lehner file {0}! chi order={1}".format(atkin_lehner_file,c.multiplicative_order())
                try:
                    meta = load(self._db.meta(atkin_lehner_file))
                except IOError:
                    meta = {}
                fname = 'atkin_lehner_evs-{0:>5}-{1:>3}-{2:>3}-{3:>3}'.format(N,k,i,d)
                fid = fs.put(dumps(atkin_lehner),filename=fname,
                             N=int(N),k=int(k),chi=int(i),newform=int(d),cchi=int(ci),
                             cputime = meta.get("cputime",""),
                             sage_version = meta.get("version",""))
                al_in_mongo.append(fid)
        return al_in_mongo
        
        
    def compute_ambient(self,N,k,i,**kwds):
        r"""
        Compute the ambient space with the given parameters and insert it in mongo if it is not there 
        """
        files_ms = self._modular_symbols
        fs_ms = gridfs.GridFS(self._mongodb, 'Modular_symbols')
        verbose = kwds.get('verbose',0)
        ci = conrey_from_sage_character_number(N,i).number()
        save_in_file = kwds.get('save_in_file',True)
        compute = kwds.get('compute',self._do_computations)
        # We first see if this space is already in the mongo database.
        rec = files_ms.find_one({'N':int(N),'k':int(k),'chi':int(i)})
        ambient_in_mongo = 0
        ambient_in_file = 0
        ambient = None
        cputime = None
        clogger.debug("In compute ambient: {0}".format((N,k,i)))
        if rec<>None:
            ambient_in_mongo = rec['_id']
            # Check factors 
            clogger.debug("Have ambient space!")
            clogger.debug("ambient_in_mongo={0}".format(rec))
        # Next see if we have it in the files database
        try:
            ambient = self._db.load_ambient_space(N,k,i)
            ambient_in_file = 1
        except ValueError:
            pass
        if ambient_in_mongo <> 0 and ambient_in_file==0:
            ambient = loads(fs_ms.get(ambient_in_mongo).read())
            cputime = rec.get('cputime')

            clogger.debug("Space {0},{1},{2} is in mongo but not in file!".format(N,k,i))
            clogger.debug("ambient={0}".format(ambient))
        if ambient_in_file == 0: 
            if not compute:
                clogger.debug("Space {0},{1},{2} not computed at all!".format(N,k,i))
                return None
            clogger.debug("Compute and/or save to file!")
            self._computedb.compute_ambient_space(N,k,i,M=ambient,tm=cputime)
            ambient_in_file = 1
#        num_factors = len(self.compute_factors(N,k,i,**kwds))
        if ambient_in_mongo == 0 and not ambient is None:
            metaname = self._db.space(N,k,i,False)+"/ambient-meta.sobj"
            clogger.debug("metaname={0}".format(metaname))
            try:
                meta = load(metaname)
            except (OSError,IOError):
                meta={}
            fname = "gamma0-ambient-modsym-{0}".format(self._db.space_name(N,k,i))
            ## Note that we have to update the number of orbits.
            fid = fs_ms.put(dumps(ambient),filename=fname,
                            N=int(N),k=int(k),chi=int(i),orbits=0,
                            cchi=int(ci),
                            cputime = meta.get("cputime",""),
                            sage_version = meta.get("version",""))
        else:
            fid = ambient_in_mongo
        if fid == 0:
            fid = None
        return fid

    def compute_factors(self,N,k,i,**kwds):
        r"""
        Compute / store newform factors of parameters N,k,i
        """
        verbose = kwds.get('verbose',0)
        ambient_id = kwds.get('ambient_id',None)
        compute = kwds.get('compute',self._do_computations)
        if ambient_id is None:
            ambient_id = self.compute_ambient(N,k,i,**kwds)
        ci = conrey_from_sage_character_number(N,i).number()
        if ambient_id is None:
            clogger.debug("No ambient space!")
            ambient_id = self.compute_ambient(N,k,i,**kwds)
        files_fact = self._newform_factors
        files_ms = self._modular_symbols
        fs_fact = gridfs.GridFS(self._mongodb, 'Newform_factors')
        factors_in_file = self._db.number_of_known_factors(N,k,i)
        factors_in_mongo = files_fact.find({'ambient_id':ambient_id}).distinct('_id')
        clogger.debug("factors_in_file={0}".format(factors_in_file))
        clogger.debug("factors_in_mongo={0}".format(factors_in_mongo))
        if factors_in_file == 0 and  len(factors_in_mongo)==0:
            if not compute:
                clogger.debug("No factors exist and we do not compute anything!")
                return None
                # Else compute and insert into the files.
            factors_in_file = self._computedb.compute_decompositions(N,k,i)
            clogger.debug("Computing factors! m={0}".format(factors_in_file))

        if len(factors_in_mongo)==0:
            fname = "gamma0-factors-{0}".format(self._db.space_name(N,k,i))
            clogger.debug("Inserting factors into mongo! fname={0}".format(fname))
            num_factors_in_file = self._db.number_of_known_factors(N,k,i)
            ambient = self.get_ambient(N,k,i,ambient_id=ambient_id)
            for d in range(num_factors_in_file):
                try: 
                    factor = self._db.load_factor(N,k,i,d,M=ambient)
                except RuntimeError:
                    ## We probably need to recompute the factors
                    clogger.debug("The factors from file was corrupt / empty. Need to compute them anew!")
                        
                    factors_in_file = self._computedb.compute_decompositions(N,k,i)
                    try:
                        factor = self._db.load_factor(N,k,i,d,M=ambient)
                    except RuntimeError:
                        raise ArithmeticError,"Could not get factors for {0}".format((N,k,i))
                metaname = self._db.space(N,k,i,False)+"/decomp-meta.sobj"
                clogger.debug("metaname={0}".format(metaname))
                try:
                    meta = load(metaname)
                except (OSError,IOError):
                    meta={}
                if factor==None:
                    clogger.debug("Factor {0},{1},{2},{3} not computed!".format(N,k,i,d))
                    continue
                fname1 = "{0}-{1:0>3}".format(fname,d)
                facid = fs_fact.put(dumps(factor),filename=fname1,
                                    N=int(N),k=int(k),chi=int(i),
                                    cchi=int(ci),
                                    newform=int(d),
                                    cputime = meta.get("cputime",""),
                                    sage_version = meta.get("version",""),
                                    ambient_id=ambient_id,
                                    v=int(1))
                if not facid is None:
                    factors_in_mongo.append(facid)
                clogger.debug("inserted factor: {0},{1}".format(d,facid))
        ambient_files = self._modular_symbols
        ambient_files.update({'_id':ambient_id},{"$set":{'orbits':len(factors_in_mongo)}})
        if factors_in_file == 0:
            D = []
            for fid in factors_in_mongo:
                D.append(loads(fs_fact.get(fid).read()))
                self._computedb.compute_decompositions(N,k,i,D=D)
                    
        return factors_in_mongo


 
        
    def compute_aps(self,N,k,i,pprec=None,**kwds):
        r"""
        Compute & store aps
        """       
        if pprec is None:
            pprec = precision_needed_for_L(N,k)
        ambient_id = kwds.get('ambient_id',None)
        if ambient_id is None:
            ambient_id = self.compute_ambient(N,k,i,**kwds)
            kwds['ambient_id']=ambient_id
        ambient = self.get_ambient(N,k,i,ambient_id=ambient_id,verbose=0)
        num_factors = len(kwds.get('factors_ids',self.compute_factors(N,k,i,**kwds)))

        compute = kwds.get('compute',self._do_computations)
        verbose = kwds.get('verbose')
        c = conrey_from_sage_character_number(N,i)
        ci = c.number()        
        if N > 1: # There is a bug in Conrey character mod 1
            orbit = [int(x.number()) for x in c.galois_orbit()]
            orbit.sort()
        else:
            orbit = [int(1)]
        fs_ap = gridfs.GridFS(self._mongodb, 'ap')
        fs_v = gridfs.GridFS(self._mongodb, 'vector_on_basis')              
        key = {'N':int(N),'k':int(k),'chi':int(i),'prec' : {"$gt": int(pprec -1) }}
        clogger.debug("key={0}".format(key))
        aps_in_mongo = self._aps.find(key).distinct('_id')
        aps_in_file = 0 
        clogger.debug("Already have {0} ap lists in mongodb! Need {1}".format(len(aps_in_mongo),num_factors))
        def insert_aps_into_mongodb(aps):
            r"""
            Insert a dictionary of aps into the mongo database.
            The format of aps is a dictionary with key/values:
            aps = {(N,k,i,d) : (E,v,meta) }
            where E*v gives the actual list of ap's and meta is a dictionary with cputime and sage version. 
            """
            for key,value in aps.iteritems():
                N,k,i,d = key
                E,v,meta = value
                clogger.debug("E={0}".format(E))
                clogger.debug("v=vector of length {0}".format(len(v)))
                clogger.debug("meta={0}".format(meta))
                fname = "gamma0-aplists-{0}".format(self._db.space_name(N,k,i))
                fname1 = "{0}-{1:0>3}-{2:0>5}".format(fname,d,pprec)
                apid = fs_ap.put(dumps( (E,v)),filename=fname1,
                                 N=int(N),k=int(k),chi=int(i),cchi=int(ci),
                                 character_galois_orbit=orbit,
                                 newform=int(d),
                                 prec = int(pprec),
                                 cputime = meta.get("cputime",""),
                                 sage_version = meta.get("version",""),
                                 ambient_id=ambient_id)
                aps_in_mongo.append(apid)
                clogger.debug("inserted aps :{0} ".format((num_factors,apid)))
                # Also insert the corresponding v separately (to protect against changes in sage)
                fnamev = "gamma0-ambient-v-{0}-{1:0>3}".format(self._db.space_name(N,k,i),d)
                clogger.debug("fnamev={0}".format(fnamev))
                vid = fs_v.put(dumps(v),filename=fnamev,
                               newform=int(d),
                               character_galois_orbit=orbit,
                               N=int(N),k=int(k),chi=int(i),cchi=int(ci),
                               prec = int(pprec),
                               sage_version = meta.get("version",""),
                               ambient_id=ambient_id)

        def insert_aps_into_filesdb(aps):
            r"""
            Insert aps (of the same format as above) into the files database.
            """
            for (N,k,i,d),(E,v,meta) in aps.itervalues():
                aplist_file = self._db.factor_aplist(N, k, i, d, False, pprec)
                apdir = join(aplist_file.split("/")[0:-1],"/")
                if not self._db.isdir(apdir):
                    self._db.makedirs(apdir)
                clogger.debug("aplist_file={0}, meta = {1}".format(aplist_file,meta))
                save((E,v), aplist_file)
                save(meta, self._computedb._db.meta(aplist_file))
                    
        # See if we have the coefficients in a file or not
        q = self._db.known("N={0} and k={1} and i={2}".format(N,k,i))
        for r in q:
            if r[4]>=pprec:
                aps_in_file = 1
                break
        clogger.debug("Have ap lists in filesdb : {0}".format(aps_in_file))

        # If we have coefficients both in mongo and files we don't do anything.
        if len(aps_in_mongo) == num_factors and aps_in_file==1:
            return aps_in_mongo
        elif len(aps_in_mongo) <> num_factors:
            if aps_in_file==0:
                if not compute:
                    return []
                ## No coefficients in either mongo or files => we compute and save if desired
                clogger.debug("Computing aplist! m={0}".format(num_factors))
                aps = self._computedb.compute_aplists(N,k,i,0,pprec,ambient=ambient,save=self._save_to_file)
                if aps == 0:
                    aps = {}
                    for d in range(num_factors):
                        aps[(N,k,i,d)] = self._db.load_aps(N,k,i,d,ambient=ambient,numc=pprec)   
            else:
                aps = {}
                for d in range(num_factors):
                    aps[(N,k,i,d)] = self._db.load_aps(N,k,i,d,ambient=ambient,numc=pprec)
            if not isinstance(aps,dict):
                clogger.critical("APS = {0}".format(aps))
            if aps == None or (isinstance(aps,dict) and len(aps.values()[0])<>3) or aps==-1:
                clogger.critical("APS: {0},{1},{2},{3} could not be computed!".format(N,k,i,d))
                return aps
            aps_in_mongo = insert_aps_into_mongodb(aps)
        elif len(aps_in_mongo) == num_factors and aps_in_file==0 and self._save_to_file:
            ### We have coefficients in mongo, need to save them to file
            aps = self.get_aps(N,k,i)
            clogger.debug("Need to insert aps into the files! num_Factors={0}".format(num_factors))
            insert_aps_into_filesdb(aps)
        return aps_in_mongo


 

 
    

    def check_records(self,nrange,krange,irange='all',check_content=False):
        r"""
        We check if the records corresponding to M(N_i,k_i,i_i) is complete or not.
        
        """
        args = []
        if not isinstance(nrange,(list,tuple)):
            nranfge = [nrange]
        if not isinstance(krange,(list,tuple)):
            krange = [krange]        
        for n in nrange:
            if irange == 'all':
                irange = range(euler_phi(n))
            elif not isinstance(irange,(list,tuple)):
                irange = [irange]
            for k in krange:                                    
                for i in irange:
                    args.append((n,k,i,check_content))
                        
        return list(self.check_record(args))
    
    @parallel(ncpus=8)        
    def check_record(self,N,k,i,check_content=False,recheck=False):
        r"""

        We check if the record corresponding to M(N,k,i) is complete or not.
        
        """
        res = {}
        s = {'N':int(N),'k':int(k),'chi':int(i)}
        res = {}
        ### Check the ambient space
        if not recheck:
            if self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(i),'complete':True}).count()>0:
                return  {'modular_symbols':True,'aps':True,'factors':True}
                
        M = self.get_ambient(N,k,i,compute=False)
        if M is None:
            res['modular_symbols']=False
            numf = 0
        else:
            if check_content and not 'modsym.ambient' in str(M.__class__):
                clogger.warning("Space is reconstructed with wrong class!")
                clogger.warning("type(M)={0}".format(type(M)))
                clogger.warning("type(MS)={0}".format(ModularSymbols(1,2).__class__))
                res['modular_symbols']=False
            else:
                res['modular_symbols']=True
            rec = self._modular_symbols.find_one({'N':int(N),'k':int(k),'chi':int(i)})
            numf = rec['orbits']
            ambient_id = rec['_id']
        ## Check factors
        numf1 = self.number_of_factors(N,k,i)
        res['factors'] = False
        facts = {}
        clogger.debug(" num facts in db={0} and in the ms record:{1}".format(numf1,numf))
        if numf1 == numf:
            res['factors'] = True
            if check_content:
                facts = self.get_factors(N,k,i)
                d = 0
                for f in facts.values():
                    d+=f.dimension()
                clogger.debug("Sum of dimensions of factors: {0}".format(d))
                #if M is None:
                if i <> 0:
                    d1 = dimension_new_cusp_forms(sage_character_galois_orbit_rep_from_number(N,i),k)
                else:
                    d1 = dimension_new_cusp_forms(N,k)
                clogger.debug("Dimension of space is: {0}".format(d1))                    
                if d <> d1:
                    res['factors'] = False
                    clogger.warning("Dimensions of all factors do not sum up to the total dimension!")
        ### Check ap's
        if check_content:
            aps = self.get_aps(N,k,i)
        else:
            newforms = self._aps.find({'N':int(N),'k':int(k),'chi':int(i)}).distinct('newform')
            aps = {(N,k,i,x) :True  for x in newforms}
            
        # Necessary for L-function computations (rounded to nearest 100).        
        pprec = 22 + int(RR(5) * RR(k) * RR(N).sqrt())
        pprec = max(pprec,100)
        pprec = ceil(RR(pprec)/RR(100))*100
        res['aps']=False
        clogger.debug("facts={0}, numf={1}".format(facts,numf))
        clogger.debug("aps.keys={0}".format(aps.keys()))
        if res['factors'] is True:
            res['aps'] = len(aps.keys())==numf
            if check_content:
                for t in facts.keys():
                    clogger.debug("t={0}".format(t))
                    apd = aps.get(t,[])
                    clogger.debug("APs[{0}]= len: {1}".format(t,len(apd)))
                    if apd == None:
                        res['aps']=False
                    else:
                        if len(apd)>0:
                            res['aps']=False
                            if check_content:
                                #clogger.debug("checking coefficients! len(apd[0])={0}".format(len(apd[0])))                                
                                for E,v,meta in apd:
                                    clogger.debug("checking coefficients! len(v)=")
                                    #a = E*v
                                    ##  Check that we have the correct number of primes.
                                    
                                    if E[0,0]<>0 and len(v)==E.ncols() and    E.nrows()>=prime_pi(pprec):
                                        res['aps'] = True
                                        break
                                    clogger.debug("done checking coeffs!")
        else:
            if len(aps.keys())==numf:
                res['aps'] = True
        if res.values().count(False)==0:
            # Record is complete so we mark it as such
            self._modular_symbols.update({'_id':ambient_id},{"$set":{'complete':True}})
        clogger.debug("done checking!")
        return res



    def find_records_needing_completion(self,nrange=[],krange=[],chi=None,check_content=False,recheck=False):
        r"""
        Check all records within a specified bound.
        """
        clogger.debug("Find records needing completion!")
        s = {}
        res = {}
        if not isinstance(krange,(list,tuple)):
            krange = [krange]
        if not isinstance(nrange,(list,tuple)):
            nrange = [nrange]            
        if nrange <> []:
            s['N'] = {"$lt":int(nrange[-1]+1), "$gt":int(nrange[0]-1)}
        if krange <> []:
            s['k'] = {"$lt":int(krange[-1]+1), "$gt":int(krange[0]-1)}
        if chi=='trivial' or chi==0:
            s['chi'] = int(0)
        args = []
        for r in self._modular_symbols.find(s):
            N = r['N']; k=r['k']; i = r['chi']
            clogger.debug("r = {0}".format((N,k,i)))
            args.append((N,k,i,check_content,recheck))

        check = self.check_record(args)
        for arg,val in check:
            if val.values().count(False)>0:
                res[arg[0:2]] = val
        return res

    def complete_records(self,nrange=[],krange=[],chi=None,ncpus=1,check_content=False):
        r"""
        Check all records within a specified bound and update / compute the incomplete ones.

        INPUT:

        - nrange -- list/tuple : give upper and lower bound for the levels we complete
        - krange -- list/tuple : give upper and lower bound for the weights we complete
        - chi    -- integer    : if not None we only look at this character (e.g. for trivial character chi=0)
        """
        recs = self.find_records_needing_completion(nrange,krange,check_content,chi=chi)
        args = []
        for N,k,i in recs.keys():
            if not chi is None and i<>chi:
                continue
            args.append((N,k,i))
        self.get_or_compute_spaces(args,ncpus=ncpus)
        return True
        
    
def precision_needed_for_L(N,k,**kwds):
    r"""
    Returns the precision (number of coefficients) needed to compute the first zero
    of the L-function (as on the LMFDB pages). This bound is taken from there and is probably heuristic.
    """
    pprec = 22 + int(RR(5) * RR(k) * RR(N).sqrt())
    pprec = max(pprec,kwds.get('pprec',100))
    ## Get even hundreds of primes to look nicer.
    return ceil(RR(pprec)/RR(100))*100
