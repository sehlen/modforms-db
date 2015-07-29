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
import os,re,sys
from string import join
import pymongo
import gridfs
from compmf.filesdb import FilenamesMFDBLoading
from compmf.compute import ComputeMFData
from compmf.character_conversions import (
    dirichlet_character_conrey_from_sage_character_number,
    dirichlet_character_conrey_galois_orbit_numbers_from_character_number,
    dirichlet_character_sage_galois_orbit_rep_from_number,
    dirichlet_character_sage_galois_orbits_reps,
    sage_character_to_sage_galois_orbit_number,
    conrey_character_number_from_sage_galois_orbit_number
)
from sage.all import nth_prime,prime_pi,parallel,loads,dimension_new_cusp_forms,RR,ceil,load,dumps,save,euler_phi,floor,QQ
from utils import are_compatible
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
        if pymongo.version_tuple[0] < 3:
            from pymongo import Connection
            _C = Connection(port=port)
            self._mongodb = pymongo.Connection('{0}:{1}'.format(host,port))[db]
            self._mongo_conn = pymongo.Connection('{0}:{1}'.format(host,port))
        else:
            from pymongo.mongo_client import MongoClient
            self._mongodb = MongoClient('{0}:{1}'.format(host,port))[db]
            self._mongo_conn = MongoClient('{0}:{1}'.format(host,port))


        ## Our databases
        self._modular_symbols_collection = 'Modular_symbols'
        self._newform_factors_collection = 'Newform_factors'
        self._aps_collection = 'ap'
        self._atkin_lehner_collection = 'Atkin_Lehner'
        self._modular_symbols = self._mongodb["{0}.files".format(self._modular_symbols_collection)]
        self._newform_factors = self._mongodb["{0}.files".format(self._newform_factors_collection)]
        self._aps = self._mongodb["{0}.files".format(self._aps_collection)]
        self._atkin_lehner = self._mongodb["{0}.files".format(self._atkin_lehner_collection)]
        self._twists = self._mongodb["twists"]
        self._file_collections = [self._modular_symbols_collection,self._newform_factors_collection,self._aps_collection,self._atkin_lehner_collection]


        ## The indices we use for the collections given above
        self._indices = {'Modular_symbols.files' :  {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING)],
                                               'unique':True,'name':"N_k_cchi"},
                         'Newform_factors.files' : {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
                                             ("newform",pymongo.ASCENDING)],
                                               'unique':True,'name':"N_k_cchi_d"},
                         'ap.files': {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
                                                ("newform",pymongo.ASCENDING),("prec",pymongo.ASCENDING)],
                                         'unique':True,'name':"N_k_cchi_d"}, 
                         'Atkin_Lehner.files' : {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
                                                    ("newform",pymongo.ASCENDING)],
                                           'unique':True,'name':"N_k_cchi_d"}
                         }
        self._sage_version = sage.version.version

    def __repr__(self):
        r"""
        String representation of self.
        """
        s="Modular forms database at mongodb: {0}".format(self._mongo_conn)
        return s

    def create_indices(self,noremove=1):
        r"""
        Creates indices for our databases
        """
        for col in self._file_collections:
            if 'files' not in col:
                col = "{0}.files".format(col)
            ix = self._indices[col]
            print "Creating index for collection ",col
            try: 
                self._mongodb[col].create_index(ix['keys'],unique=ix['unique'],name=ix['name'])
            except pymongo.errors.DuplicateKeyError:
                print "Removing duplicates!"
                keys =  [x[0] for x in ix['keys']]
                self.remove_duplicates(col,keys,dryrun=noremove)
                try: 
                    self._mongodb[col].create_index(ix['keys'],unique=ix['unique'],name=ix['name'])
                except pymongo.errors.DuplicateKeyError:
                    clogger.critical("We could not remove all duplicates!")
            # remove duplicates


    def remove_duplicates(self,col,keys,dryrun=1):
        r"""
        Remove duplicates from mongo db
        

        """
        from sage.all import deepcopy
        if 'files' not in col:
            ccol='{0}.files'.format(col)
        else:
            ccol = col
        clogger.debug("ccol={0}".format(ccol))
        clogger.debug("keys={0}".format(keys))
        flds = deepcopy(keys); flds.extend(['uploadDate','filename','_id','chi'])
        
        clogger.debug("flds={0}".format(flds))
        nnmax = max(self._mongodb[ccol].find().distinct('N'))
        args = []
        if 'ap' in col:
            step=50
        else:
            step = 100
            flds.append('prec')
        #h = RR(nnmax)/32.0
        for j in range(RR(nnmax)/RR(step)):
            nmin = j*step; nmax = (j+1)*step
            args.append((col,keys,flds,dryrun,nmin,nmax))
        return list(self.remove_duplicates32(args))
        
    @parallel(ncpus=32)
    def remove_duplicates32(self,col,keys,flds,dryrun=1,nmin=0,nmax=10000):
        r"""
        Remove duplicate records in collection col for which keys are not unique.
        """

        #keys = [x[0] for x in keys]
        fs = gridfs.GridFS(self._mongodb,col.split(".")[0])
        if 'files' not in col:
            ccol='{0}.files'.format(col)
        else:
            ccol = col
        s = {}
        clogger.debug("nmin = {0} \t nmax= {1} \t col={2} \t ccol={3}".format(nmin,nmax,col,ccol))
        #if ccol=='Newform_factors.files':
        s = {'N':{"$lt":int(nmax),"$gt":int(nmin)-1}}
        for r in self._mongodb[ccol].find(s,projection=flds).sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING),('uploadDate',pymongo.DESCENDING)]):
            id=r['_id']
            s = {}
            for k in keys:
                try:
                    s[k]=r[k]
                except KeyError as e:
                    if k=='cchi':
                        clogger.warning("rec without cchi: r={0}".format(r))
                        c = dirichlet_character_conrey_from_sage_character_number(r['N'],r['chi'])
                        ci = c.number()
                        self._mongodb[ccol].update({'_id':r['_id']},{"$set":{'cchi':ci}})
                        clogger.debug("Added cchi!")
                    #raise KeyError,e.message
            #print "s=",s
            q = self._mongodb[ccol].find(s,projection=flds).sort('uploadDate',1)
            if q.count()==1:
                continue
            for rnew in q: 
                if rnew['_id']==id:
                    continue
                clogger.debug("s = {0}".format(s))
                clogger.debug("Removing record {0} in collection {1}".format(rnew,col))
                clogger.debug("Duplicate of {0}".format(r))
                if dryrun:
                    clogger.debug("Not really deleting!")
                else:
                    clogger.debug("We are really deleting {0} from {1}".format(rnew['_id'],fs._GridFS__collection))
                    fs.delete(rnew['_id'])
        
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

    def existing_records_mongo(self,nrange=[],krange=[],complete=1):
        r"""

        Return a list of tuples (N,k,i) corresponding to (complete) records in the mongo database
        and where N is in nrange, k in krange.
        """
        s = {'complete':{"$gt":int(complete)-1}}
        if nrange <> []:
            s['N'] = {"$gt": int(nrange[0]-1), "$lt": int(nrange[1])}
        if krange <> []:
            s['k'] = {"$gt": int(krange[0]-1), "$lt": int(krange[1])}
        q = self._modular_symbols.find(s).sort('N',1).sort('k',1)
        res = [ (r['N'],r['k'],r['cchi']) for r in q]
        
        return res
        
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

        keywords:
         - compute (True) -- if True compute an ambient space if it is not in the database.
         - get_record (False) -- if True return the database record instead of the space.
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
        if kwds.get('get_record',False):
            return self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(i)})
        return self.load_from_mongo('Modular_symbols',ambient_id)

    def get_dimc(self,N,k,i):
        r"""
        Get dimension of cusp forms in S_k(N,i).

        """
        r = self._modular_symbols.find_one({'N':int(N),'k':int(k),'chi':int(i)},projection=['dimc'])
        if r is None:
            return -1
        return r.get('dimc',-1)
    
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
    
    def get_aps(self,N,k,i,d=None,character_naming='sage',prec_needed=0,coeffs=False):
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
        clogger.debug("find aps with s={0}".format(s))
        res = {}
        for r in self._aps.find(s):
            fid=r['_id']
            d = r['newform']
            prec = r['prec']
            meta = {'cputime':r.get('cputime'),'version':r.get('sage_version')}
            E,v = self.load_from_mongo('ap',fid)
            clogger.debug("id={0} and E={1}".format(fid,E))
            t = (int(N),int(k),int(i),int(d))
            if not res.has_key(t):
                res[t]={}
            if prec_needed == 0 or coeffs == False:
                res[t][prec]=(E,v,meta)
            else:
                if prec >= prec_needed and coeffs:
                    return E*v
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
        nrange = self._db.known_hosts_levels()
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
        ncpus = kwds.get('ncpus',1)
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
        clogger.debug("Compute and/or Insert {0}".format((N,k,i)))
        clogger.debug("Computing ambient modular symbols")
        if kwds.get('Nmax',0)<>0 and kwds.get('Nmax')>N:
            return 
        ambient_fid = self.compute_ambient(N,k,i,**kwds)
        # if the space is empty we return None and exit;
        if ambient_fid is None: ## If something went wrong or there isn't even Eisenstein series.
            clogger.debug("Ambient is None!")            
            return True
        if isinstance(ambient_fid,int) and ambient_fid ==1: ## If we have no cusp forms
            clogger.debug("Ambient space of cusp forms is 0 dimensional!")
            self._db.update_known_db((N,k,i,0,0))
            return True
        # Insert ambient modular symbols
        kwds['ambient_id']=ambient_fid
        clogger.debug("Getting factors!")
        factor_ids = self.compute_factors(N,k,i,**kwds)
        kwds['factor_ids']=factor_ids 
        # Necessary number of coefficients for L-function computations.
        if not factor_ids is None and factor_ids <> []: 
            pprec = precision_needed_for_L(N,k,pprec=100)
            clogger.debug("Computing aps to prec {0} for i = {1}".format(pprec,i))
            ap = self.compute_aps(N,k,i,pprec,**kwds)
            clogger.debug("Got {0}".format(ap))
            try:
                atkin_lehner = self.compute_atkin_lehner(N,k,i,**kwds)
            except:
                atkin_lehner = None
            if not (ambient_fid is None or factor_ids is None or ap is None or atkin_lehner is None):
                completeness = int(1)
            else:
                print "Factor=",factor_ids
                print "AL=",atkin_lehner
                print "ap=",ap
                print "ambient=",ambient_fid
                completeness = int(0)
        else:
            completeness = int(1)
        self._modular_symbols.update({'_id':ambient_fid},{"$set":
                                                          {'complete':completeness}})                                                              
        return True    

    def compute_atkin_lehner(self,N,k,i,**kwds):
        r"""
        Compute the Atkin-Lehner eigenvalues of space (N,k,i).
        """
        verbose = kwds.get('verbose',0)
        c = dirichlet_character_conrey_from_sage_character_number(N,i)
        ci = c.number()
        if not c.is_trivial(): #or c.multiplicative_order()==2):
            return []
        al_in_mongo = self._atkin_lehner.find({'N':int(N),'k':int(k),'chi':int(i),'cchi':int(ci)}).distinct('_id')
        fs = gridfs.GridFS(self._mongodb, 'Atkin_Lehner')
        orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
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
                fname = 'atkin_lehner_evs-{0:0>5}-{1:0>3}-{2:0>3}-{3:0>3}'.format(N,k,i,d)
                fid = fs.put(dumps(atkin_lehner),filename=fname,
                             N=int(N),k=int(k),chi=int(i),newform=int(d),cchi=int(ci),
                             character_galois_orbit = orbit,
                             cputime = meta.get("cputime",""),
                             sage_version = meta.get("version",""))
                al_in_mongo.append(fid)
        return al_in_mongo
        
        
    def compute_ambient(self,N,k,i,**kwds):
        r"""
        Compute the ambient space with the given parameters and insert it in mongo if it is not there. 
        """
        if not are_compatible(N,k,i):
            return None
        files_ms = self._modular_symbols
        fs_ms = gridfs.GridFS(self._mongodb, 'Modular_symbols')
        verbose = kwds.get('verbose',0)
        ci = dirichlet_character_conrey_from_sage_character_number(N,i).number()
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
        # Next see if we have it in the files database. If not we will compute it.
        try:
            ambient = self._db.load_ambient_space(N,k,i)
            clogger.debug("Loaded ambient={0}".format(ambient))
            ambient_in_file = 1
        except ValueError:
            clogger.debug("Could not load ambient!")
            ambient = None
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
        if ambient_in_mongo == 0 and ambient_in_file == 1: #not ambient is None:
            metaname = self._db.space(N,k,i,False)+"/ambient-meta.sobj"
            clogger.debug("metaname={0}".format(metaname))
            try:
                meta = load(metaname)
            except (OSError,IOError):
                meta={}
            fname = "gamma0-ambient-modsym-{0}".format(self._db.space_name(N,k,i).split("/")[1])
            ## Note that we have to update the number of orbits.
            if ambient is None:
                ambient = self._db.load_ambient_space(N,k,i)
            dima = int(ambient.dimension())
            dimc = int(ambient.cuspidal_submodule().dimension())
            clogger.debug("Save ambient to mongodb! ambient={0}:{1}".format((N,k,i),ambient))
            orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
            fid = fs_ms.put(dumps(ambient),filename=fname,
                            N=int(N),k=int(k),chi=int(i),orbits=int(0),
                            dima=dima,dimc=dimc,
                            character_galois_orbit=orbit,
                            cchi=int(ci),
                            cputime = meta.get("cputime",""),
                            sage_version = meta.get("version",""))
        else:
            fid = ambient_in_mongo
        clogger.debug("fid={0}".format(fid))
        if fid == 0 or ambient is None:
            fid = None
        return fid

    def compute_factors(self,N,k,i,**kwds):
        r"""
        Compute / store newform factors of parameters N,k,i
        """
        from utils import orbit_label
        verbose = kwds.get('verbose',0)
        ambient_id = kwds.get('ambient_id',None)
        compute = kwds.get('compute',self._do_computations)
        if ambient_id is None:
            ambient_id = self.compute_ambient(N,k,i,**kwds)
        ci = dirichlet_character_conrey_from_sage_character_number(N,i).number()
        if ambient_id is None:
            clogger.debug("No ambient space!")
            ambient_id = self.compute_ambient(N,k,i,**kwds)
        files_fact = self._newform_factors
        files_ms = self._modular_symbols
        fs_fact = gridfs.GridFS(self._mongodb, 'Newform_factors')
        factors_in_file = self._db.number_of_known_factors(N,k,i)
        factors_in_mongo = files_fact.find({'ambient_id':ambient_id}).distinct('_id')
        dimc = self.get_dimc(N,k,i)
        if dimc == 0:
            return []  # There are no factors to compute
        clogger.debug("factors: in_file={0} in mongo={1}".format(factors_in_file,factors_in_mongo))
        if factors_in_file == 0 and  len(factors_in_mongo)==0:
            if not compute:
                clogger.debug("No factors exist in db and we do not compute anything!")
                return None
                # Else compute and insert into the files.
            factors_in_file = self._computedb.compute_decompositions(N,k,i)
            clogger.debug("Computing factors! m={0}".format(factors_in_file))

        if len(factors_in_mongo)==0:
            fname = "gamma0-factors-{0}".format(self._db.space_name(N,k,i))
            clogger.debug("Inserting factors into mongo! fname={0}".format(fname))
            num_factors_in_file = self._db.number_of_known_factors(N,k,i)
            ambient = self.get_ambient(N,k,i,ambient_id=ambient_id)
            orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)

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
                label = orbit_label(d)
                facid = fs_fact.put(dumps(factor),filename=fname1,
                                    N=int(N),k=int(k),chi=int(i),
                                    cchi=int(ci),
                                    character_galois_orbit=orbit,
                                    newform=int(d),
                                    cputime = meta.get("cputime",""),
                                    sage_version = meta.get("version",""),
                                    ambient_id=ambient_id,
                                    hecke_orbit_label='{0}.{1}.{2}{3}'.format(N,k,ci,label),
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
        from utils import orbit_label
        if pprec is None:
            pprec = precision_needed_for_L(N,k)
        pprec = int(pprec)
        ambient_id = kwds.get('ambient_id',None)
        if ambient_id is None:
            ambient_id = self.compute_ambient(N,k,i,**kwds)
            kwds['ambient_id']=ambient_id
        ambient = self.get_ambient(N,k,i,ambient_id=ambient_id,verbose=0)
        num_factors = len(kwds.get('factors_ids',self.compute_factors(N,k,i,**kwds)))

        compute = kwds.get('compute',self._do_computations)
        verbose = kwds.get('verbose')
        c = dirichlet_character_conrey_from_sage_character_number(N,i)
        ci = c.number()        
        orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
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
            
            if not isinstance(aps,dict):
                clogger.warning("Trying to insert non-dict aps:{0}".format(aps))
                return
            res = []
            for key,value in aps.iteritems():
                N,k,i,d = key
                E,v,meta = value
                if isinstance(E,tuple):
                    E,v = E
                clogger.debug("E={0}".format(E))
                clogger.debug("v=vector of length {0}".format(len(v)))
                clogger.debug("meta={0}".format(meta))
                fname = "gamma0-aplists-{0}".format(self._db.space_name(N,k,i))
                fname1 = "{0}-{1:0>3}-{2:0>5}".format(fname,d,pprec)
                label = orbit_label(d)
                clogger.debug("label={0}".format((N,k,i,d)))
                # delete if exists
                r = self._aps.find_one({'filename':fname1})
                fs_ap.delete(r['_id'])
                apid = fs_ap.put(dumps( (E,v)),filename=fname1,
                                 N=int(N),k=int(k),chi=int(i),cchi=int(ci),
                                 character_galois_orbit=orbit,
                                 newform=int(d),
                                 prec = int(pprec),
                                 cputime = meta.get("cputime",""),
                                 sage_version = meta.get("version",""),
                                 hecke_orbit_label='{0}.{1}.{2}{3}'.format(N,k,ci,label),
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
                res.append(vid)
            return res
        def insert_aps_into_filesdb(aps):
            r"""
            Insert aps (of the same format as above) into the files database.
            """
            res = []
            for key,val in aps.iteritems():
                clogger.debug("key={0}".format(key))
                clogger.debug("type(val)={0}".format(type(val)))
                if isinstance(val,dict):
                    clogger.debug("val.keys()={0}".format(val.keys()))
                N,k,i,d = key
                for prec in val.keys():
                    E,v,meta = val[prec]
                    aplist_file = self._db.factor_aplist(N, k, i, d, False, prec)
                    apdir = join(aplist_file.split("/")[0:-1],"/")
                    if not self._db.isdir(apdir):
                        self._db.makedirs(apdir)
                    clogger.debug("aplist_file={0}, meta = {1}".format(aplist_file,meta))
                    save(E, aplist_file)
                    
                    vname = self._db.factor_dual_eigenvector(N, k, i, d)
                    if not self._db.path_exists(vname):
                        save(v,vname)
                    save(meta, self._db.meta(aplist_file))
                    res.append(vname)
            return res
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
                        E,v,meta  = self._db.load_aps(N,k,i,d,ambient=ambient,numc=pprec)   
                        aps[(N,k,i,d)] = E,v,meta
            else:
                aps = {}
                for d in range(num_factors):
                    E,v,meta = self._db.load_aps(N,k,i,d,ambient=ambient,numc=pprec)
                    aps[(N,k,i,d)] = E,v,meta
                    
            if not isinstance(aps,dict):
                clogger.critical("APS = {0}".format(aps))
            if aps == None or (isinstance(aps,dict) and len(aps.values()[0])<>3) or aps==-1:
                clogger.critical("APS: {0},{1},{2},{3} could not be computed!".format(N,k,i,d))
                return aps
            return insert_aps_into_mongodb(aps)
        elif len(aps_in_mongo) == num_factors and aps_in_file==0 and self._save_to_file:
            ### We have coefficients in mongo, need to save them to file
            aps = self.get_aps(N,k,i)
            clogger.debug("Need to insert aps into the files! num_Factors={0}".format(num_factors))
            return insert_aps_into_filesdb(aps)
        clogger.critical("aps for: {0},{1},{2},{3} could not be computed!".format(N,k,i,d))
        return aps_in_mongo


 

 
    

    def check_records(self,nrange,krange,irange='all',check_content=False,recheck=False,ncpus=8):
        r"""
        We check if the records corresponding to M(N_i,k_i,i_i) is complete or not.
        
        """
        args = []
        if not isinstance(nrange,(list,tuple)):
            nranfge = [nrange]
        if not isinstance(krange,(list,tuple)):
            krange = [krange]
        if check_content:
            check_level = 2
        else:
            check_level = 1
        clogger.debug("check_record with level: {0}".format(check_level))
        s = {'projection':['N','k','chi']}
        if recheck is False:
            s['complete']={"$lt":check_level+int(1)}
        clogger.debug("s = {0}".format(s))
        for r in self._modular_symbols.find(s).sort([('N',pymongo.ASCENDING),('chi',pymongo.ASCENDING),('k',pymongo.ASCENDING)]):
            n = r['N']; k=r['k']; chi=r['chi']
            if nrange<>[] and (n < nrange[0] or n > nrange[-1]):
                continue
            if krange<> [] and (k < krange[-1] or k > krange[0]):
                continue
            if irange <> 'all' and irange<>[] and (chi < irange[0] or chi >irange[-1]):
                continue
            args.append((n,k,chi,check_content,recheck))
        clogger.debug("args={0}".format(args))
        if ncpus >= 32:
            return list(self.check_record32(args))
        elif ncpus >= 16:
            return list(self.check_record16(args))
        elif ncpus >= 8:
            return list(self.check_record8(args))                    
        return list(self.check_record(args))
    
    @parallel(ncpus=8,verbose=True)           
    def check_record8(self,N,k,i,check_content=False,recheck=False):
        return self.check_record(N,k,i,check_content,recheck)
    
    @parallel(ncpus=16,verbose=True)                   
    def check_record16(self,N,k,i,check_content=False,recheck=False):
        return self.check_record(N,k,i,check_content,recheck)
    @parallel(ncpus=32,verbose=True)            
    def check_record32(self,N,k,i,check_content=False,recheck=False):
        return self.check_record(N,k,i,check_content,recheck)
        
    @parallel(ncpus=1,verbose=True)            
    def check_record(self,N,k,i,check_content=False,recheck=False):
        r"""

        We check if the record corresponding to M(N,k,i) is complete or not.
        
        """
        
        res = {}
        s = {'N':int(N),'k':int(k),'chi':int(i)}
        res = {}
        ### Check the ambient space
        check_level = int(3) if check_content else int(1)
        if not recheck:
            if self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(i),'complete':{"$gt":check_level-int(1)}}).count()>0:
                return  {'modular_symbols':True,'aps':True,'factors':True}
        clogger.debug("Checking N,k,i={0}".format((N,k,i)))
        ambient_id = None
        M = self.get_ambient(N,k,i,compute=False)
        if M is None:
            res['modular_symbols']=False
            numf = 0
        else:
            self.check_character(N,k,i)
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
        clogger.debug("will check number of factors! for M={0}".format(M))
        numf1 = self.number_of_factors(N,k,i)
        facts = {}
        clogger.debug(" num facts in db={0} and in the ms record:{1}".format(numf1,numf))
        if numf1 == 0 and numf == 0:
            self._modular_symbols.update({'_id':ambient_id},{"$set":{'complete':int(3)}})
            return res
        res['factors'] = False
        d = -1; d1= -1
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
                    d1 = dimension_new_cusp_forms(dirichlet_character_sage_galois_orbit_rep_from_number(N,i),k)
                else:
                    d1 = dimension_new_cusp_forms(N,k)
                clogger.debug("Dimension of space is: {0}".format(d1))                    
                if d <> d1:
                    res['factors'] = False
                    clogger.warning("Dimensions of all factors do not sum up to the total dimension! n,k,chi={0}".format((N,k,i)))
                if not M is None:
                    dim = M.dimension()
                    self._modular_symbols.update({'_id':ambient_id},{"$set":{'dim_new_cusp':int(d1),'dimension':int(dim)}})            
        ### Check ap's
        # Necessary for L-function computations (rounded to nearest 100).        
        pprec = 22 + int(RR(5) * RR(k) * RR(N).sqrt())
        pprec = max(pprec,100)
        pprec = ceil(RR(pprec)/RR(100))*100
        #if check_content:
        #    #aps = self.get_aps(N,k,i)
        #else:
        #    newforms = self._aps.find({'N':int(N),'k':int(k),'chi':int(i)}).distinct('newform')
        #    #aps = {(N,k,i,x) : {pprec: True}  for x in newforms}
        if d==d1 and d==0:
            res['aps']=True
        else:
            res['aps']=False
        if not check_content and res['factors']:
            newforms_with_aps = self._aps.find({'N':int(N),'k':int(k),'chi':int(i)}).distinct('newform')
            res['aps'] = len(newforms_with_aps) == numf
        clogger.debug("facts={0}, numf={1}".format(facts,numf))
        fs_ap = gridfs.GridFS(self._mongodb,self._aps_collection)
        if res['factors'] is False:
            facts = {}
        for t in facts.keys():
            clogger.debug("t={0}".format(t))
            N,k,i,d=t
            s = {'N':int(N),'k':int(k),'chi':int(i),'newform':int(d)}
            q = self._aps.find(s)
            if q.count()==0:
                res['aps']=False
                break
            precs = []
            clogger.debug("Number of coefficient records of this t={0}".format(q.count()))
            for r in q:
                id =r['_id']; prec=r['prec']
                E,v = loads(fs_ap.get(id).read())                    
                #clogger.debug("type(E)={0}".format(type(E)))
                clogger.debug("rec={0}".format(r))
                if isinstance(E,tuple):
                    print "r=",r
                    raise ValueError,"Wrong format of E!"
                res['aps']=False
                clogger.debug("checking coefficients! len(v)={0} E.nrows={1}, E.ncols={2}, E[0,0]==0:{3}, pi(pprec)={4} assumed prec={5}".format(len(v),E.nrows(),E.ncols(),E[0,0] is 0,prime_pi(pprec),prec))
                
                nprimes_in_db = E.nrows()
                nprimes_assumed = prime_pi(prec)
                prec_in_db = int(nth_prime(nprimes_in_db+1)-1) # observe that we can get all coefficients up to the next prime - 1
                precs.append(prec_in_db)
                clogger.debug("Precision: claimed ={0}  actually in db: {1}".format(prec,prec_in_db))
                clogger.debug("Nprimes assumed: {0} in db: {1}".format(nprimes_assumed,nprimes_in_db))
                if nprimes_in_db <> nprimes_assumed:  ### The coefficients in the database are not as many as assumed!
                    clogger.debug("Have {0} aps in the database and we claim that we have {1}".format(E.nrows(),prime_pi(prec)))
                    #int(ceil(RR(nth_prime(E.nrows()))/RR(100))*100)
                    if  self._aps.find({'N':r['N'],'k':r['k'],'chi':r['chi'],'newform':r['newform'],'prec':prec_in_db}).count()>0:
                        clogger.debug("We already have this prec in the database so we remove this record!")
                        self._aps.remove({'_id':r['_id']})
                    else:
                        fname = "gamma0-aplists-{0:0>5}-{1:0>3}-{2:0>3}-{3:0>3}-{4:0>3}".format(N,k,i,d,prec_in_db)
                        clogger.debug("updating record {0} with prec = prec_in_db={1}".format(id,prec_in_db))                                  
                        q = self._aps.update({'_id':id},
                                         {"$set":{'prec':prec_in_db,'filename':fname}},w=int(1),multi=False,upsert=True)
                        clogger.debug("Updated : {0}".format(q))
                    ## Also check the file system:
                    
                    ##  We now check that E,v is consistent: i.e. E is non-zero E*v exists and that we have the correct number of primes.
                if (not (E[0,0] is 0)) and len(v)==E.ncols() and  prec_in_db >= prec:
                    res['aps'] = True
            maxprec = max(precs)
            if maxprec < pprec:
                clogger.debug("have coefficients but not sufficiently many! Need {0} and got {1}".format(pprec,maxprec))
                res['aps'] = False
            clogger.debug("done checking coeffs! t={0}".format(t))
        if res.values().count(False)==0:
            # Record is complete so we mark it as such
            self._modular_symbols.update({'_id':ambient_id},{"$set":{'complete':check_level}})
        elif ambient_id is not None:
            self._modular_symbols.update({'_id':ambient_id},{"$set":{'complete':int(0)}})
        else:
            res['aps']=False
            res['factors']=False
            res['modular_symbols']=False
            clogger.debug("Space {0},{1},{2} is not in database".format(N,k,i))
        clogger.debug("done checking record!")
        sys.stdout.flush()
        return res



    def find_records_needing_completion(self,nrange=[],krange=[],chi=None,check_content=False,recheck=False,ncpus=1):
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
        elif isinstance(chi,list):
            s['chi'] = {"$in":map(int,chi)}
        args = []
        clogger.debug("search  pattern :{0}".format(s))
        for r in self._modular_symbols.find(s):
            N = r['N']; k=r['k']; i = r['chi']
            clogger.debug("r = {0}".format((N,k,i)))
            args.append((N,k,i,check_content,recheck))
#        clogger.debug("args={0}".format(args))
        if ncpus >= 32:
            check = list(self.check_record32(args))
        elif ncpus >= 16:
            check = list(self.check_record16(args))
        elif ncpus >= 8:
            check = list(self.check_record8(args))                    
        else:
            check = list(self.check_record(args))
#        clogger.debug("check={0}".format(check))            
        #check = self.check_record(args)
        for arg,val in check:
            try: 
                if val.values().count(False)>0:
                    res[arg[0][0:3]] = val
            except AttributeError:
                clogger.critical("arg = {0} val={1}".format(arg,val))
        ## Then add the spaces not yet in the database
        for n in nrange:
            norbits=len(dirichlet_character_sage_galois_orbits_reps(n))
            for k in krange:
                for i in range(norbits):
                    if self._modular_symbols.find({'N':int(n),'k':int(k),'chi':int(i)}).count()==0:
                        res[(n,k,i)]=[False]

        return res

    def complete_records(self,nrange=[],krange=[],chi=None,ncpus=1,check_content=False,recheck=False):
        r"""
        Check all records within a specified bound and update / compute the incomplete ones.

        INPUT:

        - nrange -- list/tuple : give upper and lower bound for the levels we complete
        - krange -- list/tuple : give upper and lower bound for the weights we complete
        - chi    -- integer    : if not None we only look at this character (e.g. for trivial character chi=0)

        - 'check_content' -- bool : set to True to make a more detailed check of completeness of records
        - 'recheck' -- bool: set to True if you want to recheck records already marked as complete.
        """
        recs = self.find_records_needing_completion(nrange,krange,chi=chi,check_content=check_content,recheck=recheck,ncpus=ncpus)
        args = []
        for N,k,i in recs.keys():
            if k==1:
                clogger.debug("Weight 1 is not implemented!")
                continue
            if not chi is None and i not in chi:
                continue
            if not are_compatible(N,k,i):
                clogger.debug("N,k,i={0} is incompatible!".format((N,k,i)))
                continue
            args.append((N,k,i))
        clogger.debug("Completing {0} spaces!".format(len(args)))
        self.get_or_compute_spaces(args,ncpus=ncpus)
        return True

    
    def check_character(self,N,k,chi,remove=1,files_separately=0):
        #if N % 10 == 0:
        clogger.debug("Checking N={0}, k={1}, chi={2}".format(N,k,chi))
        
        problems=[]
        i = 0 
        for r in self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(chi)}):
            if i>1:
                clogger.warning("Multiple records for {0}!".format((N,k,chi)))
            i=i+1
            sage.modular.modsym.modsym.ModularSymbols_clear_cache()
            ms = gridfs.GridFS(self._mongodb,'Modular_symbols')
            aps = gridfs.GridFS(self._mongodb,'aps')
            factors = gridfs.GridFS(self._mongodb,'Newform_factors')
            al = gridfs.GridFS(self._mongodb,'Atkin_Lehner')
            id = r['_id']
            cchi = r.get('cchi')
            clogger.debug("checking record: {0}".format(r))
            if cchi is None:
                clogger.debug("We don't have a Conrey character for r={0}".format(r))
                c = dirichlet_character_conrey_from_sage_character_number(N,chi)
                cchi = c.number()
                for col in self._file_collections:
                    if col=='Modular_symbols':
                        self._mongodb['{0}.files'.format(col)].update({'_id':id},{"$set":{'cchi':cchi}})
                    else:
                        self._mongodb['{0}.files'.format(col)].update({'ambient_id':id},{"$set":{'cchi':cchi}})            
            #clogger.debug("Get Modular symbols from Mongo! col={0} id={1}".format(self._modular_symbols_collection,id))
            M = self.load_from_mongo(self._modular_symbols_collection,id)
            #clogger.debug("Got Modular symbols from Mongo!")
            x = M.character()
            if N == 1:
                si = 0
                ci = 1
            else:
                si = sage_character_to_sage_galois_orbit_number(x)
                ci = conrey_character_number_from_sage_galois_orbit_number(N,si)
            if si <> chi or ci<>cchi:
                problems.append((N,k,chi))
                clogger.debug("Chi is wrong! Should be {0} not {1}".format(si,chi))
                clogger.debug("CChi is wrong! Should be {0} not {1}".format(ci,cchi))
                clogger.debug("r={0}".format(r))
                raise ArithmeticError()
                if remove == 1:
                    # First delete from mongo
                    ms.delete(id)
                    # delete aps
                    for rf in  self._aps.find({'_ambient_id':id},projection=['_id']):
                        aps.delete(rf['_id'])
                    # delete factors
                    for rf in  self._newform_factors.find({'_ambient_id':id},projection=['_id','newform']):
                        fid = rf['_id']
                        factors.delete(fid) # Delete factor
                    #    dname = self._db.factor(N,k,chi,rf['newform'])
                    #    for fname in self._db.listdir(dname):
                    #        self._db.delete_file(fname)
                    #aname = self._db.ambient(N,k,chi)
                    #dname = join(s.split("/")[0:-1],"/")
                    #for fname in self._db.listdir(dname):
                    #    self._db.delete_file(fname)
                    #os.removedirs(aname)
            else:
                clogger.debug("Characters are correct!")
            #    r = self._modular_symbols.find_one({'_id':id})
            #    #if r.get('complete') is None or r.get('complete')<2:
            #self.check_record(N,k,chi,check_content=True)
            #    self._modular_symbols.update({'_id':id},{"$set":{'complete':int(3)}})
            s = "N={0} and k={1} and i={2}".format(N,k,chi)
            clogger.debug("searching files for: {0}".format(s))
            
            for N1,k1,chi1,d,prec in self._db.known(s):
                #if N < minn or N>maxn or k<mink or k>maxk:
                #    continue
                clogger.debug("  in file with : {0}".format((N1,k1,chi1,d,prec)))
                sage.modular.modsym.modsym.ModularSymbols_clear_cache()
                if not files_separately==1 or (N1,k1,chi1) in problems:
                    continue
                if not self._db.path_exists(self._db.ambient(N1,k1,chi1)):
                    continue
                M = self._db.load_ambient_space(N1,k1,chi1)
                x = M.character()
                if N1 == 1:
                    si = 0
                else:
                    si = sage_character_to_sage_galois_orbit_number(x)
                clogger.debug("si={0}, chi={1}".format(si,chi))
                if si <> chi:
                    if (N1,k1,chi1) not in problems:
                        problems.append((N1,k1,chi1))
                    clogger.debug("in file: Chi is wrong! Got: {0} Should be {1}".format(si,chi))
                    if remove == 1:
                        dname = self._db.factor(N1,k1,chi1,d)
                        for fname in self._db.listdir(dname):
                            self._db.delete_file(fname)
                        aname = self._db.ambient(N1,k1,chi1)
                        dname = join(aname.split("/")[0:-1],"/")
                        for fname in self._db.listdir(dname):
                            self._db.delete_file(fname)
                        os.removedirs(dname)
                        clogger.debug("removed directory {0}".format(dname))
            if remove == 1 and problems<>[]:
                clogger.debug("Removed {0} records!".format(len(problems)))
        clogger.debug("Finished checking  character: N={0}, k={1}, chi={2}".format(N,k,chi))
        return problems            



    def check_if_twist(self,N=1,k=2,ci=0,d=0,fid=None):
        r"""
        Check if the newform factor given by factor_id is a twist of another form.

        Recall that if chi is a character modulo Q and F has level N then
        F_chi has level lcm(N,Q^2)

        """
        from sage.all import is_squarefree,ZZ,divisors,lcm,sturm_bound,primes
        from character_conversions import sage_character_from_number,conrey_character_from_number,dirichlet_group_sage,sage_character_to_number
        if not fid is None:
            q = self._newform_factors.find_one({'_id':factor_id})
        else:
            q = self._newform_factors.find_one({'N':int(N),'k':int(k),'chi':int(ci),'newform':int(d)})
        if q is None:
            raise ValueError,"This newform is not in the database!"
        
        N=q['N']; k=q['k']; cchi=q['cchi']; fid = q['_id']; ci=q['chi']
        d = q['newform']
            
        chi = sage_character_from_number(N,ci)
        DG = dirichlet_group_sage(N)
        is_twist = True
        possible_twists = []
        if is_squarefree(N): ## Twists always have levels with square factors.
            is_twist = False
            ### TODO: make a similar check for the character.
        else:
            #x = conrey_char
            ### Possible conductors of character we can twist with to get here
            for q in ZZ(N).divisors():
                if q == 1:
                    continue
                if not (q**2).divides(N):
                    continue
                divs = ZZ(N).divisors()
                divs.sort()
                for M in divs:
                    # Can we twist a form of level Q to get the one we have?
                    if not lcm(M,q**2).divides(N):
                        continue
                    print "possible M,q=",M,q,lcm(M,q**2)
                    for x in dirichlet_group_sage(M):
                        if x(-1) <> (-1)**k:
                            continue
                        m,xi = sage_character_to_number(x) # get number for later reference
                        ## Now check the characters
                        for y in dirichlet_group_sage(q):
                            if y.is_trivial():
                                continue
                            if DG(x)*DG(y**2)==chi:
                                print "possible character=",y
                                q,yi = sage_character_to_number(y)
                                possible_twists.append((M,xi,q,yi)) # it is possible that we have a twist of F \in M_k(M,x) with character y of modulus q.
        bd = max(sturm_bound(N,k),2)
        print "Sturm bound=",bd
        aps = self.get_aps(N,k,ci,d,prec_needed=bd,coeffs=True)
        K0 = aps[0].base_ring()
        for M,xi,q,yi in possible_twists:
            print "Checking:",M,xi,q,yi
            y = sage_character_from_number(q,yi)
            for d1 in range(self.number_of_factors(M,k,xi)):
                apsf = self.get_aps(M,k,xi,d1,prec_needed=bd,coeffs=True)
                print "aps=",apsf
                twisted_aps = [ apsf[prime_pi(p)-1]*y(p) for p in primes(bd+1)]
                K1 = twisted_aps[0].base_ring()
                if K0 == QQ:
                    K = K1
                elif K1 == QQ:
                    K = K0
                else:
                    K = K1.extend(K0.polymomial)
                print "twisted_aps=",twisted_aps
                print "K=",K
                test = [K(twisted_aps[i]) - K(aps[i]) for i in range(len(twisted_aps))]
                print "test=",test
                if test.count(0)==len(test):
                    print "We have a twist of {0} by {1}!".format((M,k,xi,d1),(q,yi))
                    return (M,k,xi,d1),(q,yi)
        return []

    def check_all_twists(self,verbose=0):
#        fs_twist = gridfs.GridFS(self._mongodb,self._twists_collection)
        for r in self._newform_factors.find().sort('N'):
            N=r['N']; k=r['k']; i=r['chi']; d=r['newform']
            if is_squarefree(N):
                continue
            t = self.check_if_twist(N,k,i,d)
            if len(t)>0:
                labela = utils.newform_label(N,k,i,d)
                labelb = utils.newform_label(t[0][0],t[0][1],t[0][2],t[0][3])
                clabel = "{0}.{1}".format(t[1][0],t[1][1])
                twist_rec = {'hecke_orbit_label':labela,'twist_of':labelb,'by_char':clabel}
                if verbose>0:
                    print "inserting; ",twist_rec
                self._twists.update(twist_rec)
                


                
def precision_needed_for_L(N,k,**kwds):
    r"""
    Returns the precision (number of coefficients) needed to compute the first zero
    of the L-function (as on the LMFDB pages). This bound is taken from there and is probably heuristic.
    """
    pprec = 22 + int(RR(5) * RR(k) * RR(N).sqrt())
    pprec = max(pprec,kwds.get('pprec',100))
    ## Get even hundreds of primes to look nicer.
    return ceil(RR(pprec)/RR(100))*100


