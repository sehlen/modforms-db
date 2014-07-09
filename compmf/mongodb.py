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
from string import join
import pymongo
import gridfs
from compmf.filesdb import FilenamesMFDBLoading
from compmf.compute import ComputeMFData
from compmf.character_conversions import dirichlet_character_conrey_from_sage_character_number,dirichlet_character_conrey_galois_orbit_numbers_from_character_number,dirichlet_character_sage_galois_orbit_rep_from_number
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

        self._file_collections = [self._modular_symbols_collection,self._newform_factors_collection,self._aps_collection,self._atkin_lehner_collection]


        ## The indices we use for the collections given above
        self._indices = {'Modular_symbols.files' :  {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING)],
                                               'unique':True,'name':"N_k_cchi"},
                         'Newform_factors.files' : {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
                                             ("newform",pymongo.ASCENDING)],
                                               'unique':True,'name':"N_k_cchi_d"},
                         'ap': {'keys.files':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
                                                  ("newform",pymongo.ASCENDING)],
                                         'unique':True,'name':"N_k_cchi_d"}, 
                         'Atkin_lehner.files' : {'keys':  [("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
                                                    ("newform",pymongo.ASCENDING),("prec",pymongo.ASCENDING)],
                                           'unique':True,'name':"N_k_cchi_d"}
                         }
        self._sage_version = sage.version.version

    def __repr__(self):
        r"""
        String representation of self.
        """
        s="Modular forms database at mongodb: {0}".format(self._mongo_conn)
        return s

    def create_indices(self):
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
                self.remove_duplicates(col,keys)
            # remove duplicates
            

    def remove_duplicates(self,col,keys,dryrun=1):
        r"""
        Remove duplicate records in collection col for which keys are not unique.
        """
        from sage.all import deepcopy
        #keys = [x[0] for x in keys]
        print "keys=",keys
        fs = gridfs.GridFS(self._mongodb,col)
        flds = deepcopy(keys); flds.extend(['uploadDate','filename','_id','chi'])
        print "flds=",flds
        if 'files' not in col:
            ccol='{0}.files'.format(col)
        else:
            ccol = col
        print "ccol=",ccol
        for r in self._mongodb[ccol].find({},fields=flds).sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING),('uploadDate',pymongo.DESCENDING)]):
            id=r['_id']
            s = {}
            for k in keys:
                s[k]=r[k]
            #print "s=",s
            for rnew in self._mongodb[col].find(s,fields=['uploadDate','filename','_id','N','k','chi','cchi']).sort('uploadDate',1):
                if rnew['_id']==id:
                    continue
                if dryrun:
                    print "Removing record {0} in collection {1}".format(rnew,col)
                    print "Duplicate of {0}".format(r)
                else:
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

    def existing_records_mongo(self,nrange=[],krange=[],complete=True):
        r"""

        Return a list of tuples (N,k,i) corresponding to (complete) records in the mongo database
        and where N is in nrange, k in krange.
        """
        s = {'complete':{"$gt":int(0)}}
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
        clogger.debug("find aps with s={0}".format(s))
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
                fname = 'atkin_lehner_evs-{0:>5}-{1:>3}-{2:>3}-{3:>3}'.format(N,k,i,d)
                fid = fs.put(dumps(atkin_lehner),filename=fname,
                             N=int(N),k=int(k),chi=int(i),newform=int(d),cchi=int(ci),
                             character_galois_orbit = orbit,
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
            orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
            fid = fs_ms.put(dumps(ambient),filename=fname,
                            N=int(N),k=int(k),chi=int(i),orbits=0,
                            character_galois_orbit=orbit,
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
        ci = dirichlet_character_conrey_from_sage_character_number(N,i).number()
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
            orbit = dirichlet_character_conrey_galois_orbit_numbers(N,ci)

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
                                    character_galois_orbit=orbit,
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
        c = dirichlet_character_conrey_from_sage_character_number(N,i)
        ci = c.number()        
        orbit = dirichlet_character_conrey_galois_orbit_numbers(N,xi)
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
            for key,val in aps.iteritems():
                clogger.debug("key={0}".format(key))
                clogger.debug("val={0}".format(val))
                N,k,i,d = key
                E,v,meta = val[0]
                aplist_file = self._db.factor_aplist(N, k, i, d, False, pprec)
                apdir = join(aplist_file.split("/")[0:-1],"/")
                if not self._db.isdir(apdir):
                    self._db.makedirs(apdir)
                clogger.debug("aplist_file={0}, meta = {1}".format(aplist_file,meta))
                save((E,v), aplist_file)
                save(meta, self._db.meta(aplist_file))
                    
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
            aps_in_mongo = insert_aps_into_mongodb(aps)
        elif len(aps_in_mongo) == num_factors and aps_in_file==0 and self._save_to_file:
            ### We have coefficients in mongo, need to save them to file
            aps = self.get_aps(N,k,i)
            clogger.debug("Need to insert aps into the files! num_Factors={0}".format(num_factors))
            insert_aps_into_filesdb(aps)
        return aps_in_mongo


 

 
    

    def check_records(self,nrange,krange,irange='all',check_content=False,ncpus=8):
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
        for r in self._modular_symbols.find({'complete':{"$lt":check_level+int(1)}},fields=['N','k','chi']).sort([('N',pymongo.ASCENDING),('chi',pymongo.ASCENDING),('k',pymongo.ASCENDING)]):
            n = r['N']; k=r['k']; chi=r['chi']
            if nrange<>[] and (n < nrange[0] or n > nrange[-1]):
                continue
            if krange<> [] and (k < krange[-1] or k > krange[0]):
                continue
            if irange <> 'all' and irange<>[] and (chi < irange[0] or chi >irange[-1]):
                continue
            args.append((n,k,chi,check_content))
        clogger.debug("args={0}".format(args))
        if ncpus >= 32:
            return list(self.check_record32(args))
        elif ncpus >= 16:
            return list(self.check_record16(args))
        elif ncpus >= 8:
            return list(self.check_record8(args))                    
        return self.check_record(args)
    
    @parallel(ncpus=8)        
    def check_record8(self,N,k,i,check_content=False,recheck=False):
        return self.check_record(N,k,i,check_content,recheck)
    @parallel(ncpus=16)        
    def check_record16(self,N,k,i,check_content=False,recheck=False):
        return self.check_record(N,k,i,check_content,recheck)
    @parallel(ncpus=32)        
    def check_record32(self,N,k,i,check_content=False,recheck=False):
        return self.check_record(N,k,i,check_content,recheck)
        
    @parallel(ncpus=1)            
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
        ### Might as well check the character here as well.
        #x = M.character()
        #si = sage_character_to_galois_orbit_number(x)
        #if si <> i:
        #    clogger.warning("Character for this record is wrong!")
        M = self.get_ambient(N,k,i,compute=False)
        clogger.debug("Checking N,k,i={0}".format((N,k,i)))
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
                    d1 = dimension_new_cusp_forms(dirichlet_character_sage_galois_orbit_rep_from_number(N,i),k)
                else:
                    d1 = dimension_new_cusp_forms(N,k)
                clogger.debug("Dimension of space is: {0}".format(d1))                    
                if d <> d1:
                    res['factors'] = False
                    clogger.warning("Dimensions of all factors do not sum up to the total dimension! n,k,chi={0}".format((N,k,i)))
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
                                    clogger.debug("checking coefficients! len(v)={0} E.nrows={1}, E.ncols={2}, E[0,0]={3}".format(len(v),E.ncols(),E.nrows(),E[0,0] is 0)
                                    #a = E*v
                                    ##  Check that we have the correct number of primes.
                                    
                                    if (not (E[0,0] is 0)) and len(v)==E.ncols() and   E.nrows()>=prime_pi(pprec):
                                        res['aps'] = True
                                        break
                                    clogger.debug("done checking coeffs!")
        else:
            if len(aps.keys())==numf:
                res['aps'] = True
        if res.values().count(False)==0:
            # Record is complete so we mark it as such
            self._modular_symbols.update({'_id':ambient_id},{"$set":{'complete':check_level}})
        clogger.debug("done checking!")
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
        args = []
        for r in self._modular_symbols.find(s):
            N = r['N']; k=r['k']; i = r['chi']
            #clogger.debug("r = {0}".format((N,k,i)))
            args.append((N,k,i,check_content,recheck))
        if ncpus >= 32:
            check = list(self.check_record32(args))
        elif ncpus >= 16:
            check = list(self.check_record16(args))
        elif ncpus >= 8:
            check = list(self.check_record8(args))                    
        else:
            check = list(self.check_record(args))
        #check = self.check_record(args)
        for arg,val in check:
            if val.values().count(False)>0:
                res[arg[0][0:3]] = val
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
            if not chi is None and i<>chi:
                continue
            args.append((N,k,i))
        clogger.debug("Completing {0} spaces!".format(len(args)))
        self.get_or_compute_spaces(args,ncpus=ncpus)
        return True

    
    def check_character(self,N,k,chi,remove=1,files_separately=0):
        from compmf.character_conversions import sage_character_to_galois_orbit_number,conrey_character_number_from_sage_galois_orbit_number
        if N % 10 == 0:
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
            al = gridfs.GridFS(self._mongodb,'Atkin-Lehner')
            id = r['_id']
            cchi = r.get('cchi')
            if cchi is None:
                print "We don't have a Conrey character for r={0}".format(r)
                c = dirichlet_character_conrey_from_sage_character_number(N,chi)
                cchi = c.number()
                for col in self._file_collections:
                    if col=='Modular_symbols':
                        self._mongodb['{0}.files'.format(col)].update({'_id':id},{"$set":{'cchi':cchi}})
                    else:
                        self._mongodb['{0}.files'.format(col)].update({'ambient_id':id},{"$set":{'cchi':cchi}})            
            M = self.load_from_mongo(self._modular_symbols_collection,id)
            x = M.character()
            if N == 1:
                si = 0
                ci = 1
            else:
                si = sage_character_to_galois_orbit_number(x)
                ci = conrey_character_number_from_sage_galois_orbit_number(N,si)
            if si <> chi or ci<>cchi:
                problems.append((N,k,chi))
                print "Chi is wrong! Should be {0}".format(si)
                print "CChi is wrong! Should be {0}".format(ci)
                print "r=",r
                if remove == 1:
                    # First delete from mongo
                    ms.delete(id)
                    # delete aps
                    for rf in  self._aps.find({'_ambient_id':id},fields=['_id']):
                        aps.delete(rf['_id'])
                    # delete factors
                    for rf in  self._newform_factors.find({'_ambient_id':id},fields=['_id','newform']):
                        fid = rf['_id']
                        factors.delete(fid) # Delete factor
                        dname = self._db.factor(N,k,chi,rf['newform'])
                        for fname in self._db.listdir(dname):
                            self._db.delete_file(fname)
                    aname = self._db.ambient(N,k,chi)
                    for fname in self._db.listdir(aname):
                        self._db.delete_file(fname)
                    os.removedirs(aname)
            #else:
            #    r = self._modular_symbols.find_one({'_id':id})
            #    #if r.get('complete') is None or r.get('complete')<2:
            #self.check_record(N,k,chi,check_content=True)
            #    self._modular_symbols.update({'_id':id},{"$set":{'complete':int(3)}})
            for N1,k1,chi1,d,prec in self._db.known("N={0} and k={1} and i={2}".format(N,k,chi)):
                #if N < minn or N>maxn or k<mink or k>maxk:
                #    continue
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
                    si = sage_character_to_galois_orbit_number(x)
                if si <> chi:
                    if (N1,k1,chi1) not in problems:
                        problems.append((N1,k1,chi1))
                    print "in file: Chi is wrong! Should be {0}".format(si)
                    if remove == 1:
                        dname = self._db.factor(N1,k1,chi1,d)
                        for fname in self._db.listdir(dname):
                            self._db.delete_file(fname)
                        aname = self._db.ambient(N1,k1,chi1)
                        for fname in self._db.listdir(aname):
                            self._db.delete_file(fname)
                        os.removedirs(aname)
                        if verbose>0:
                            print "removed directory {0}".format(aname)
            if remove == 1 and problems<>[]:
                print "Removed {0} records!".format(len(problems))
        return problems            
    
def precision_needed_for_L(N,k,**kwds):
    r"""
    Returns the precision (number of coefficients) needed to compute the first zero
    of the L-function (as on the LMFDB pages). This bound is taken from there and is probably heuristic.
    """
    pprec = 22 + int(RR(5) * RR(k) * RR(N).sqrt())
    pprec = max(pprec,kwds.get('pprec',100))
    ## Get even hundreds of primes to look nicer.
    return ceil(RR(pprec)/RR(100))*100
