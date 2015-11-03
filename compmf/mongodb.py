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

NOTE: All future interaction with the database should use Conrey's character numbering scheme!

"""


import sage.all
import os,re,sys
from string import join
import pymongo
import gridfs

import multiprocessing
from multiprocessing import Pool

from compmf.filesdb import FilenamesMFDBLoading
from compmf.compute import ComputeMFData
from compmf import character_conversions
from compmf.character_conversions import (
    dirichlet_character_conrey_from_sage_character_number,
    dirichlet_character_conrey_galois_orbit_numbers_from_character_number,
    dirichlet_character_conrey_galois_orbits_reps,
#    dirichlet_character_sage_galois_orbit_rep_from_number,
    dirichlet_character_sage_from_conrey_character_number,
    dirichlet_character_sage_galois_orbits_reps,
    sage_character_to_sage_galois_orbit_number,
    conrey_character_number_from_sage_galois_orbit_number,
    dirichlet_character_conrey_from_sage_galois_orbit_number,
    sage_galois_orbit_number_from_conrey_character_number,
    conrey_character_number_to_conrey_galois_orbit_number,
    dirichlet_group_conrey,
    conrey_character_from_number
)
from sage.all import nth_prime,prime_pi,parallel,loads,dimension_new_cusp_forms,RR,ceil,load,dumps,save,euler_phi,floor,QQ,Integer
from utils import are_compatible
from compmf import clogger

class MongoMF(object):
    def __init__(self,host='localhost',port=37010,db='modularforms2',user='editor',password='',verbose=0,**kwds):
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
        self._db_raw = kwds.get('db_raw','modularforms_raw') # is really a reflection of data which is in the files.
        self._host = host
        self._port = int(port)
        self._verbose = int(verbose)
        self._db_name = db
        self._user = user
        self._password = password
        if pymongo.version_tuple[0] < 3:
            from pymongo import Connection
            _C = Connection(port=port)
            self._mongodb = pymongo.Connection('{0}:{1}'.format(host,port))[db]
            self._mongo_conn = pymongo.Connection('{0}:{1}'.format(host,port))
        else:
            from pymongo.mongo_client import MongoClient
            self._mongodb = MongoClient('{0}:{1}'.format(host,port))[db]
            self._mongo_conn = MongoClient('{0}:{1}'.format(host,port))
        self._mongodb.authenticate(user,password)
        
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
        self._computations = self._mongodb['computations']
        self._galois_orbits = self._mongodb['galois_orbits']
        
        ## The indices we use for the collections given above
        self._collections_indexes = [
        { 'name': 'Modular_symbols.files', 'index':[
            ("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING)],
          'unique':True},
        { 'name': 'Modular_symbols.files', 'index':[
            ("hecke_orbit_label",pymongo.ASCENDING)],
          'unique':True},
        { 'name': 'Newform_factors.files', 'index':[
            ("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
            ("newform",pymongo.ASCENDING)],
          'unique':True},
        { 'name': 'Newform_factors.files', 'index':[
            ("hecke_orbit_label",pymongo.ASCENDING)],
            'unique':True},
        { 'name' : 'ap.files', 'index' : [
            ("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
            ("newform",pymongo.ASCENDING),("prec",pymongo.ASCENDING)],
          'unique':True},
        { 'name' : 'Atkin_Lehner.files', 'index': [
            ("N",pymongo.ASCENDING),("k",pymongo.ASCENDING),("cchi",pymongo.ASCENDING),
            ("newform",pymongo.ASCENDING)],
          'unique':True}
                         ]
        self._sage_version = sage.version.version

    def __repr__(self):
        r"""
        String representation of self.
        """
        s="Modular forms database at mongodb: {0}".format(self._mongo_conn)
        return s

    def create_indices(self,noremove=1,only=None):
        r"""
        Creates indices for our databases
        """
        for r in self._collections_indexes:
            if not only is None:
                if r['name'] <> only:
                    continue
            if r['name'] in self._mongodb.collection_names():
                try:
                    self._mongodb[r['name']].create_index(r['index'],unique=r['unique'])
                except pymongo.errors.DuplicateKeyError as e:
                    print "Error: {0}!".format(e)
#                    print "Removing duplicates!"
#                    keys =  [x[0] for x in r['index']]
#                    self.remove_duplicates(r['name'],keys,dryrun=noremove)
#                    try:
#                        self._mongodb[r['name']].create_index(r['index'],unique=r['unique'])
#                    except pymongo.errors.DuplicateKeyError:
#                        clogger.critical("We could not remove all duplicates!")
            # remove duplicates


#     def remove_duplicates(self,col,keys,dryrun=1):
#         r"""
#         Remove duplicates from mongo db


#         """
#         from sage.all import deepcopy
#         if 'files' not in col:
#             ccol='{0}.files'.format(col)
#         else:
#             ccol = col
#         clogger.debug("ccol={0}".format(ccol))
#         clogger.debug("keys={0}".format(keys))
#         flds = deepcopy(keys); flds.extend(['uploadDate','filename','_id','chi'])

#         clogger.debug("flds={0}".format(flds))
#         nnmax = max(self._mongodb[ccol].find().distinct('N'))
#         args = []
#         if 'ap' in col:
#             step=50
#         else:
#             step = 100
#             flds.append('prec')
#         #h = RR(nnmax)/32.0
#         for j in range(RR(nnmax)/RR(step)):
#             nmin = j*step; nmax = (j+1)*step
#             args.append((col,keys,flds,dryrun,nmin,nmax))
#         return list(self.remove_duplicates32(args))

#     @parallel(ncpus=32)
#     def remove_duplicates32(self,col,keys,flds,dryrun=1,nmin=0,nmax=10000):
#         r"""
#         Remove duplicate records in collection col for which keys are not unique.
#         """

#         #keys = [x[0] for x in keys]
#         fs = gridfs.GridFS(self._mongodb,col.split(".")[0])
#         if 'files' not in col:
#             ccol='{0}.files'.format(col)
#         else:
#             ccol = col
#         s = {}
#         clogger.debug("nmin = {0} \t nmax= {1} \t col={2} \t ccol={3}".format(nmin,nmax,col,ccol))
#         #if ccol=='Newform_factors.files':
#         s = {'N':{"$lt":int(nmax),"$gt":int(nmin)-1}}
#         for r in self._mongodb[ccol].find(s,projection=flds).sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING),('uploadDate',pymongo.DESCENDING)]):
#             id=r['_id']
#             s = {}
#             for k in keys:
#                 try:
#                     s[k]=r[k]
#                 except KeyError as e:
#                     if k=='cchi':
#                         clogger.warning("rec without cchi: r={0}".format(r))
#                         ci = conrey_character_number_from_sage_galois_orbit_number(r['N'],r['chi'])
# #
# #                        c = dirichlet_character_conrey_from_sage_character_number(r['N'],r['chi'])
# #                        ci = c.number()
#                         self._mongodb[ccol].update({'_id':r['_id']},{"$set":{'cchi':ci}})
#                         clogger.debug("Added cchi!")
#                     #raise KeyError,e.message
#             #print "s=",s
#             q = self._mongodb[ccol].find(s,projection=flds).sort('uploadDate',pymongo.ASCENDING)
#             if q.count()==1:
#                 continue
#             for rnew in q:
#                 if rnew['_id']==id:
#                     continue
#                 clogger.debug("s = {0}".format(s))
#                 clogger.debug("Removing record {0} in collection {1}".format(rnew,col))
#                 clogger.debug("Duplicate of {0}".format(r))
#                 if dryrun:
#                     clogger.debug("Not really deleting!")
#                 else:
#                     clogger.debug("We are really deleting {0} from {1}".format(rnew['_id'],fs._GridFS__collection))
#                     fs.delete(rnew['_id'])


#     def remove_duplicates1(self,col,key,dryrun=1):
#         r"""
#         Remove duplicate records in collection col for which keys are not unique.
#         """
#         fs = gridfs.GridFS(self._mongodb,col.split(".")[0])
#         if 'files' not in col:
#             ccol='{0}.files'.format(col)
#         else:
#             ccol = col
#         s = {}
#         for r in self._mongodb[ccol].find().sort('uploadDate',pymongo.DESCENDING):
#             id=r['_id']
#             val = r[key]
#             q = self._mongodb[ccol].find({key:val})
#             if q.count()==1:
#                 continue
#             for rnew in q:
#                 if rnew['_id']==id:
#                     continue
#                 clogger.debug("s = {0}".format(s))
#                 clogger.debug("Removing record {0} in collection {1}".format(rnew,col))
#                 clogger.debug("Duplicate of {0}".format(r))
#                 if dryrun:
#                     clogger.debug("Not really deleting!")
#                 else:
#                     clogger.debug("We are really deleting {0} from {1}".format(rnew['_id'],fs._GridFS__collection))
#                     fs.delete(rnew['_id'])

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


    def show_how_many_complete(self,nrange=[0,50],krange=[2,12]):
        r"""
        Check up to which levels and weights we have a complete set of records...
        """
        from sage.all import flatten,dimension_cusp_forms,Gamma1
        files = self._modular_symbols
        factors = self._newform_factors
        levels = files.distinct('N')
        levels.sort()
        weights = files.distinct('k')
        completed={}
        missing=[]
        for N in range(nrange[0],nrange[1]):
            gal_orbits = character_conversions.dirichlet_group_conrey_galois_orbits_numbers(N)
            ##
            
            #print "Galois orbits for {0}: {1}".format(N,gal_orbits)
            even_orbits = flatten(filter(lambda x:conrey_character_from_number(N,x[0]).is_even(),gal_orbits))
            even_orbits.sort()
            odd_orbits = flatten(filter(lambda x: not conrey_character_from_number(N,x[0]).is_even(),gal_orbits))
            odd_orbits.sort()
            #print "even orbits for {0} = {1}".format(N,even_orbits)
            for k in range(krange[0],krange[1]):
                orbits = files.find({'N':N,'k':k,'complete':{"$gt":int(0)}}).distinct('character_galois_orbit')
                if orbits == []:
                    if dimension_cusp_forms(Gamma1(N),k)>0:
                        print "missing space: N={0}, k={1}".format(N,k)
                    else:
                        continue
                orbits.sort()
                # This returnes all character numbers from all orbits.
                if k % 2 == 0:
                    if orbits <> even_orbits:
                        print "For N={0} and k={1}:".format(N,k)
                        print "orbits in db=",orbits
                        print "even orbits=",even_orbits
                else:
                    if orbits <> odd_orbits:
                        print "For N={0} and k={1}:".format(N,k)                        
                        print "orbits in db=",orbits
                        print "odd orbits=",odd_orbits                        
                    #print N,orbits==even_orbits
                    
                #for r in files.find({'N':N,'complete':{"$gt":int(0)}}):
                
        
    def view_latest(self):
        r"""
        Make a list of the latest additions to the database (all collections except .chunks)
        """

        names = [col for col in self._mongodb.collection_names() if '.chunks' not in col and 'indexes' not in col]
        names.sort()
        longest = max(map(len,names))
        t =  len(self._mongodb[names[0]].find_one()['uploadDate'].ctime())
        print "{0:{width}} \t {1:{dwidth}} \t {2:{dwidth}} \t {3}\n".format("Collection","First","Last","Number of records",width=longest,dwidth=t)

        for col in names:
            no = self._mongodb[col].count()
            if no == 0:
                print "{0:{width}} \t {1:{dwidth}} \t {2:{dwidth}} ".format(col,"Empty","",width=longest,dwidth=t)
                continue
            try:
                r = self._mongodb[col].find().limit(int(1)).sort('uploadDate',int(-1)).next()
                last = r['uploadDate'].ctime()
                r = self._mongodb[col].find().limit(int(1)).sort('uploadDate',int(-1)).next()
                first = r['uploadDate'].ctime()
            except KeyError:
                last="?"
                first = "?"
            print "{0:{width}} \t {1:{dwidth}} \t {2:{dwidth}} \t {3}".format(col,first,last,no,width=longest,dwidth=t)

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

    def get_parameters_from_input(self,*args,**kwds):
        r"""
        Extract the parameters: Level,weight,conrey_character_no
        from the input.
        """
        if isinstance(args[0],basestring):
            N,k,ci = param_from_label(args[0])
        elif len(args)==3:
            N,k,ci = args
        elif len(args)==4:
            N,k,ci,d = args
        #    
        return int(N),int(k),int(ci)
    
    def get_ambient(self,N,k,ci,**kwds):
        r"""
        Return an ambient space M(N,k,ci') for some i' in the Galois orbit of i

        INPUT:
         - 'N' -- integer (level) or string of the form 'N.k.i' (label)
         - 'k' -- integer (weight)
         - 'i' -- integer (Character number modulo N in Conrey's numbering)
        keywords:
            - get_record (False) -- if True return the database record instead of the space.
        """
        ambient_id = kwds.get('ambient_id',None)
        if isinstance(N,basestring):
            N,k,ci = param_from_label(N)
        N = int(N); k = int(k); ci = int(ci)
        if ambient_id is None:
            s = {'N':N,'k':k,'cchi':ci}
            f = self._modular_symbols.find_one(s)
            if f is None:
                s = {'N':N,'k':k,'character_galois_orbit':{"$in":[ci]}}
                f = self._modular_symbols.find_one(s)
            if f is None:
                return None
            else:
                if kwds.get('get_record',False):
                    return f
                ambient_id = f['_id']
        return self.load_from_mongo('Modular_symbols',ambient_id)

    def get_dimc(self,N,k,ci):
        r"""
        Get dimension of cusp forms in S_k(N,i).

        """
        if isinstance(N,basestring):
            N,k,ci = param_from_label(N)
        r = self._modular_symbols.find_one({'N':int(N),'k':int(k),'cchi':int(ci)},projection=['dimc'])
        if r is None:
            return -1
        return r.get('dimc',-1)

    def number_of_factors(self,N,k,ci,**kwds):
        r"""
        Return the number of factors.
        """
        if isinstance(N,basestring):
            N,k,ci = param_from_label(N)
        s = {'N':N,'k':k,'cchi':int(ci)}
        nf = self._newform_factors.find(s).count()
        if nf > 0:
            return nf
        # make sure we count the number of factors for a given element in the Galois orbit
        s= {'N':N,'k':k,'character_galois_orbit':{"$in":[int(ci)]}}
        q = self._newform_factors.find_one(s)
        if q is None:
            return 0
        s= {'N':N,'k':k,'cchi':q['cchi']}
        return self._newform_factors.find(s).count()

    def get_factors(self,N,k,ci,d=None,sources=['mongo','files']):
        r"""
        Get factor nr. d of the space M(N,k,i)
        ci should be in the conrey naming scheme. Any conversions should take place before calling this function. 
        """
        from utils import orbit_index_from_label,param_from_label
        if isinstance(N,basestring):
            N,k,ci,d = param_from_label(N)
        if isinstance(d,basestring):
            d = orbit_index_from_label(d)
        res = None
        if sources == ['mongo']:
            ## Fetch the record in the database corresponding to cchi=i or the Galois orbit of
            ## chi_N,i  (in Conrey's character naming scheme)
            s = {'N':int(N),'k':int(k),'cchi':int(ci)}            
            if self._newform_factors.find(s).count()==0:
                s = {'N':int(N),'k':int(k),'character_galois_orbit':{"$all":[int(ci)]}}
            if not d is None:
                s['newform']=int(d)
            for r in self._newform_factors.find(s):
                if res is None:
                    res = {}
                newform = r['newform']
                fid = r['_id']
                f = self.load_from_mongo('Newform_factors',fid)
                t = (int(N),int(k),int(ci),int(newform))
                if not d is None:
                    res = f
                else:
                    res[t] = f
            return res
        elif sources == ['files']:
            # The files are named according to Galois orbits.
            on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)
            if not d is None:
                res = self._db.load_factor(N,k,on,d)
            else:
                if res is None:
                    res = {}
                nfacts = self._db.number_of_known_factors(N,k,on)
                for d in range(nfacts):
                    F = self._db.load_factor(N,k,on,d)
                    if not F is None:
                        res[d]=F
            return res
        elif len(sources)>1:
            for ss in sources:
                res = self.get_factors(N,k,ci,d,sources=[ss])
                if not res is None:
                    return res
        else:
            clogger.critical("Can not load factors, unknown source:{0}".format(sources))
        return res

    def get_aps(self,N,k,ci,d=None,prec_needed=0,coeffs=False,sources=['mongo','mongo_raw','file']):
        r"""
        Get the lists of Fourier coefficients for the space M(n,k,i) and orbit nr. d
        i is in the Conrey naming scheme.

        INPUT:

        - sources -- describe a list of sources to check, including:
            * 'mongo' -- high level sage objects in mongodb
            * 'mongo_raw' -- low-level objects in mongodb (will need some construction of objects)
            * 'files' -- low-level objects from files.


        """

        clogger.debug("find aps with N,k,ci,d={0} from sources:{1}".format((N,k,ci,d),sources))
        res = None
        if sources == ['mongo']:
            s = {'N':int(N),'k':int(k),'cchi':int(ci)}
            if self._aps.find(s).count()==0:
                ## Fetch the record in the database corresponding to the Galois orbit of
                ## chi_N,i  (in Conrey's character naming scheme)
                s = {'N':int(N),'k':int(k),'character_galois_orbit':{"$all":[int(ci)]}}
            if not d is None:
                s['newform']=int(d)
            clogger.debug("find aps with s={0} from sources:{1}".format(s,sources))

            for r in self._aps.find(s):
                if res is None:
                    res = {}
                fid=r['_id']; newform = r['newform']; prec = r['prec']
                meta = {'cputime':r.get('cputime'),'version':r.get('sage_version')}
                try:
                    E,v = self.load_from_mongo('ap',fid)
                except ValueError as e:
                    clogger.debug("Can not load ap's: {0}".format(e.message))
                    clogger.debug("Removing these ap's from database!: match={0}".format(s))
                    self._aps.delete_one({'_id':fid})
                    return None
                #clogger.debug("id={0} and E={1}".format(fid,E))
                t = (int(N),int(k),int(ci),int(newform))
                if not res.has_key(newform):
                    res[newform]={}
                if prec_needed == 0 or coeffs == False:
                    res[newform][prec]=(E,v,meta)
                else:
                    if prec >= prec_needed and coeffs:
                        res[newform][prec] = E*v
                if not d is None:
                    res = res[newform]
        elif sources == ['files']:
            # The files are named according to Galois orbits but here we simply request
            # the character (even though we might get another character in the same orbit)
            #on = character_conversions.conrey_character_number_to_conrey_galois_orbit_number(N,ci)[1]
            try:
                res = self._db.load_aps(N,k,ci,d,numc=prec_needed)
            except Exception as e:
                clogger.critical("Could not get ap's from file:{0}".format(e))
                res = None
        elif len(sources)>1:
            for ss in sources:
                res = self.get_aps(N,k,ci,d,prec_needed,coeffs,[ss])
                #print "res=",res
                if not res is None:
                    return res
        return res


    ### We want to see which computations are going on currently.
    def register_computation(self,level,weight,cchi,typec='mf'):
        r"""
        Insert a record in the database registring the start of a computation.
        """
        import datetime
        import time
        import os
        r = {'startTime':datetime.datetime.fromtimestamp(time.time()),
             'server':os.uname()[1],
             'pid':os.getpid(),
             'type':typec,
             'N':int(level), 'k':int(weight),'cchi':int(cchi)}
        if self._computations.find({'N':r['N'],'k':r['k'],'cchi':int(cchi)}).count()>0:
            return None
        fid = self._computations.insert(r)
        return fid 

    def register_computation_closed(self,cid):
        import datetime
        import time
        now = datetime.datetime.fromtimestamp(time.time())
        self._computations.update({"_id":cid},{"$set":{"stopTime":now}})

    def find_running_computations(self,typec='mf'):
        import datetime
        import time
        now = datetime.datetime.fromtimestamp(time.time())
        res = []
        for t in ['mf','wmf']:
            if t == 'mf':
                print "Modular forms computations"
            else:
                print "WebModularForms/NewForms computations"
            for r in self._computations.find({'stopTime':{"$exists":False},'type':t}):
                duration = str(now - r['startTime']).split(".")[0]
                print "{0},{1},{2} \t\t {3} \t\t {4} \t {5}".format(r['N'],r['k'],r['chi'],r['startTime'],duration,r['pid'])
            
    def clear_running_computations(self,typec='mf'):
        res = self._computations.delete_many({"type":typec})
        print "Removed {0} computations from db!".format(res.deleted_count)

        
                    



# def unwrap_compute_space(D,*args,**kwds):
#     r"""
#     To overcome some unpickling problems with the builtin parallel decorators.
#     """
#     args,kwds = args
#     clogger.debug("in unwrap: args={0} kwds={1}".format(args,kwds))
#     r = CompMF.compute_and_insert_one_space(*args,**kwds)
#     return r

    
class CompMF(MongoMF):
    r"""
    Class for computing modular forms and inserting records in a mongodb as well as a file system database.
    """

    def __init__(self,datadir='',host='localhost',port=37010,db='modularforms2',user='',password='',verbose=0,**kwds):
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
        super(CompMF,self).__init__(host,port,db,user=user,password=password,verbose=verbose,**kwds)
        self._datadir = datadir
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

    def __reduce__(self):
        return (type(self),(self._datadir,self._host,self._port,self._db_name,self._user,self._password,self._verbose))
    

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

    # def convert_to_mongo_all(self,par=0,**kwds):
    #     r"""
    #     Convert all records in the files to MongoDB
    #     """
    #     nrange = self._db.known_hosts_levels()
    #     nrange.sort()
    #     clogger.debug("nrange={0}".format(nrange))
    #     return self.convert_to_mongo(nrange)

    # def convert_to_mongo(self,N=None,k=None,ncpus=1,**kwds):
    #     r"""
    #     Converting records which exists (possibly partly) in the files for a given level and insert into the mongodb.
    #     """
    #     if isinstance(N,basestring):
    #         s = N
    #     else:
    #         s = ""
    #         if s == "":
    #             if N<>None:
    #                 s = "N={0}".format(N)
    #             if k<>None:
    #                 s+= " k={0}".format(k)
    #     args = []
    #     clogger.debug("Converting records matching pattern {0}".format(s))
    #     for (N,k,i,newforms,nap) in self._db.known(s):
    #         if not are_compatible(N,k,i):
    #             continue
    #         print "Check ",N,k,i
    #         if kwds.get('force',False):
    #             q = self._newform_factors.find({'N':int(N),'k':int(k),'cchi':int(i)})
    #             if q.count() == 0:
    #                 args.append((N,k,i))
    #         else:
    #             args.append((N,k,i))
    #     self.get_or_compute_spaces(args,**kwds)

    def get_or_compute_spaces(self,args,**kwds):
        r"""
        Get or compute a list of spaces in parallel.

        """
        ncpus = kwds.pop('ncpus',1)
        pool = Pool(processes=ncpus)
        clogger.debug("ncpus={0}".format(ncpus))        
        n = len(args)
        if n > 100:
            chunksize = 10
        else:
            chunksize = 1
        args = [(x,kwds) for x in args]
        clogger.debug("args={0}".format(args))
        results = pool.imap_unordered(self.unwrap,args,chunksize)
#        results = pool.imap_unordered(unwrap_compute_space,args,chunksize)
#        results = pool.map_async(self.unwrap,args)        
        return results
        for res in results:
            res.get()
        pool.close()
        pool.join()
#        if ncpus==8:
#            return list(self._compute_and_insert_one_space8(args,**kwds))
#        elif ncpus==16:
#            return list(self._compute_and_insert_one_space16(args,**kwds))
#        elif ncpus==32:
#            return list(self._compute_and_insert_one_space32(args,**kwds))
#        else:
#            return [self.compute_and_insert_one_space(x[0],x[1],x[2],**kwds) for x in args]

    ## Different levels of parallelization
    def unwrap(self,args):
        args,kwds = args
        clogger.debug("in self.unwrap: args={0} kwds={1}".format(args,kwds))
        return self.compute_and_insert_one_space(args[0],args[1],args[2],**kwds)
        clogger.debug("unwrapped")

# @parallel(ncpus=8)
# #    @classmethod
# #    @staticmethod
#     def _compute_and_insert_one_space8(self,N,k,i,**kwds):
#         return self.compute_and_insert_one_space(N,k,i,**kwds)

#     @parallel(ncpus=16)
# #    @staticmethod
#     def _compute_and_insert_one_space16(self,N,k,i,**kwds):
#         return self.compute_and_insert_one_space(N,k,i,**kwds)

#     @parallel(ncpus=32)
# #    @staticmethod
#     def _compute_and_insert_one_space32(self,N,k,i,**kwds):
#         return self.compute_and_insert_one_space(N,k,i,**kwds)

    def compute_and_insert_one_space(self,N,k,ci,**kwds):
        r"""
        If data for the given space M(N,k,i) exists in either mongo or files database we fetch this data (and possible complete if e.g. more coefficients are needed) and then insert the result into both databases unless explicitly told not to.

        ci is a Conrey character number.
        
        """
        N = int(N); k = int(k); ci = int(ci)
        clogger.debug("In Compute and/or Insert N,k,i= {0} kwds={1}".format((N,k,ci),kwds))        
        sage.modular.modsym.modsym.ModularSymbols_clear_cache()
        #
        #ci = conrey_character_number_from_sage_galois_orbit_number(N,i)
        #c = dirichlet_character_conrey_from_sage_character_number(N,i)
        #ci = c.number()        
        cid = self.register_computation(level=N,weight=k,cchi=ci,typec='mf')
        if cid == None:
            clogger.critical("Computation with N,k,i={0} is already underway!!".format((N,k,ci)))
        clogger.debug("Compute and/or Insert {0}".format((N,k,ci)))
        clogger.debug("Computing ambient modular symbols")
        if kwds.get('Nmax',0)<>0 and kwds.get('Nmax')>N:
            return
        ambient_fid = self.compute_ambient(N,k,ci,**kwds)
        # if the space is empty we return None and exit;
        if ambient_fid is None: ## If something went wrong or there isn't even Eisenstein series.
            clogger.debug("Ambient is None!")
            return True
        if isinstance(ambient_fid,int) and ambient_fid ==1: ## If we have no cusp forms
            clogger.debug("Ambient space of cusp forms is 0 dimensional!")
            self._db.update_known_db((N,k,ci,0,0))
            return True
        # Insert ambient modular symbols
        kwds['ambient_id']=ambient_fid
        clogger.debug("Getting factors!")
        factor_ids = self.compute_factors(N,k,ci,**kwds)
        kwds['factor_ids']=factor_ids
        # Necessary number of coefficients for L-function computations.
        if not factor_ids is None and factor_ids <> []:
            pprec = precision_needed_for_L(N,k,pprec=100)
            clogger.debug("Computing aps to prec {0} for ci = {1}".format(pprec,ci))
            ap = self.compute_aps(N,k,ci,pprec,**kwds)
            clogger.debug("Got {0}".format(ap))
            try:
                atkin_lehner = self.compute_atkin_lehner(N,k,ci,**kwds)
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
        sage.modular.modsym.modsym.ModularSymbols_clear_cache()
        self.register_computation_closed(cid)
        return True

    def compute_atkin_lehner(self,N,k,ci,**kwds):
        r"""
        Compute the Atkin-Lehner eigenvalues of space (N,k,i).
        Only for trivial character!
        """
        verbose = kwds.get('verbose',0)
        from character_conversions import (
            dirichlet_character_conrey_galois_orbit_numbers_from_character_number,
            conrey_character_number_to_conrey_galois_orbit_number,
            sage_galois_orbit_number_from_conrey_character_number)
        #c = dirichlet_character_conrey_from_sage_galois_orbit_number(N,i)
        #ci = c.number()
        #ci = conrey_character_number_from_sage_galois_orbit_number(N,i)
        
        if not ci == 1: #c.is_trivial(): #or c.multiplicative_order()==2):
            return []
        al_in_mongo = self._atkin_lehner.find({'N':int(N),'k':int(k),'cchi':int(ci)}).distinct('_id')
        fs = gridfs.GridFS(self._mongodb, 'Atkin_Lehner')
        orbit =  dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
        on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)
        if len(al_in_mongo)==0:
            ambient = self.get_ambient(N,k,ci,**kwds)
            number_of_factors = self.number_of_factors(N,k,ci)
            self._computedb.compute_atkin_lehner(N,k,ci,M=ambient,m=number_of_factors,verbose=verbose)
            m = self._computedb._db.number_of_known_factors(N,k,on)
            for d in range(m):

                atkin_lehner_file = self._computedb._db.factor_atkin_lehner(N,k,ci,d,True)
                #atkin_lehner = load(atkin_lehner_file)
                try:
                    atkin_lehner = open(atkin_lehner_file,'r').read()
                except IOError:
                    raise ArithmeticError,"Error with opening the Atkin-Lehner file {0}! chi order={1}".format(atkin_lehner_file,c.multiplicative_order())
                try:
                    meta = load(self._db.meta(atkin_lehner_file))
                except IOError:
                    meta = {}
                fname = 'atkin_lehner_evs-{0:0>5}-{1:0>3}-{2:0>3}-{3:0>3}'.format(N,k,ci,d)
                sage_i = sage_galois_orbit_number_from_conrey_character_number(N,ci)
                fid = fs.put(dumps(atkin_lehner),filename=fname,
                             N=int(N),k=int(k),chi=int(sage_i),newform=int(d),cchi=int(ci),
                             character_galois_orbit = orbit,
                             cputime = meta.get("cputime",""),
                             sage_version = meta.get("version",""))
                al_in_mongo.append(fid)
        return al_in_mongo


    def compute_ambient(self,N,k,ci,**kwds):
        r"""
        Compute the ambient space with the given parameters and insert it in mongo if it is not there.
        """
        if not are_compatible(N,k,ci):
            return None
        files_ms = self._modular_symbols
        fs_ms = gridfs.GridFS(self._mongodb, 'Modular_symbols')
        verbose = kwds.get('verbose',0)
        #ci = conrey_character_number_from_sage_galois_orbit_number(N,i)
#        ci = dirichlet_character_conrey_from_sage_character_number(N,i).number()
        save_in_file = kwds.get('save_in_file',True)
        compute = kwds.get('compute',self._do_computations)
        # We first see if this space is already in the mongo database.
        rec = files_ms.find_one({'N':int(N),'k':int(k),'cchi':int(ci)})
        ambient_in_mongo = 0
        ambient_in_file = 0
        ambient = None
        cputime = None
        clogger.debug("In compute ambient: {0}".format((N,k,ci)))
        if rec<>None:
            ambient_in_mongo = rec['_id']
            # Check factors
            clogger.debug("Have ambient space!")
            clogger.debug("ambient_in_mongo={0}".format(rec))
        # Next see if we have it in the files database. If not we will compute it.
        try:
            ambient = self._db.load_ambient_space(N,k,ci)
            clogger.debug("Loaded ambient={0}".format(ambient))
            ambient_in_file = 1
        except ValueError:
            clogger.debug("Could not load ambient!")
            ambient = None
            pass
        if ambient_in_mongo <> 0 and ambient_in_file==0:
            ambient = loads(fs_ms.get(ambient_in_mongo).read())
            cputime = rec.get('cputime')

            clogger.debug("Space {0},{1},{2} is in mongo but not in file!".format(N,k,ci))
            clogger.debug("ambient={0}".format(ambient))
        if ambient_in_file == 0:
            if not compute:
                clogger.debug("Space {0},{1},{2} not computed at all!".format(N,k,ci))
                return None
            clogger.debug("Compute and/or save to file!")
            self._computedb.compute_ambient_space(N,k,ci,M=ambient,tm=cputime)
            ambient_in_file = 1
        if ambient_in_mongo == 0 and ambient_in_file == 1: #not ambient is None:
            metaname = self._db.space(N,k,ci,False)+"/ambient-meta.sobj"
            clogger.debug("metaname={0}".format(metaname))
            try:
                meta = load(metaname)
            except (OSError,IOError):
                meta={}                
            ## Note that we have to update the number of orbits.
            if ambient is None:
                ambient = self._db.load_ambient_space(N,k,ci)
            dima = int(ambient.dimension())
            dimc = int(ambient.cuspidal_submodule().dimension())
            dimn = int(ambient.cuspidal_submodule().new_submodule().dimension())
            clogger.debug("Save ambient to mongodb! ambient={0}:{1}".format((N,k,ci),ambient))
            on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)
            orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
            sage_orbit_no = sage_galois_orbit_number_from_conrey_character_number(N,ci) 
            fname = "gamma0-ambient-modsym-{0}".format(self._db.space_name(N,k,ci).split("/")[1])
            fid = None
            try:
                dump_ambient = dumps(ambient)
            except Exception as e:
                clogger.debug("Could not dump the ambient space with {0}! : {1}".format((N,k,i),e))
            try:
                fid = fs_ms.put(dumps(ambient),filename=fname,
                                N=int(N),k=int(k),nfactors=int(0),
                                sage_orbit_no=int(sage_orbit_no[1]),
                                orbits=int(0),
                                space_label="{0}.{1}.{2}".format(N,k,ci),
                                space_orbit_label="{0}.{1}.{2}".format(N,k,on),
                                conrey_galois_orbit_number=int(on[1]),
                                dima=dima,dimc=dimc,dimn=dimn,
                                character_galois_orbit=orbit,
                                cchi=int(ci),
                                cputime = meta.get("cputime",""),
                                sage_version = meta.get("version",""))
            except gridfs.errors.FileExists as e:
                clogger.debug("We can not insert the same record twice! Error:{0}".format(e.message))
                rec = files_ms.find_one({'N':int(N),'k':int(k),'chi':int(i)})
                if rec is None:
                    clogger.critical("We could nt find the double record!")
                else:
                    fid = rec['_id']
        else:
            fid = ambient_in_mongo
        clogger.debug("fid={0}".format(fid))
        if fid == 0 or ambient is None:
            fid = None
        return fid

    def compute_factors(self,N,k,ci,**kwds):
        r"""
        Compute / store newform factors of parameters N,k,i
        """
        from utils import orbit_label
        verbose = kwds.get('verbose',0)
        ambient_id = kwds.get('ambient_id',None)
        compute = kwds.get('compute',self._do_computations)
        if ambient_id is None:
            ambient_id = self.compute_ambient(N,k,ci,**kwds)
        on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)

        #ci = conrey_character_number_from_sage_galois_orbit_number(N,i)
        if ambient_id is None:
            clogger.debug("No ambient space!")
            ambient_id = self.compute_ambient(N,k,ci,**kwds)
        files_fact = self._newform_factors
        files_ms = self._modular_symbols
        fs_fact = gridfs.GridFS(self._mongodb, 'Newform_factors')
        factors_in_file = self._db.number_of_known_factors(N,k,ci)
        factors_in_mongo = files_fact.find({'ambient_id':ambient_id}).distinct('_id')
        dimc = self.get_dimc(N,k,ci)
        if dimc == 0:
            return []  # There are no factors to compute
        clogger.debug("factors: in_file={0} in mongo={1}".format(factors_in_file,factors_in_mongo))
        if factors_in_file == 0 and  len(factors_in_mongo)==0:
            if not compute:
                clogger.debug("No factors exist in db and we do not compute anything!")
                return None
                # Else compute and insert into the files.
            factors_in_file = self._computedb.compute_decompositions(N,k,ci)
            clogger.debug("Computing factors! m={0}".format(factors_in_file))

        if len(factors_in_mongo)==0:
            fname = "gamma0-factors-{0}".format(self._db.space_name(N,k,ci))
            clogger.debug("Inserting factors into mongo! fname={0}".format(fname))
            num_factors_in_file = self._db.number_of_known_factors(N,k,ci)
            ambient = self.get_ambient(N,k,ci,ambient_id=ambient_id)
            orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)

            for d in range(num_factors_in_file):
                try:
                    factor = self._db.load_factor(N,k,ci,d,M=ambient)
                except RuntimeError:
                    ## We probably need to recompute the factors
                    clogger.debug("The factors from file was corrupt / empty. Need to compute them anew!")

                    factors_in_file = self._computedb.compute_decompositions(N,k,ci)
                    try:
                        factor = self._db.load_factor(N,k,ci,d,M=ambient)
                    except RuntimeError:
                        raise ArithmeticError,"Could not get factors for {0}".format((N,k,ci))
                metaname = self._db.space(N,k,ci,False)+"/decomp-meta.sobj"
                clogger.debug("metaname={0}".format(metaname))
                try:
                    meta = load(metaname)
                except (OSError,IOError):
                    meta={}
                if factor==None:
                    clogger.debug("Factor {0},{1},{2},{3} not computed!".format(N,k,ci,d))
                    continue
                fname1 = "{0}-{1:0>3}".format(fname,d)
                label = orbit_label(d)
                sage_i = sage_galois_orbit_number_from_conrey_character_number(N,ci)      
                clogger.debug("{filename},{N},{k},{chi},{cchi},{character_galois_orbit},{conrey_galois_orbit_number},{newform},{cputime},{sage_version},{ambient_id},{hecke_orbit_label},{v}".format(
                    filename=fname1,N=int(N),k=int(k),chi=int(sage_i[1]),
                    cchi=int(ci),
                    character_galois_orbit=orbit,
                    conrey_galois_orbit_number=int(on[1]),
                    newform=int(d),
                    cputime = meta.get("cputime",""),
                    sage_version = meta.get("version",""),
                    ambient_id=ambient_id,
                    hecke_orbit_label='{0}.{1}.{2}{3}'.format(N,k,ci,label),
                    v=int(1)))

                facid = fs_fact.put(dumps(factor),filename=fname1,
                                    N=int(N),k=int(k),chi=int(sage_i[1]),
                                    cchi=int(ci),
                                    character_galois_orbit=orbit,
                                    conrey_galois_orbit_number=int(on[1]),
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
                self._computedb.compute_decompositions(N,k,ci,D=D)

        return factors_in_mongo




    def compute_aps(self,N,k,ci,pprec=None,**kwds):
        r"""
        Compute & store aps
        """
        from utils import orbit_label
        if pprec is None:
            pprec = precision_needed_for_L(N,k)
        pprec = int(pprec)
        ambient_id = kwds.get('ambient_id',None)
        if ambient_id is None:
            ambient_id = self.compute_ambient(N,k,ci,**kwds)
            kwds['ambient_id']=ambient_id
        ambient = self.get_ambient(N,k,ci,ambient_id=ambient_id,verbose=0)
        num_factors = len(kwds.get('factors_ids',self.compute_factors(N,k,ci,**kwds)))

        compute = kwds.get('compute',self._do_computations)
        verbose = kwds.get('verbose')
        #        c = dirichlet_character_conrey_from_sage_character_number(N,i)
        #        ci = c.number()
        #ci = conrey_character_number_from_sage_galois_orbit_number(N,i)
        orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
        fs_ap = gridfs.GridFS(self._mongodb, 'ap')
        fs_v = gridfs.GridFS(self._mongodb, 'vector_on_basis')
        key = {'N':int(N),'k':int(k),'cchi':int(ci),'prec' : {"$gt": int(pprec -1) }}
        clogger.debug("key={0}".format(key))
        aps_in_mongo = self._aps.find(key).distinct('_id')
        aps_in_file = 0
        clogger.debug("Already have {0} ap lists in mongodb! Need at least {1}".format(len(aps_in_mongo),num_factors))
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
                N,k,ci,d = key
                on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)
                sage_i = sage_galois_orbit_number_from_conrey_character_number(N,ci)      
                E,v,meta = value
                if isinstance(E,tuple):
                    E,v = E
                clogger.debug("E={0}".format(E))
                clogger.debug("v=vector of length {0}".format(len(v)))
                clogger.debug("meta={0}".format(meta))
                fname = "gamma0-aplists-{0}".format(self._db.space_name(N,k,ci).split("/")[-1])
                fname1 = "{0}-{1:0>3}-{2:0>5}".format(fname,d,pprec)
                label = orbit_label(d)
                clogger.debug("label={0}".format((N,k,ci,d)))
                # delete if exists
                r = self._aps.find_one({'filename':fname1})
                if not r is None:
                    fs_ap.delete(r['_id'])
                try:
                    # check again if we have this record in the gridfs db
                    clogger.debug("ambient id: {0} prec={1}".format(ambient_id,pprec))             
                    apid = fs_ap.put(dumps( (E,v)),filename=fname1,
                                     N=int(N),k=int(k),chi=int(sage_i[1]),cchi=int(ci),
                                     character_galois_orbit=orbit,
                                     conrey_galois_orbit_number=int(on[1]),
                                     newform=int(d),
                                     hecke_orbit_label='{0}.{1}.{2}{3}'.format(N,k,ci,label),
#                                     prec = int(pprec))
                                     cputime = meta.get("cputime",""),
                                     sage_version = meta.get("version",""),
                                     ambient_id=ambient_id,
                                     prec = int(pprec))
                    aps_in_mongo.append(apid)
                    clogger.debug("We could insert {0} fname={1}".format(apid,fname1))
                except ValueError as e: #gridfs.errors.FileExists as e:
                    clogger.critical("Could not insert coefficients for fname={0}: Error:{1}".format(fname1,e.message))
                    q = self._aps.find({'hecke_orbit_label':'{0}.{1}.{2}{3}'.format(N,k,ci,label)})
                    clogger.critical("We have {0} records in the database!".format(q.count()))
                clogger.debug("inserted aps :{0} ".format((num_factors,apid)))
                # Also insert the corresponding v separately (to protect against changes in sage)
                fnamev = "gamma0-ambient-v-{0}-{1:0>3}".format(self._db.space_name(N,k,ci),d)
                clogger.debug("fnamev={0}".format(fnamev))
                vid = fs_v.put(dumps(v),filename=fnamev,
                               newform=int(d),
                               character_galois_orbit=orbit,
                               N=int(N),k=int(k),chi=int(sage_i[1]),cchi=int(ci),
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
            for d,val in aps.iteritems():
                clogger.debug("d={0}".format(d))
                clogger.debug("type(val)={0}".format(type(val)))
                if isinstance(val,dict):
                    clogger.debug("val.keys()={0}".format(val.keys()))
                #N,k,i,d = key
                #d = key
                for prec in val.keys():
                    E,v,meta = val[prec]
                    aplist_file = self._db.factor_aplist(N, k, ci, d, False, prec)
                    apdir = join(aplist_file.split("/")[0:-1],"/")
                    if not self._db.isdir(apdir):
                        self._db.makedirs(apdir)
                    clogger.debug("aplist_file={0}, meta = {1}".format(aplist_file,meta))
                    save(E, aplist_file)

                    vname = self._db.factor_dual_eigenvector(N, k, ci, d)
                    if not self._db.path_exists(vname):
                        save(v,vname)
                    save(meta, self._db.meta(aplist_file))
                    res.append(vname)
            return res
        # See if we have the coefficients in a file or not
        q = self._db.known("N={0} and k={1} and i={2}".format(N,k,ci))
        for r in q:
            if r[4]>=pprec:
                aps_in_file = 1
                break
        clogger.debug("Have ap lists in filesdb : {0}".format(aps_in_file))

        # If we have coefficients both in mongo and files we don't do anything.
        if len(aps_in_mongo) == num_factors and aps_in_file==1:
            return aps_in_mongo
        elif len(aps_in_mongo) < num_factors:

            if aps_in_file==0:
                if not compute:
                    return []
                ## No coefficients in either mongo or files => we compute and save if desired
                clogger.debug("Computing aplist! m={0} with pprec={1}".format(num_factors,pprec))
                aps = self._computedb.compute_aplists(N,k,ci,0,pprec,ambient=ambient,save=self._save_to_file)
                if aps == 0:
                    aps = {}
                    for d in range(num_factors):
                        E,v,meta  = self._db.load_aps(N,k,ci,d,ambient=ambient,numc=pprec)
                        aps[(N,k,ci,d)] = E,v,meta
            else:
                aps = {}
                for d in range(num_factors):
                    E,v,meta = self._db.load_aps(N,k,ci,d,ambient=ambient,numc=pprec)
                    aps[(N,k,ci,d)] = E,v,meta

            if not isinstance(aps,dict):
                clogger.critical("APS = {0}".format(aps))
            if aps == None or (isinstance(aps,dict) and len(aps.values()[0])<>3) or aps==-1:
                clogger.critical("APS: {0},{1},{2},{3} could not be computed!".format(N,k,ci,d))
                return aps
            return insert_aps_into_mongodb(aps)
        elif len(aps_in_mongo) >= num_factors and len(aps_in_mongo)>aps_in_file and self._save_to_file:
            ### We have coefficients in mongo, need to save them to file
            aps = self.get_aps(N,k,ci,sources=['mongo'])
            clogger.debug("Need to insert aps into the files! num_Factors={0}".format(num_factors))
            return insert_aps_into_filesdb(aps)
        elif len(aps_in_mongo)>=num_factors: #  we are ok anyway
            clogger.debug("We have anough coefficients in mongo and do not want to write to file!")
        else:
            clogger.critical("aps for: {0},{1},{2} could not be computed!".format(N,k,ci))
        return aps_in_mongo


    






   


    # def check_character(self,N,k,chi,remove=1,files_separately=0):
    #     #if N % 10 == 0:
    #     clogger.debug("Checking N={0}, k={1}, chi={2}".format(N,k,chi))

    #     problems=[]
    #     i = 0
    #     for r in self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(chi)}):
    #         if i>1:
    #             clogger.warning("Multiple records for {0}!".format((N,k,chi)))
    #         i=i+1
    #         sage.modular.modsym.modsym.ModularSymbols_clear_cache()
    #         ms = gridfs.GridFS(self._mongodb,'Modular_symbols')
    #         aps = gridfs.GridFS(self._mongodb,'aps')
    #         factors = gridfs.GridFS(self._mongodb,'Newform_factors')
    #         al = gridfs.GridFS(self._mongodb,'Atkin_Lehner')
    #         id = r['_id']
    #         cchi = r.get('cchi')
    #         clogger.debug("checking record: {0}".format(r))
    #         if cchi is None:
    #             clogger.debug("We don't have a Conrey character for r={0}".format(r))
    #             #c = dirichlet_character_conrey_from_sage_character_number(N,chi)
    #             #cchi = c.number()
    #             cchi =conrey_character_number_from_sage_galois_orbit_number(N,chi)
    #             for col in self._file_collections:
    #                 if col=='Modular_symbols':
    #                     self._mongodb['{0}.files'.format(col)].update({'_id':id},{"$set":{'cchi':cchi}})
    #                 else:
    #                     self._mongodb['{0}.files'.format(col)].update({'ambient_id':id},{"$set":{'cchi':cchi}})
    #         #clogger.debug("Get Modular symbols from Mongo! col={0} id={1}".format(self._modular_symbols_collection,id))
    #         M = self.load_from_mongo(self._modular_symbols_collection,id)
    #         #clogger.debug("Got Modular symbols from Mongo!")
    #         x = M.character()
    #         if N == 1:
    #             si = 0
    #             ci = 1
    #         else:
    #             si = sage_character_to_sage_galois_orbit_number(x)
    #             ci = conrey_character_number_from_sage_galois_orbit_number(N,si)
    #         if si <> chi or ci<>cchi:
    #             problems.append((N,k,chi))
    #             clogger.debug("Chi is wrong! Should be {0} not {1}".format(si,chi))
    #             clogger.debug("CChi is wrong! Should be {0} not {1}".format(ci,cchi))
    #             clogger.debug("r={0}".format(r))
    #             raise ArithmeticError()
    #             if remove == 1:
    #                 # First delete from mongo
    #                 ms.delete(id)
    #                 # delete aps
    #                 for rf in  self._aps.find({'_ambient_id':id},projection=['_id']):
    #                     aps.delete(rf['_id'])
    #                 # delete factors
    #                 for rf in  self._newform_factors.find({'_ambient_id':id},projection=['_id','newform']):
    #                     fid = rf['_id']
    #                     factors.delete(fid) # Delete factor
    #                 #    dname = self._db.factor(N,k,chi,rf['newform'])
    #                 #    for fname in self._db.listdir(dname):
    #                 #        self._db.delete_file(fname)
    #                 #aname = self._db.ambient(N,k,chi)
    #                 #dname = join(s.split("/")[0:-1],"/")
    #                 #for fname in self._db.listdir(dname):
    #                 #    self._db.delete_file(fname)
    #                 #os.removedirs(aname)
    #         else:
    #             clogger.debug("Characters are correct!")
    #         #    r = self._modular_symbols.find_one({'_id':id})
    #         #    #if r.get('complete') is None or r.get('complete')<2:
    #         #self.check_record(N,k,chi,check_content=True)
    #         #    self._modular_symbols.update({'_id':id},{"$set":{'complete':int(3)}})
    #         s = "N={0} and k={1} and i={2}".format(N,k,chi)
    #         clogger.debug("searching files for: {0}".format(s))

    #         for N1,k1,chi1,d,prec in self._db.known(s):
    #             #if N < minn or N>maxn or k<mink or k>maxk:
    #             #    continue
    #             clogger.debug("  in file with : {0}".format((N1,k1,chi1,d,prec)))
    #             sage.modular.modsym.modsym.ModularSymbols_clear_cache()
    #             if not files_separately==1 or (N1,k1,chi1) in problems:
    #                 continue
    #             if not self._db.path_exists(self._db.ambient(N1,k1,chi1)):
    #                 continue
    #             M = self._db.load_ambient_space(N1,k1,chi1)
    #             x = M.character()
    #             if N1 == 1:
    #                 si = 0
    #             else:
    #                 si = sage_character_to_sage_galois_orbit_number(x)
    #             clogger.debug("si={0}, chi={1}".format(si,chi))
    #             if si <> chi:
    #                 if (N1,k1,chi1) not in problems:
    #                     problems.append((N1,k1,chi1))
    #                 clogger.debug("in file: Chi is wrong! Got: {0} Should be {1}".format(si,chi))
    #                 if remove == 1:
    #                     dname = self._db.factor(N1,k1,chi1,d)
    #                     for fname in self._db.listdir(dname):
    #                         self._db.delete_file(fname)
    #                     aname = self._db.ambient(N1,k1,chi1)
    #                     dname = join(aname.split("/")[0:-1],"/")
    #                     for fname in self._db.listdir(dname):
    #                         self._db.delete_file(fname)
    #                     os.removedirs(dname)
    #                     clogger.debug("removed directory {0}".format(dname))
    #         if remove == 1 and problems<>[]:
    #             clogger.debug("Removed {0} records!".format(len(problems)))
    #     clogger.debug("Finished checking  character: N={0}, k={1}, chi={2}".format(N,k,chi))
    #     return problems



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
            q = self._newform_factors.find_one({'N':int(N),'k':int(k),'cchi':int(ci),'newform':int(d)})
        if q is None:
            raise ValueError,"This newform is not in the database!"

        N=q['N']; k=q['k']; cchi=q['cchi']; fid = q['_id']
        d = q['newform']
        chi = dirichlet_character_sage_from_conrey_character_number(N,ci)
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
                        #m,xi = sage_character_to_number(x) # get number for later reference
                        xi = sage_character_to_conrey_character(x).number()
                        ## Now check the characters
                        for y in dirichlet_group_sage(q):
                            if y.is_trivial():
                                continue
                            if DG(x)*DG(y**2)==chi:
                                print "possible character=",y
                                q = y.modulus()
                                yi = sage_character_to_conrey_character(y).number()
                                #q,yi = sage_character_to_number(y)
                                possible_twists.append((M,xi,q,yi)) # it is possible that we have a twist of F \in M_k(M,x) with character y of modulus q.
        bd = max(sturm_bound(N,k),2)
        print "Sturm bound=",bd
        aps = self.get_aps(N,k,ci,d,prec_needed=bd,coeffs=True)
        K0 = aps[0].base_ring()
        for M,xi,q,yi in possible_twists:
            print "Checking :",M,xi,q,yi
            y = dirichlet_character_sage_from_conrey_character_number(q,yi)
            for d1 in range(self.number_of_factors(M,k,xi)):
                apsf = self.get_aps(M,k,xi,d1,prec_needed=bd,coeffs=True)
                print "aps=",apsf
                twisted_aps = [ apsf[prime_pi(p)-1]*y(p) for p in primes(bd+1)]
                K1 = twisted_aps[0].base_ring()
                # if K0 == QQ:
                #     K = K1
                # elif K1 == QQ:
                #     K = K0
                # else:
                #     K = K1.extend(K0.polymomial)
                # print "twisted_aps=",twisted_aps
                # print "K=",K
                # test = [K(twisted_aps[i]) - K(aps[i]) for i in range(len(twisted_aps))]
                ## Since we don't know if we had the correct Galois conjugates
                ## we check as complex numbers
                if K0 == QQ:
                    aps_t = aps; aps_n = aps
                else:
                    aps_t = [aps[i].trace() for i in range(len(twisted_aps))]
                    aps_n = [aps[i].norm() for i in range(len(twisted_aps))]
                if K1 == QQ:
                    twisted_aps_t = twisted_aps
                    twisted_aps_n = twisted_aps
                else:
                    twisted_aps_t = [twisted_aps[i].trace() for i in range(len(twisted_aps))]
                    twisted_aps_n = [twisted_aps[i].norm() for i in range(len(twisted_aps))]
                    

                test_t = [twisted_aps_t[i] - aps_t[i] for i in range(len(twisted_aps))]
                test_n = [twisted_aps_n[i] - aps_n[i] for i in range(len(twisted_aps))]                
                print "test_t=",test_t
                print "test_n=",test_n                
                if test_t.count(0)==len(test_t) and test_n.count(0)==len(test_n):
                    print "We have a twist of {0} by {1}!".format((M,k,xi,d1),(q,yi))
                    return (M,k,xi,d1),(q,yi)
        return []

    def check_all_twists(self,verbose=0):
#        fs_twist = gridfs.GridFS(self._mongodb,self._twists_collection)
        for r in self._newform_factors.find().sort('N'):
            N=r['N']; k=r['k']; ci=r['cchi']; d=r['newform']
            if is_squarefree(N):
                continue
            t = self.check_if_twist(N,k,ci,d)
            if len(t)>0:
                labela = utils.newform_label(N,k,ci,d)
                labelb = utils.newform_label(t[0][0],t[0][1],t[0][2],t[0][3])
                clabel = "{0}.{1}".format(t[1][0],t[1][1])
                twist_rec = {'hecke_orbit_label':labela,'twist_of':labelb,'by_char':clabel}
                if verbose>0:
                    print "inserting; ",twist_rec
                self._twists.update(twist_rec)



    def insert_raw_data(self,pattern=""):
        r"""
        Go through modular forms data in the files and insert everything into mongo database in primitive / raw format, i.e. as a dict
        instead of a higher level object.
        pattern is a search pattern for the sql database
        """
        from utils import orbit_label
        
        mdb = self._mongo_conn[self._db_raw]
        mdb_ambient_files = mdb["ambient_data.files"]
        fs_a = gridfs.GridFS(mdb,'ambient_data')
        mdb_factor_files = mdb["factor_data.files"]
        fs_f = gridfs.GridFS(mdb,'factor_data')
        s = pattern
        ### Assuming we have changed the 'known' db to work with conrey galois orbits
        for (N,k,on,newforms,nap) in self._db.known(s):
            print N,k,on,newforms,nap
            orbit = dirichlet_group_conrey_galois_orbits_numbers(N)[on]
#            ii+=1
#            if ii>10:
#                return
            # insert the ambient space
            #sage_label = "{0}.{1}.{2}".format(N,k,ci)
            ci = conrey_character_from_galois_orbit_number(N,on).number()
            s = {'N':int(N),'k':int(k),'cchi':int(ci)}
            q = mdb_ambient_files.find_one(s)
            space_label = q['space_label']
            #conrey_galois_number = on
            #N, ci = conrey_character_number_to_conrey_galois_orbit_number(N,ci)
            #conrey_char = dirichlet_character_conrey_from_sage_character_number(N,i)
            #conrey_char = conrey_character_from_number(N,ci)
            #conrey_character_number = ci
            #conrey_char.number() #character_conversions.dirichlet_character_conrey_used_in_computation(N,conrey_char.number())
            if q is None: # insert it
                #print N,i
                #conrey_label = "{0}.{1}.{2}".format(N,k,conrey_galois_number)
                ambient_fname = self._db.ambient(N,k,on)
                try:
                    ambient = load(ambient_fname)
                    F = ambient['rels'].parent().base_ring() 
                    if F == QQ:
                        eps = DirichletGroup(N)(ambient['eps'])
                    else:
                        eps = DirichletGroup(N,F)(ambient['eps'])
                    ci = sage_character_to_conrey_character(eps).number()
                    t = ambient['space'] # Space should contain the correct character.
                    if t[2] <> ci:
                        clogger.critical("Space {0} not yet fixed!".format(ambient_name))
                    fname = ambient_fname.split("/")[-2]
                    space_label = "{0}.{1}.{2}".format(N,k,ci)
                    orbit_label = "{0}.{1}.{2}".format(N,k,on)
                    if ambient <> None:
                        fs_a.put(dumps(ambient),filename=fname,
                                 space_label=space_label,
                                 space_orbit_label = orbit_label,
                                 level=int(N),
                                 weight=int(k),
                                 character_galois_orbit = orbit,
                                 #chi=int(i),
                                 cchi=int(conrey_character_number))
                        clogger.debug("inserted ambient: {0} / {1}".format(sage_label,conrey_label))
                except IOError:
                    clogger.debug("Space {0} is not in files at {1}".format((N,k,no),ambient_fname))
            for newform in range(newforms):
                newform_label = "{0}.{1}".format(space_label,orbit_label(newform))

                q = mdb_factor_files.find_one({'newform_label':newform_label})
                if q is None: # insert newform
                    #conrey_newform_label="{0}.{1}.{2}{3}".format(N,k,conrey_galois_number,orbit_label(newform))
                    factor_fname = self._db.factor(N,k,on,newform)
                    try:
                        B = load(self._db.factor_basis_matrix(N, k, on, newform))
                        Bd = load(self._db.factor_dual_basis_matrix(N, k, on, newform))
                        v = load(self._db.factor_dual_eigenvector(N, k, on, newform))
                        nz = load(self._db.factor_eigen_nonzero(N, k, on, newform))

                        if B._cache is None:
                            B._cache = {}
                        if Bd._cache is None:
                            Bd._cache = {}
                        B._cache['in_echelon_form'] = True
                        Bd._cache['in_echelon_form'] = True
                        factor = {'B':B,'Bd':Bd,'v':v,'nz':nz}
                        fname = factor_fname.split("/")[-2]
                        if factor <> None:
                            fs_f.put(dumps(factor),filename=fname,
                                     newform_label=newform_label,
#                                     conrey_newform_label=conrey_newform_label,
                                     level=int(N),
                                     weight=int(k),
                                     character_galois_orbit = orbit,
#                                     chi=int(i),
                                     cchi=int(conrey_character_number),
                                     newform=int(newform))
                        clogger.debug("inserted newform:{0} / {1}".format(sage_newform_label,conrey_newform_label))

                    except IOError:
                        clogger.debug("Data is incomplete for factor ({0}) at {1}".format((N,k,on,newform),factor_fname))


    def check_all_characters2(self,Nmin=4):
        r"""

        """
        for N in self._newform_factors.distinct('N'):
            if N < Nmin:
                continue
            for k in self._newform_factors.find({'N':N}).distinct('k'):
                l = character_conversions.dirichlet_character_sage_galois_orbits_reps(N)
                l1 = filter(lambda x : x(-1)==(-1)**k,l)
                o = len(l1)
                chis = self._newform_factors.find({'N':N,'k':k}).distinct('cchi')
                if len(chis) >= o:
                    print N,k,o,chis
                else: print N,k,o,chis

    def check_all_characters(self,typec='ambient',dry_run=0):
        r"""
        Chck all characters in the database of ambient modular symbol spaces.
        Note that the update only works if there is currently no inddexes on the
        sollection self._modular_symbols but this should be recreated afterwards.
        """
        from utils import label_from_param
        from character_conversions import dirichlet_character_sage_galois_orbits_reps,conrey_character_number_to_conrey_galois_orbit_number,sage_character_to_sage_galois_orbit_number, sage_character_to_conrey_character,conrey_character_number_to_conrey_galois_orbit_number
        cnt = 0
        if typec=='ambient':
            for r in self._modular_symbols.find({'space_label':{"$exists":False}}):
                t = self.check_characters_ambient(r['N'],r['k'],r['chi'],dry_run=dry_run,verbose=0)
                if t > 0:
                    cnt+=t
        else:
            for r in self._modular_symbols.find():
                aid = r['_id']
                N = r['N']
                for f in self._newform_factors.find({'ambient_id':aid}):
                    fid = f['_id']
                    # Make sure again we have the correct character....
                    factor = self.load_from_mongo('Newform_factors',fid)
                    x = factor.character()
                    assert sage_character_to_sage_galois_orbit_number(x)==f['chi']
                    assert f['cchi'] == sage_character_to_conrey_character(x).number()
                    newlabel = label_from_param(r['N'],r['k'],r['cchi'],f['newform'])
                    if f['cchi']<>r['cchi'] or f['character_galois_orbit']<>r['character_galois_orbit'] or newlabel <> f['hecke_orbit_label'] or f.get('conrey_galois_orbit_number') is None:
                        newfname = "gamma0-factors-{0}".format(f["filename"].split("/")[-1])
                        on = conrey_character_number_to_conrey_galois_orbit_number(N,f['cchi'])
                        updates = {
                            "cchi":r['cchi'],
                            "character_galois_orbit":r['character_galois_orbit'],
                            'hecke_orbit_label' : newlabel,
                            "filename":newfname,
                            "conrey_galois_orbit_number":on}
                        self._newform_factors.update({'_id':fid},{"$set":updates})
                        cnt+=1
                                                     
                            
                
        print "Updated {0} records!".format(cnt)

    def check_characters_ambient(self,N,k,i,verbose=0,dry_run=0):
        from character_conversions import dirichlet_character_sage_galois_orbits_reps,conrey_character_number_to_conrey_galois_orbit_number
        q = self._modular_symbols.find({'N':int(N),'k':int(k),'chi':int(i)})
        if q is None:
            return 0
        cnt = 0
        for r in q:
            t = self.check_characters_one_ambient(r,verbose,dry_run)
            if t is True:
                cnt +=1
        return cnt
    
    def check_characters_one_ambient(self,r,verbose=0,dry_run=0):
        from character_conversions import dirichlet_character_sage_galois_orbits_reps,conrey_character_number_to_conrey_galois_orbit_number
        #r = self._modular_symbols.find_one({'N':int(N),'k':int(k),'chi':int(i)})
        cchi = r['cchi']
        fid = r['_id']
        N = r['N']; k=r['k']; chi=r['chi']
        print "Get ambient with ",N,k,cchi
        M = self.get_ambient(N,k,cchi)
        if M is None:
            # This space is probably empty due to inconsistency.
            # Check the character we say that we have
            x =  dirichlet_character_sage_galois_orbit_rep_from_number(N,i)
            cx =  dirichlet_character_conrey_from_sage_galois_orbit_number(N,i)
            cxx = conrey_character_from_number(N,cchi)
            if cx.sage_character() == x and cx==cxx:
                return True
            ci = cx.number()
            
            orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
            label = "{0}.{1}.{2}".format(N,k,ci)
            updates =  {'character_galois_orbit':orbit,
                        'cchi':ci,
                        "hecke_orbit_label":label}
            if verbose>0:
                clogger.debug("Updating {0} -> {1}".format(r['hecke_orbit_label'],label))
                clogger.debug("Updates: {0}".format(updates))
            if dry_run == 1:
                return False
            self._modular_symbols.update({'_id':fid},{"$set":updates})

            return False
            #raise ValueError,"Space with N,k,i={0} does not exist in the database!".format((N,k,i))
        # If the space is non-empty and in the datbase we make an extra check that
        # we indeed have the character which was used in the space
        #fname = 'gamma0-ambient-modsym-00101-002-006'
        x1 = M.character()
        if verbose>0:
            clogger.debug("x1={0}".format(x1))
        reps = dirichlet_character_sage_galois_orbits_reps(N)
        si = None
        for j in range(len(reps)):
            if x1 in reps[j].galois_orbit():
                si = int(j)
                break
        if si <> r['chi']:
            clogger.debug("Wrong character number! chi={0} and orbit of x={1}".format(r['chi'],si))
            if si is None:
                raise ArithmeticError,"Could not find correct character for N,k,i={0},{1},{2}".format(N,k,i)
        # Need to get the Conrey character number
        
        for x2 in dirichlet_group_conrey(N):
            if x1 == x2.sage_character():
                ci = x2.number()
                on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)
                clogger.debug("on: {0} and r={1}".format(on,r))
                if ci == cchi and si == r['chi'] and r.get('conrey_galois_orbit_number',-1)==on:
                    return True
                orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(N,ci)
                # we might need to update the filename as well
                label = "{0}.{1}.{2}".format(N,k,ci)
                orbit_label = "{0}.{1}.{2}".format(N,k,on)
                updates = {
                    'character_galois_orbit':orbit,
                    'cchi':ci,
                    'chi':si,
                    'nfactors':r['orbits'],
                    "space_label":label,
                    "space_orbit_label":orbit_label,
                    "conrey_galois_orbit_number":on}
                if verbose>0:
                    clogger.debug("Updating {0} -> {1}".format(r['hecke_orbit_label'],label))
                    clogger.debug("updates: {0}".format(updates))
                if dry_run == 1:
                    if verbose>0:
                        clogger.debug("We do not modify the database!")
                    return False
                self._modular_symbols.update({'_id':fid},{"$set":updates})
                return False
        raise ArithmeticError,"Could not find appropriate Conrey character!"


    def create_galois_orbits_maps(self,nmax):
        r"""
        Set up a database with correspondences between the two ordering of Galois orbits.
        Records are of the form:
        ## A character is given by e.g.:
        {'N':63,'values':['1','zeta6'],'conrey_number':i,
        'orbit':[],
        'orbit_no':j
        'galois_orbit_sage':k }
        {'N','sage_
        """
        raise NotImplementedError,"Haven't done this yet"
        for N in range(2,nmax):
            if self._galois_orbits.find({'N':int(N)}).count()>0:
                continue
            D = dirichlet_group_sage(N)
            D = dirichlet_group_conrey(N)            
            reps = D.galois_orbits(reps_only=True)
            for x in reps:
                for y in x.galois_orbit():
                    pass

    def get_record_from_character(self,x):
        r"""
        Input a Sage character x and output a record for database.
        """
        
        DC = dirichlet_group_conrey(x.modulus())
        #vals = map(str,x.values_on_gens())
        # The generators of (Z/NZ)^* given by sage
        gens = list(IntegerModRing(c.modulus()).unit_group().gens_values())
        gens.sort()
        # To make sure that we have the values of the character on these generators...
        vals = map(x,gens)
        clogger.debug("Gens of Z/{0}Z:{1}".format(x.modulus(),gens))
        clogger.debug("Values : {2}".format(vals))
        rec={'N':int(x.modulus()),'gens':gens,'vals':vals}

            
    def check_characters_in_files(self,nmin=1,nmax=10000,verbose=0):
        from sage.all import trivial_character,DirichletGroup
        from sage.all import dimension_new_cusp_forms
        from character_conversions import sage_character_to_conrey_galois_orbit_number,sage_character_to_conrey_character
        from dirichlet_conrey import DirichletGroup_conrey
        import os
        rename_list = []
        missing = []
        for N,k,i,d,ap in self._db.known("N<{1} AND N>{0}".format(nmin,nmax)):
            ## find the character in file...
            if verbose>1:
                print "Checking ",N,k,i
            sname = self._db.space(N,k,i)
            if not self._db.isdir(sname):
                if verbose>0:
                    clogger.debug("Space {0} is missing!".format((N,k,i)))
                continue
            mname = self._db.ambient(N,k,i)            
            try:
                modsym = load(mname)
            except IOError:
                if k % 2 == 0 and i==0:
                    missing.append((N,k,i))
                    clogger.debug("Space {0} is missing at {1}!".format((N,k,i),sname))
                continue
            except ValueError as e:
                clogger.critical("Could not load {0} at {1} due to:{2}".format((N,k,i),sname,e.message))
                f = open('needs_recomputation', 'a')
                f.write("{0},{1},{2}".format(N,k,i))
                f.close()
                continue
            rels  = modsym['rels']
            F = rels.base_ring()
            if i == 0:
                eps = trivial_character(N)
                # this is always ok
                continue
            if F==QQ:
                eps = DirichletGroup(N)(modsym['eps'])
            else:
                eps = DirichletGroup(N,F)(modsym['eps'])
            # Instead of constructing the DirichletgGroup with base_ring since this
            # messes up the conversion to Conrey characters we make sure afterwards...
            if verbose>0:
                print "eps[sage]=",eps
            conrey_eps = sage_character_to_conrey_character(eps)
            if conrey_eps.sage_character()<>eps:
                clogger.critical("We do not get the correct character for {0}".format(mname))
                continue
            conrey_i = conrey_eps.number()
            #clogger.debug("eps={0}".format(eps))
            NN,j = conrey_character_number_to_conrey_galois_orbit_number(N,conrey_i)
            #NN,j = sage_character_to_conrey_galois_orbit_number(eps)
            #if j == i:
            #    clogger.debug("File {0} is ok!".format(mname))
            #    continue
#            clogger.debug("eps={0} \t conrey_eps = {1} conrey_i={2},\t conrey_gal_nr={3}".format(eps,conrey_eps,conrey_i,j))
#            clogger.debug("Conrey eps:{0}".format(DirichletGroup_conrey(N).from_sage_character(eps)))
#            clogger.debug("Sage eps.values:{0}".format(eps.values()))
#            clogger.debug("Conrey eps.values:{0}".format(conrey_eps.values()))
            
            mnamenew = self._db.ambient(N,k,j)
            snamenew = self._db.space(N,k,j)
            snamenew = snamenew+"-c"
            t = (int(N),int(k),int(conrey_i))
            if modsym['space'] <> t:
                modsym['space'] = t
                save(modsym,mname) # save wih updated space name
            if self._db.isdir(snamenew):
                clogger.critical("\t Directory {0} already exists!".format(snamenew))
            else:
                clogger.debug("Will change filename from {0} to {1}".format(sname,snamenew))
                os.rename(sname,snamenew)
                rename_list.append([sname,snamenew])            
        print "Need to change name of {0} directories!".format(len(rename_list))
        return missing,rename_list

    def _remove_c_from_filename(self):
        
        for N,k,on,d,ap in self._db.known(""):
            mname = self._db.space(N,k,on)
            tmp = mname + "-c"
            if self._db.isdir(mname):
                clogger.critical("Could not move {0} to {1}: Already existing!",format(tmp,mname))
            if self._db.isdir(tmp):
                os.rename(tmp,mname)
                clogger.debug("Renamed {0} --> {1}".format(tmp,mname))
                  

def precision_needed_for_L(N,k,**kwds):
    r"""
    Returns the precision (number of coefficients) needed to compute the first zero
    of the L-function (as on the LMFDB pages). This bound is taken from there and is probably heuristic.
    """
    from sage.all import ceil
    pprec = 20 + int(RR(5) * ceil(RR(k) * RR(N).sqrt()))
    pprec = max(pprec,kwds.get('pprec',100))
    ## Get even hundreds of primes to look nicer.
    return ceil(RR(pprec)/RR(100))*100



class CheckingDB(CompMF):
    r"""
    Class which is used for checking and fixing data in the Mongo Database
    
    """
    def __init__(self,datadir='',host='localhost',port=37010,db='modularforms2',user='',password='',verbose=0,**kwds):
        r"""
'',
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
        super(CheckingDB,self).__init__(datadir,host,port,db,user=user,password=password,verbose=verbose,**kwds)





    def find_records_needing_completion(self,nrange=[],krange=[],cchi=None,check_content=False,recheck=False,ncpus=1):
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
        if cchi=='trivial' or cchi==1:
            s['cchi'] = int(1)
        elif isinstance(cchi,list):
            s['cchi'] = {"$in":map(int,cchi)}
        args = []
        clogger.debug("search  pattern :{0}".format(s))
        for r in self._modular_symbols.find(s):
            N = r['N']; k=r['k']; ci = r['cchi']
            clogger.debug("r = {0}".format((N,k,ci)))
            args.append((N,k,ci,check_content,recheck))
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
            norbits=len(dirichlet_character_conrey_galois_orbits_reps(n))
            for k in krange:
                for i in range(norbits):
                    if self._modular_symbols.find({'N':int(n),'k':int(k),'cchi':int(ci)}).count()==0:
                        res[(n,k,i)]=[False]

        return res

    def complete_records(self,nrange=[],krange=[],cchi=None,ncpus=1,check_content=False,recheck=False,from_files=True):
        r"""
        Check all records within a specified bound and update / compute the incomplete ones.

        INPUT:

        - nrange -- list/tuple : give upper and lower bound for the levels we complete
        - krange -- list/tuple : give upper and lower bound for the weights we complete
        - chi    -- integer    : if not None we only look at this character (e.g. for trivial character chi=0)

        - 'check_content' -- bool : set to True to make a more detailed check of completeness of records
        - 'recheck' -- bool: set to True if you want to recheck records already marked as complete.
        """
        recs = self.find_records_needing_completion(nrange,krange,cchi=cchi,check_content=check_content,recheck=recheck,ncpus=ncpus)
        args = []
        for N,k,ci in recs.keys():
            if k==1:
                clogger.debug("Weight 1 is not implemented!")
                continue
            if isinstance(cchi,(list,tuple)):
                if ci not in cchi:
                    continue
            if isinstance(cchi,(int,Integer)):
                if ci <> cchi:
                    continue
            if not are_compatible(N,k,ci):
                #clogger.debug("N,k,i={0} is incompatible!".format((N,k,i)))
                continue
            args.append((N,k,ci))
#        if from_files:
#            N,k,ci
        clogger.debug("Completing {0} spaces!".format(len(args)))
        self.get_or_compute_spaces(args,ncpus=ncpus,compute=True)
        return True
        
    
    def check_records(self,nrange,krange,irange='all',check_content=False,recheck=False,ncpus=8):
        r"""
        We check if the records corresponding to M(N_i,k_i,i_i) is complete or not.

        """
        args = []
        if not isinstance(nrange,(list,tuple)):
            nrange = [nrange]
        if not isinstance(krange,(list,tuple)):
            krange = [krange]
        if isinstance(irange,(int,Integer)):
            irange = [irange]
        if check_content:
            check_level = 2
        else:
            check_level = 1
        clogger.debug("check_record with level: {0}".format(check_level))
        s = {}
        if recheck is False:
            s['complete']={"$lt":check_level+int(1)}
        clogger.debug("s = {0}".format(s))
        for r in self._modular_symbols.find(s,projection=['N','k','cchi']).sort([('N',pymongo.ASCENDING),('cchi',pymongo.ASCENDING),('k',pymongo.ASCENDING)]):
            n = r['N']; k=r['k']; cchi=r['cchi']
            if n==7 and k==3:
                print r
            if nrange<>[] and (n < nrange[0] or n > nrange[-1]):
                continue
            if krange<> [] and (k < krange[-1] or k > krange[0]):
                continue
            if irange <> 'all' and irange<>[] and (cchi < irange[0] or cchi >irange[-1]):
                continue
            args.append((n,k,cchi,check_content,recheck))
        clogger.debug("args={0}".format(args))
        if ncpus >= 32:
            return list(self.check_record32(args))
        elif ncpus >= 16:
            return list(self.check_record16(args))
        elif ncpus >= 8:
            return list(self.check_record8(args))
        return list(self.check_record(args))

    @parallel(ncpus=8,verbose=True)
#    @classmethod
    def check_record8(self,N,k,ci,check_content=False,recheck=False):
        return self.check_record(N,k,ci,check_content,recheck)

    @parallel(ncpus=16,verbose=True)
#    @classmethod
    def check_record16(self,N,k,ci,check_content=False,recheck=False):
        return self.check_record(N,k,ci,check_content,recheck)
    @parallel(ncpus=32,verbose=True)
#    @classmethod
    def check_record32(self,N,k,ci,check_content=False,recheck=False):
        return self.check_record(N,k,ci,check_content,recheck)

    @parallel(ncpus=1,verbose=True)
#    @classmethod
    def check_record(self,N,k,ci,check_content=False,recheck=False):
        r"""

        We check if the record corresponding to M(N,k,i) is complete or not.

        """

        res = {}
        s = {'N':int(N),'k':int(k),'cchi':int(ci)}
        res = {}
        ### Check the ambient space
        check_level = int(3) if check_content else int(1)
        if not recheck:
            if self._modular_symbols.find({'N':int(N),'k':int(k),'cchi':int(ci),'complete':{"$gt":check_level-int(1)}}).count()>0:
                return  {'modular_symbols':True,'aps':True,'factors':True}
        clogger.debug("Checking N,k,ci={0}".format((N,k,ci)))
        ambient_id = None
        M = self.get_ambient(N,k,ci,compute=False)
        if M is None:
            res['modular_symbols']=False
            numf = 0
        else:
            self.check_characters_ambient(N,k,ci)
            if check_content and not 'modsym.ambient' in str(M.__class__):
                clogger.warning("Space is reconstructed with wrong class!")
                clogger.warning("type(M)={0}".format(type(M)))
                clogger.warning("type(MS)={0}".format(ModularSymbols(1,2).__class__))
                res['modular_symbols']=False
            else:
                res['modular_symbols']=True
            rec = self._modular_symbols.find_one({'N':int(N),'k':int(k),'cchi':int(ci)})
            clogger.debug("rec={0}".format(rec))
            numf = rec['orbits']
            ambient_id = rec['_id']
        ## Check factors
        clogger.debug("will check number of factors! for M={0}".format(M))
        numf1 = self.number_of_factors(N,k,ci)
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
                facts = self.get_factors(N,k,ci)
                d = 0
                for f in facts.values():
                    d+=f.dimension()
                clogger.debug("Sum of dimensions of factors: {0}".format(d))
                #if M is None:
                if ci > 1:
                    eps = dirichlet_character_sage_from_conrey_character_number(N,ci)
                    d1 = dimension_new_cusp_forms(eps,k)
                else:
                    d1 = dimension_new_cusp_forms(N,k)
                clogger.debug("Dimension of space is: {0}".format(d1))
                if d <> d1:
                    res['factors'] = False
                    clogger.warning("Dimensions of all factors do not sum up to the total dimension! n,k,cchi={0}".format((N,k,ci)))
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
        #    #aps = {(N,k,ci,x) : {pprec: True}  for x in newforms}
        if d==d1 and d==0:
            res['aps']=True
        else:
            res['aps']=False
        if not check_content and res['factors']:
            newforms_with_aps = self._aps.find({'N':int(N),'k':int(k),'cchi':int(ci)}).distinct('newform')
            res['aps'] = len(newforms_with_aps) == numf
        clogger.debug("facts={0}, numf={1}".format(facts,numf))
        fs_ap = gridfs.GridFS(self._mongodb,self._aps_collection)
        if res['factors'] is False:
            facts = {}
        for t in facts.keys():
            clogger.debug("t={0}".format(t))
            N,k,ci,d=t
            s = {'N':int(N),'k':int(k),'cchi':int(ci),'newform':int(d)}
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
                    s = {'N':r['N'],'k':r['k'],'cchi':r['cchi'],'newform':r['newform'],
                         'prec':prec_in_db}
                    if  self._aps.find(s).count()>0:
                        clogger.debug("We already have this prec in the database so we remove this record!")
                        self._aps.remove({'_id':r['_id']})
                    else:
                        fname = "gamma0-aplists-{0:0>5}-{1:0>3}-{2:0>3}-{3:0>3}-{4:0>3}".format(N,k,ci,d,prec_in_db)
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
            clogger.debug("Space {0},{1},{2} is not in database".format(N,k,ci))
        clogger.debug("done checking record!")
        sys.stdout.flush()
        return res
    

    
    def check_ramanujan(self):
        r"""
        Check Ramanujan bound for ap's in the database (simple chekc to catch trivially wrong stuff)
        """
        # first get
        from sage.all import primes_first_n
        for r in self._aps.find():
            N=r['N']; k=r['k']; cchi=r['cchi']
            fid = r['_id']
            kk = RR(k-1)/RR(2)
            aps = self.get_aps(N,k,cchi)
            for d in aps.keys():
                for prec in aps[d].keys():
                    E,v=aps[d][prec][0:2]
                    c = E*v
                    i = 0
                    for p in primes_first_n(len(c)):
                        t = c[i].abs()/p**kk
                        if abs(t)>2:
                            raise ArithmeticError,"ERROR: a({0})={t} for (N,k,chi)={1},{2},{3}, prec={4}".format(p,N,k,chi,prec,t=t)
                        i+=1
#    def _mul_nf_elts_appr(self,a,b):
#        ap = a._pari_()
#        bp = b._pari_()
#        return (ap*bp).norm().numerical_approx()

    def fix_coefficient_records(self,nlim=10,chi_in=0,check=None):
        from sage.all import next_prime,previous_prime,prime_pi
        for r in self._aps.find({"pmax":{"$exists":False},"N":{"$lt":nlim}}).sort([("N",int(1))]):
            fid = r['_id']
            N=r['N']; k=r['k']; cchi=r['cchi']
            if isinstance(check,tuple) and check==(N,k,cchi):
                pass
            else:
                if chi_in ==0 and cchi <>0:
                    continue
            kk = RR(k-1)/RR(2)
            fname = r['filename'].split("/")[-1]
            print "Checking:",N,k,chi,fname,fid
            l = fname.split("-")[2:]
            if filter(lambda x: not x.isdigit(),l)<>[] or len(l)<>5:
                print "File name in wrong format:{0}".format(fname)
                continue
            Ns,ks,ist,newform,nmax = l[-5:]
            nmax=int(nmax)
            pmax = previous_prime(nmax)
            nump = prime_pi(nmax) #- prime_pi(0)
            E,v = self.load_from_mongo('ap',fid)
            # check first and last coefficients
            #cpmin = sum([E[0,i]*v[i] for i in range(len(v))])
            #cpmax = sum([E[-1,i]*v[i] for i in range(len(v))])
            #t1 = cpmin.abs()/2**kk
            #t2 = cpmax.abs()/pmax**kk
            #if t2>2: # more seriously, the estimate of the last coefficient is wrong
            #    raise ArithmeticError,"ERROR: a({0})>2 for (N,k,chi)={1},{2},{3}, fname={4}".format(p2,N,k,chi,fname)
            if nump <> E.nrows() :
                # we probably need to change the name:
                print "fname=",fname
            else:
                print "fname=",fname," is ok! rows=",E.nrows()
                self._aps.update({'_id':fid},{"$set":{"pmax":int(pmax)}})
            print "checked ",fname

            
    def add_dimension_newforms(self,**kwds):

        ncpus = kwds.pop('ncpus',1)
        pool = Pool(processes=ncpus)
        clogger.debug("ncpus={0}".format(ncpus))     
        ambient_ids = []
        for r in self._modular_symbols.find({'dimn':{"$exists":False}}):
            ambient_ids.append({'space_label':r['space_label'],'aid':r['_id']})
        chunksize = 20
        results = pool.imap_unordered(self.add_one_dimn,ambient_ids,chunksize)
        for res in results:
            res.get()
        pool.close()
        pool.join()

    def add_one_dimn(self,r):
        M = self.load_from_mongo('Modular_symbols',r['aid'])
        if M is None:
            return
        clogger.debug("M={0}".format(r['space_label']))
        S = M.cuspidal_submodule().new_submodule()
        dimn = int(S.dimension())
        clogger.debug("dimn {0}  = {1}".format(r['space_label'],dimn))
        return self._modular_symbols.update({'_id':r['aid']},{"$set":{"dimn":dimn}})
