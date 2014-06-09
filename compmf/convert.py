r"""
Programs aimed to read Williams data and convert to/from other files or databases
"""
import sage.all
import os,re
import glob
import sqlite3
#import nosqlite
#nsql = nosqlite.Client('nsql')
import pymongo
try:
    MDB = pymongo.Connection('localhost:27017')
except:
    print "No MongoDB at 27017"
try:
    MDB1 = pymongo.Connection('localhost:37010')
except:
    print "No MongoDB at 37010"


from sage.all import (ModularSymbols, DirichletGroup, trivial_character,
                      dimension_new_cusp_forms,
                      save, load,
                      cputime,
                      fork, parallel,
                      Integer,
                      prime_range, prime_divisors,
                      version,
                      Sequence,
                      cached_function,
                      Integer)
#import mfdb.compute
#db = mfdb.compute.filenames
#import mfdb
#import mfdb.compute
import compute
from compute import *

class WDB(object):
    r"""
    Handles Williams database
    """
    def __init__(self,dir='data',dir2='',update=False,**kwds):
        #self._mfdb = mfdb
        #self._compute = mfdb.compute
        self._db = FilenamesMFDB(dir)
        self._db2 = None
        if update==True:
            self.update()
        if dir2<>'':
            print "Initing a second database at {0}".format(dir2)
            self._db2 = FilenamesMFDB(dir2) ## Can be readonly.
    def update(self):
        self._db.update_known_db()
            
    def _is_valid_sql(self):
        raise NotImplementedError
        
    def known(self,N='all',k='all',i='all',format='iterator'):
        r"""
        Check which data is available matching the parameters N,k,i
        Available formats: 'iterator', 'list' and 'unique'
        """
        if format=='iterator':
            return self._db.known(self._make_query(N,k,i))
        if format =='list':
            return list(self._db.known(self._make_query(N,k,i)))
        levels = []; weights=[]; chars=[]; n=0
        for N0,k0,i0,j0,a0 in self._db.known(self._make_query(N,k,i)):
            if N0 not in levels:
                levels.append(N0)
            if k0 not in weights:
                weights.append(k0)
            if i0 not in chars:
                chars.append(i0)
            n+=1
        if format=='data':
            return {'num_recs':n,'levels':levels,'weights':weights,'chars':chars}
        elif format=='web':
            s = "Numbers of records: {n} \n ".format(n=n)
            s+= "Levels in range {0}--{1} \n ".format(min(levels),max(levels))
            s+= "Weights in range  {0}--{1} \n ".format(min(weights),max(weights))
            if chars==[0]:
                s+= "Character: trivial"
            else:
                s+= "Character: non-trivial"
            return s
        else:
            raise ValueError,"Incorrect format:{0}".format(format)
        
    def get_spaces(self,N='all',k='all',i='all',n=None,format='web'):
        r"""
        Fetch all spaces matching the search parameters
        If init_spaces = False we only return a list of file names.

        If format = 'web' then return string representation to display.
        If format = 'data' then return modular symbols objects
        
        """
        q = self._make_query(N,k,i)
        l = self._db.known(q)
        res = []
        for N,k,i,numorbits,naps in l:
            ### We have a modular form with this data.
            res0={}
            if format=='web':
                #s = self.space_info_web(N,k,i,numorbits)
                if numorbits==0:
                    s = "S^new({N},{k},{i}) is empty!".format(N=N,k=k,i=i)
                elif numorbits==1:
                    s = "S^new({N},{k},{i}) has {d} orbit with {n} coefficients".format(N=N,k=k,i=i,d=numorbits,n=naps)
                else:
                    s = "S^new({N},{k},{i}) has {d} orbits with {n} coefficients each".format(N=N,k=k,i=i,d=numorbits,n=naps)
                res.append(s)
            else:
                fname = self._db.ambient(N,k,i)
                l = load(fname)
                if isinstance(l,dict):
                    l = [l]
                for Md in l:
                    if not isinstance(Md,dict):
                        raise ValueError,"Need to have a dict here!"
                    M = self._db.dict_to_ambient(Md)
                    ##  M is an ambient space
                    res0['ambient_dict']=Md
                    res0['ambient']=M
                    res0['num_orbits']=numorbits
                    res0['orbits']=[]
                    for d in range(numorbits):
                        A = self._db.load_factor(N,k,i,d,M)
                        
                        res0['orbits'].append(A)
                res.append(res0)
            #return res0
        return res

    
    
    def get_decomposition(self,N,k,i):
        r"""
        Fetch tuples (B,Bd,v,nz)
        """        
        res = {}
        for N,k,i,dmax,aps in self._db.known(self._make_query(N,k,i)):
            l=[]
            res[(N,k,i)]={}
            for d in range(dmax):
                B = load(self._db.factor_basis_matrix(N, k, i, d))
                Bd =load(self._db.factor_dual_basis_matrix(N, k, i, d))
                v = load(self._db.factor_dual_eigenvector(N, k, i, d))
                nz =load(self._db.factor_eigen_nonzero(N, k, i, d))
                res[(N,k,i)][d]={'B':B,'Bd':Bd,'v':v,'nz':nz}
        return res

    def get_aps(self,N0='all',k0='all',i0='all',n0=None,format='web',all_coeffs=True):
        r"""
        Fetch aps for all spaces matching the search parameters
        If init_coeffs = False we only return a list of file names.

        
        """
        q = self._make_query(N0,k0,i0)
        l = self._db.known(q)
        res = {}
        rb = "[0-9][0-9][0-9][0-9][0-9]"
        print "l=",l
        for N,k,i,numorbits,aps in l:
            #print N,k,i,numorbits,aps
            #res[(N,k,i)]={}
            for j in range(numorbits):
                if n0 <>None and j<>n0:
                    continue
                #print "aps=",aps
                v = compute.load(self._db.factor_dual_eigenvector(N,k,i,j,makedir=False))
                cdir = self._db.factor(N,k,i,j)
                print "cdir=",cdir
                s = '{dir}/aplist-*{{rb}}.sobj'.format(dir=cdir,rb=rb)
                print "glob_Str=",s
                list_of_files = glob.glob(s)
                print "list=",list_of_files
                for name in list_of_files:
                    print "name=",name
                    apn = re.findall("aplist.*",name)[0]                
                    ns = re.findall("[0-9]+",apn)
                    print "apn=",apn
                    print "ns=",ns
                    if len(ns)==1:
                        n=Integer(ns[0].lstrip('0')); m=2
                        fname = '{dir}/aplist-{n:05d}.sobj'.format(dir=cdir,n=n)      
                    elif len(ns)==2:
                        m= Integer(ns[0].lstrip('0')); n = Integer(ns[1].lstrip('0'))
                        fname = '{dir}/aplist-{m:05d}-{n:05d}.sobj'.format(dir=cdir,m=m,n=n)
                    else:
                        raise ValueError,"Got wrong filename format!"

                    aplist = load(fname)                             
                    meta = load(self._db.meta(fname))
                    if format == 'web':
                        s += self._format_aps(aplist,v)                             
                        s += "\n"+str(meta)
                        res[(N,k,i,j)]=s                         
                    else:
                        res[(N,k,i,j)]={(m,n):(aplist*v,meta)}
                    if not all_coeffs:
                        break
 ## initial list of ap's
                ## Getting the first file
                #fname = self._db.factor_aplist(N,k,i,j,False,100)
                #meta = load(self._db.meta(fname))
                #aplist = compute.load(fname)
                
                ## Do we want more coefficients or just the first?
                # if format == 'web':
                #     s = self._format_aps(aplist,v)
                #     s += "\n"+str(meta)
                #     res[(N,k,i)][j]=s
                # else:
                #     #aplist*v
                #     res[(N,k,i)][j]={(0,100):(aplist*v,meta)}
                # # if all_coeffs:
                #     cdir = self._db.factor(N,k,i,j)
                #     list_of_files = glob.glob('{dir}/aplist-{rb}-{rb}.sobj'.format(dir=cdir,rb=rb))
                #     #print "list=",list_of_files
                #     for name in list_of_files:
                #         #print "name=",name
                #         apn = re.findall("aplist.*",name)[0]                
                #         ns = re.findall("[0-9]+",apn)
                #         m,n =ns
                #         fname = '{dir}/aplist-{m:05d}-{n:05d}.sobj'.format(dir=cdir,m=m,n=n)
                #         aplist = mfdb.compute.load(fname)                             
                #         meta = load(mfdb.compute.filenames.meta(fname))
                #         if format == 'web':
                #              s += self._format_aps(aplist,v)                             
                #              s += "\n"+str(meta)
                #              res[(N,k,i)][j]=s 
                            
                #         else:

                #         res[(N,k,i,j)]={(m,n):(s,meta)}
            return res


    def insert_records_to_external_db(self,N,k,i,ext_db,format='full'):
        r"""
        Take records from self._db and insert information into external database
        Available formats: 'full' 
        """
        #q = self._make_query(N,k,i)
        f = self.known(N,k,i)
        for r in f:
            self.insert_one_record_in_external_db(N,k,i,ext_db,format)


    def insert_one_record_in_external_db(self,N,k,i,ext_db,format):
        if format=='full':
            space = self.get_spaces(N,k,i)
            aps = self.get_aps(N,k,i)[(N,k,i)]
            dec = self.get_decomposition(N,k,i)[(N,k,i)]
            #space_record = self.format_space_ext(space)
            rec['num_orbits']=len(dec.keys())
            rec['q_exp']
            
            
    ## Internal helper functions
        
    def _format_aps(self,aplist,v):
        r"""
        String representation of the ap'"""
        K = v.parent().base_ring()
        if hasattr(K, 'defining_polynomial'):
            s = "{0}".format(K.defining_polynomial().list())
        else:
            s = "1"  ## Defining polynomial for Q
        ap = str(list(aplist*v)).replace(' ','')
        if K.absolute_degree() >= 4:
            ap = ap.replace(',',',\n')
        s += '; {ap}'.format(ap=ap)
        return s

    def _make_query(self,N='all',k='all',i='all'):
        r"""
        Make_query a query string from N,k and i
        If N is a string and N<>'all' we assume that N is an SQL query
        """
        if isinstance(N,str) and N<>'all':
            return N
        qN="1";qk="1";qi="1"
        if N<>'all':
            qN = "N={N}".format(N=N)
        if k<>'all':
            qk = "k={k}".format(k=k)
        if i<>'all':
            qi = "i={i}".format(i=i)
        q = "{qN} and {qk} and {qi}".format(qN=qN,qk=qk,qi=qi)
        return q



### Testing alternative storage of number field elements as plain text

ITMAX = 10000  ## Avoid infinite recursion
def random_polynomial(degb=5,cb=100,irreducible=True,it=0):
    r"""
    Return a random (irreducible) polynomial of degree <= degb and
    with coefficients bounded by cb in absolute value.
    """

    if it>ITMAX:
        raise Arithmeticerror,"Could not find irreducible polynomial within {0} iterations!".format(ITMAX)
    d = ZZ.random_element(2,degb)
    x = ZZ['x'].gens()[0]
    p = x**d
    
    for i in range(d-1):
        c = ZZ.random_element(-cb,cb)
        if c==0:
            c=1
        p = p + c*x**i

    if irreducible:
        if not p.is_irreducible():
            return random_polynomial(degb,cb,irreducible=True,it=it+1)
    return p

def random_number_field(degb=5,discb=5000,cb =100,it=0):
    r"""
    Return random number field with discriminant bounded by discb.
    """
    if it>ITMAX:
        raise Arithmeticerror,"Could not find number field with disciminant bounded by {0} in {1} iterations!".format(discbd,ITMAX)
    p = random_polynomial(degb,cb,True)
    if p.discriminant()>discb:
        return random_number_field(degb,discb,cb,it=it+1)
    return NumberField(p,names='x')
    
def random_list_of_number_field_elts(N,degb=5,discb=5000,cb=100):
    res = []
    F = random_number_field(degb,discb,cb)
    for i in range(N):
        res.append(F.random_element())
    return res
        

def convert_list(l):
    res = []
    for x in l:
        res.append(NF_to_str(x))
    return str(res)

def convert_dict(d):
    res = dict()
    for x in d.keys():
        res[d]=append(NF_to_str())
    return str(res)



def NF_to_str(x):
    K = x.parent()
    v = x.vector()
    #p = K.polynomial().list()
    s = "{0}:{1}".format(K.polynomial().list(),x.vector())
    return s
