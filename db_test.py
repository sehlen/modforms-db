r"""
Inserting records in the database from ?

EXAMPLES::

    sage: C=ComputeMFData('/home/purem/cvzx53/Programming/git/modforms-db/data')
    # Compute ambient space of newforms of level 7, weight 12 and all characters
    sage: C.compute_ambient_space(7,12,'all',1)
    # Compute the decomposition of these spaces into Galois orbits
    sage: C.compute_decompositions(7,12,'all')
    # Compute lists of Fourier coeficients
    sage: C.compute_aplists(7,12,'all')
    
"""
import inspect,os

basedir =  os.path.dirname(inspect.getabsfile(inspect.currentframe()))
os.sys.path.append("{0}/../mfdb/".format(basedir))
os.sys.path.append("{0}".format(basedir))
os.sys.path.append("{0}/mdb/".format(basedir))

import mfdb
import mdb
from mdb import schema
from schema import ModularSymbols_ambient_DB_class,Coefficient_DB_class,NumberField_DB_class,ModularSymbols_base_field_DB_class,CoefficientField_DB_class,ModularSymbols_oldspace_factor_DB_class,ModularSymbols_newspace_factor_DB_class,ModularSymbols_base_field,CoefficientField,Coefficient_DB
from nf_schema import NumberField_DB,AlgebraicNumber_DB,AlgebraicNumber_DB_class
from conversions import extract_args_kwds

#from mdb import schema_sage 
#from schema_sage import ModularSymbols_ambient,ModularSymbols_newspace_factor,ModularSymbols_oldspace_factor,Coefficient,NumberField,ModularSymbols_base_field, CoefficientField,AlgebraicNumber
#schema.setup_all() 
#schema.create_all()
#DB=mfdb.WDB('git/mfdb/data/')
from mfdb import WDB,ComputeMFData
from mfdb import FilenamesMFDB

from mdb.db import db
#print DB.known(format='web')
#def do_computations_ranges(
from sage.all_cmdline import *   # import sage library
if not os.path.exists('data'):
    os.makedirs('data')
import glob, os, os.path,re, sys
import pymongo
from pymongo import Connection
import gridfs


class WDBtoMFDB(WDB):
    r"""
    Class to pull records from database in William's format and insert in our database.

    
    """
    def __init__(self,datadir,verbose=0,**kwds):
        self._dir = datadir
        super(WDBtoMFDB,self).__init__(dir=datadir,**kwds)
        self._dbb = mdb.db.db 
        if verbose>0:
            self._dbb.session.bind.echo = True
        else:
            self._dbb.session.bind.echo = False
        #self._dbb.setup_all() 
        self._dbb.create_all()       

    def source_db(self):
        r"""
        Return a description of the source database.
        """
        print "Source DB = Files at {0}".format(self._dir)

    def target_db(self):
        r"""
        Return a description of the target database.
        """
        print "Target DB= {0}".format(self._dbb.metadata.bind)


        
    def insert_spaces(self,q="N=1 and k=12"):
        r"""
        Insert spaces matching query.
        """
        if q=='all':
            q = ""
        for level,weight,character,numo,nap in self.known(q):
            ## First check if record already exists.
            level=int(level); weight=int(weight); characer=int(character)
            if ModularSymbols_ambient.query.filter_by(level=level,weight=weight,character=character).count()>0:
                continue
            print "Inserting {0},{1},{2},{3},{4}".format(level,weight,character,numo,nap)
            d = self.get_spaces(level,weight,character,format='data')[0]
            M = d['ambient']
            d['level']=level; d['weight']=weight; d['character']=character
            orbits = d['orbits']
            assert d['num_orbits']==numo
            d['orbits_dict'] = self.get_decomposition(level,weight,character)[(level,weight,character)]    
            if character == 0 or character == 1:
                d['dimension_modular_forms']=dimension_modular_forms(level,weight)
                d['dimension_cusp_forms']=dimension_cusp_forms(level,weight)
                d['dimension_new_cusp_forms']=dimension_new_cusp_forms(level,weight)
            else:
                x = DirichletGroup(level)[character]
                d['dimension_modular_forms']=dimension_modular_forms(x,weight)
                d['dimension_cusp_forms']=dimension_cusp_forms(x,weight)
                d['dimension_new_cusp_forms']=dimension_new_cusp_forms(x,weight)                
            print "d=",d

            self.insert_space_into_new_db(d)

    def insert_coefficients(self,q="N=1 and k=12"):
        r"""
        Insert coefficients for the spaces matching the query
        """
        if q=='all':
            q = ""
        for level,weight,character,numo,nap in self.known(q):
            if nap == 0:
                continue
            q = ModularSymbols_ambient.query
            q = q.filter_by(level=level,weight=weight,character=character)
            if q.count()==0:
                continue
            A = q[0]
            numfac = len(A.newspace_factors)
            assert numfac == numo
            for i in range(numfac):
                aps = A.get_aps(level,weight,character,i)[(level,weight,character)][0]
                ## Make an's from ap's by builtin methods
                N = A.newspace_factors[i]
                if not hasattr(N,_HeckeModule_free_module__eigenvalues):
                    evs = {}
                else:
                    evs = N._HeckeModule_free_module__eigenvalues
                evs[1]={'alpha':1}
                list_of_primes = primes_first_n(len(aps))
                for i in len(list_of_primes):
                    if not has_key(evs,p):
                        evs[list_of_primes[i]]={'alpha': ap[i]}
                ## Now we have ap's and want to set an's
                N._HeckeModule_free_module__eigenvalues = evs
#                for n in range(100):
#                    if not is_prime(n):                   
                for n in range(100):
                    c = Coefficient(N.eigenvalue(n))
                    A.coefficients.append()
            
            

    def insert_ambient_spac_into_new_db(self,M):
        r"""
        Insert the ambient space into the database.
        """
        
            

    def insert_space_into_new_db(self,M):
        r"""
        Insert M into the database.
        M can be either a ModularSymbol_ambient or a dictionary.
        """
        if not isinstance(M,dict):
            raise NotImplementedError("Method needs to be called with dictionary")
        orbits = M['orbits']
        num_orbits = len(orbits)
        level=M.get('level',0); weight=M.get('weight',-1); character = M.get('character',-1)
        level=int(level); weight=int(weight); characer=int(character)
        if ModularSymbols_ambient.query.filter_by(level=level,weight=weight,character=character).count()>0:
            return 
        Md = M['ambient_dict']
        basis = Md['basis']; manin=Md['manin']
        rels = Md['rels']; mod2term=Md['mod2term']
        d_mf = M['ambient'].dimension()
        d_cf = M['ambient'].cuspidal_subspace().dimension()        
        d_new_cf = M['ambient'].new_subspace().cuspidal_subspace().dimension()
        A = ModularSymbols_ambient(level=level,weight=weight,character=character,dimension_modular_forms=d_mf,dimension_cusp_forms=d_cf,dimension_new_cusp_forms=d_new_cf)
        A.set_basis(basis)
        A.set_manin(manin)
        A.set_rels(rels)
        A.set_mod2term(mod2term)
        print "Inserted ambient space"
        for i in range(M.get('num_orbits',0)):
            orbit = orbits[i]
            d = int(orbit.dimension())
            B=str(M['orbits_dict'][i]['B'])
            Bd=str(M['orbits_dict'][i]['Bd'])
            v=unicode(M['orbits_dict'][i]['v'])
            nz=unicode(M['orbits_dict'][i]['nz'])
            print "v=",v
            print "nz=",nz
            Anew = ModularSymbols_newspace_factor(dimension=d)
            Anew.B = B
            Anew.Bd = Bd
            Anew.v = v
            Anew.nz = nz
            A.newspace_factors.append(Anew)
        self._dbb.session.commit()

    def number_of_records(self):
        r"""
        Find how many records are in the new DB.
        """        
        return ModularSymbols_ambient.query.count()            

    def levels_in_DB(self):
        r"""
        Return a list of levels in the database
        """
        level = ModularSymbols_ambient.level
        res = []
        for x in self._dbb.session.query(level).distinct().all():
            res.append(x[0])
        return res
    
    def weights_in_DB(self,level=None):
        weight = ModularSymbols_ambient.weight
        res = []
        query = self._dbb.session.query(weight)
        if level==None:
            query = query.distinct().all()
        else:
            query = query.filter_by(level=int(level))
            query = query.all()
        for x in self._dbb.session.query(weight).distinct().all():
            res.append(x[0])
        return res


    def view_db(self):        
        return ModularSymbols_ambient.query.all()








class WDBtoMongo(WDBtoMFDB):
    r"""

    """
    def __init__(self,datadir,host='localhost',port=37010,verbose=0,**kwds):
        super(WDBtoMongo,self).__init__(datadir,**kwds)
        self._mongo_conn = pymongo.Connection('{0}:{1}'.format(host,port))
        
        self._mongod = pymongo.Connection('{0}:{1}'.format(host,port))['modularforms']
        # Old collection
        self._ms_collection = self._mongod['Modular_symbols']
        

    

    def show_existing_mongo(self):
        files = self._collection['Modular_symbols'].files
        levels = files.distinct('N')
        weights = files.distinct('k')
        print "{0} records with levels in range {1} -- {2}".format(files.count(),min(levels),max(levels))

    def N_k_i_in_files(self,N,k,i):
        s = "N={0} and k={1} and i={2}".format(N,k,i)
        q = self._db.known(s).fetchall()
        # if we have more than one file database
        if self._db2<>None:
            q2 = self._db2.known(s).fetchall()
        q = q + q2
        if len(q)==0:
            return None
        else:
            N,k,i,o,nap=q[0]
            return nap
        return
    
    def convert_mongo_rec_ambient_to_file(self,rec,compute=False):
        fid = rec['_id']
        fs = gridfs.GridFS(self._mongod, 'Modular_symbols')
        f = fs.get(fid)
        ambient = loads(f.read())  # Ambient modular symbols space
        #d = self._db.ambient_to_dict()
        #fname = 'gamma0-ambient-'+self._db.space_name
        N = rec['N']; k=rec['k']; i=rec['chi']        
        filename = self._db.ambient(N, k, i)
        if self.N_k_i_in_files(N,k,i)<>None:
            return 
        self._db.save_ambient_space(ambient,i)
        ## We do not want to compute stuff so we only add orbits
        ## if present in ambient
        if ambient.__dict__.has_key('_HeckeModule_free_module__decomposition') or compute==True:
            self._db.compute_decompositions(N,k,i)
        # if fs.exists({"ambient-filename":fname}):
        #     return False # Record already exists
        # # Update record
        # fs.remove({'_id':fid})
        # N = rec['N']; k=rec['k']; chi=rec['chi']
        # old_filename = rec['filename']
        # dim = int(ambient.dimensionsion())
        # version = str(version())
        # t = (int(N),int(k),int(0))
        #     fs.put(dumps(g[0]),filename=filename,N=int(N), k=int(k),chi=int(0),t=t,version=version)
        #return True


    def convert_records(self,N,k='all'):
        files = self._ms_collection.files
        i=0
        v = [(N,k) for N in rangify(N) for k in rangify(krange)]
        print v
        if k=='all':
            s = {'N':int(N)}
        else:
            s = {'N':N,'k':k}
        print "s=",s
        #print "files=",files
        for f in files.find(s):
            #print "f=",f
            self.convert_mongo_rec_ambient_to_file(f)
        #print "Converted {0} records!"
        return True

    def convert_all_records(self):
        nrange = self._ms_collection.files.distinct('N')
        self.convert_all_records_N(nrange)
        
    @parallel(ncpus=8)
    def convert_all_records_N(self,N):
        if N % 100 == 1:
            print "Converting N={0}".format(N)
        files = self._ms_collection.files        
        for f in files.find({'N':int(N)}):
            self.convert_mongo_rec_ambient_to_file(f)
            
        
    def insert_one_set(self,fs,N,k):
        """
        INPUT:
        N -- positive integer
        k -- even integer >= 2
        OUTPUT:
        (N,k) -- level and weight
        (t0,t1,t2,t3,tall)  -- timings
        Modular symbols space -- with new cuspidal subspace decomposed.
        """
        G = DirichletGroup(N)
        G_orbits = G.galois_orbits(reps_only=True)
        for i in range(len(G_orbits)):
            name = 'data/gamma0-%s-%s-%s.sobj'%(N,k,i)
            chi = G_orbits[i]
            #if chi.is_even() and not is_even(k):
            #   continue
            d=dimension_modular_forms(chi,k)
            if d>200:
                print "Dimension too large! d=",d
                continue
            if d==0:
                continue
            t = (int(N),int(k),int(i))
            if not fs.exists({"t":t}):
                ModularSymbols_clear_cache()
                print name," dim=",d
                #fname = 'gamma0-modsym-%s-%s-%s'%(N,k,i)
                try:
                    #M=load(fname)
                    spaces = self.get_spaces(N,k,i)
                    M = spaces['ambient']
                    M._HeckeModule_free_module__decomposition = {(None, True): Sequence(spaces['orbits'],check=False)}
                except:
                   print "{0},{1},{2}, Not in the old database!".format(N,k,i)
                   continue
                try:
                    try: 
                        fs.put(dumps(M),filename=fname,N=int(N), k=int(k),sign=int(1),chi=i,t=t)
                    except MemoryError:
                        # restart connection and see if this helps...
                        fs = gridfs.GridFS(connection.modularforms,'Modular_symbols') 
                        fs.put(dumps(M),filename=fname,N=int(N), k=int(k),sign=int(1),chi=i,t=t)
                except Exception as e:
                    save(M,fname)
                    print e
                    print "saved to ",fname
                    return
            else:
                print name," is already in db! dim=",d 
            ModularSymbols_clear_cache()
        return N,k

    def get_data(N0,N1,k0=2,k1=12): 
        v = (ellipsis_range(Integer(N0),Ellipsis,Integer(N1)))
        kv = (ellipsis_range(Integer(k0),Ellipsis,Integer(k1+1)))
        step = 1 #len(v)//chunk
        for i in range(step):
            for k in kv: #range(2,13):
                for N in v: #v[chunk*i:chunk*(i+Integer(1))]:
                    fs = gridfs.GridFS(connection.modularforms,'Modular_symbols')	
                    insert_one_set(fs,N,k)




















    
def my_get(dict, key, default, f=None):
    r"""
    Improved version of dict.get where an empty string also gives default.
    and before returning we apply f on the result.
    """
    x = dict.get(key, default)
    if x == '':
        x = default
    if f is not None:
        try:
            x = f(x)
        except:
            pass
    return x

def find_modular_form_data(data={},**kwds):
    level = data.get('level',kwds.get('level',0))
    weight = data.get('weight',kwds.get('weight',0))
    character = data.get('character',kwds.get('character',0))
    q = ModularSymbols_ambient.query.filter(level==level,weight==weight,character==character)
    return q #.all()

## Methods to find or insert objects

def get_input(M=None,**kwds):
    res = {}
    if isinstance(M,sage.modular.modsym.ambient.ModularSymbolsAmbient):
        level = M.level(); weight = M.weight()
        character = int(dirichlet_character_to_int(M.character(), convention='Conrey'))
    elif isinstance(M,db.Model):
        pass
    else:
        res.update(**kwds)
    return res

def create_query(obj,data={},**kwds):
    r"""
    
    """
    q = obj.query
    search = {}
    search.update(data,**kwds)
    #DictLowerJson(data,**kwds)
    keys_num = ['level','weight','character','dimension_modular_forms','number_of_orbits']
    keys_in_db = {key:key for key in keys_num}
    for key in keys_num:
        val = search.get(key,None)        
        if val<>None:
            q=q.filter(obj.__dict__[keys_in_db[key]]==val)
        else: # check range arguments
            val = search.get("{0}_range".format(key),[])
            if len(val)>0:
                q=q.filter(obj.__dict__[key].in_(val))
#    for key in keys_bool:
#        val = search.get(key,None)
#        print "search for {0}=={1}".format(key,val)
#        if val<>None:
#            q=q.filter(obj.__dict__[keys_in_db[key]]==val)
    return q

import sage.modular.modsym.ambient
from conversions import dirichlet_character_to_int,sage_ambient_to_dict,dict_to_factor_sage,factor_to_dict_sage

def ModularSymbols_ambient(M=None,**kwds): #ModularSymbols_ambient_DB, SageObject_DB):
    r"""

    """
    data = get_input(M,**kwds)
    s = create_query(ModularSymbols_ambient_DB_class,data,**kwds)     
    if s.count()>0:
        return s    
    # Else insert a new instance.
    if isinstance(M,sage.modular.modsym.ambient.ModularSymbolsAmbient):
        if M.sign()==0:
            raise ValueError,"Set sign to 1 or -1! otherwise factors wil not be simple!"
        M = sage_ambient_to_dict(M,**kwds) #return ModularSymbols_ambient_from_sage(M,**kwds)
    elif not isinstance(M,dict):
        raise ValueError,"Need to call this with a ModularSymbolsAmbient or dict!!"
    return ModularSymbols_ambient_from_dict(M,**kwds)

## ef ModularSymbols_ambient_from_sage(M,**kwds):
 ##    r"""
 ##    Construct a ModularSymbols_ambient_DB_class instance from sage.
 ##    """
 ##    if not isinstance(M,sage.modular.modsym.ambient.ModularSymbolsAmbient):
 ##        raise ValueError,"Need to call this with a ModularSymbolsAmbient!"
 ##        # From sage object
 ##    new_rec.level = int(M.level())
 ##    new_rec.weight = int(M.weight())
 ##    new_rec.character = int(dirichlet_character_to_int(M.character(), convention='Conrey'))
 ##    new_rec = ModularSymbols_ambient_DB_class(level=level,weight=weight,character=character)
 ##    d=sage_ambient_to_dict(M)
 ##    return update_ModularSymbols_ambient(new_rec,d,**kwds)

def ModularSymbols_ambient_from_dict(data={},**kwds):
    r"""
    Construct a ModularSymbols_ambient_DB_class instance from a dict.
    """
    if not isinstance(data,dict):
        raise ValueError,"Need to call this with a ModularSymbolsAmbient!"
        # From sage object
#    print "data=",data
    new_rec = ModularSymbols_ambient_DB_class()
    new_rec.level=data['space'][0]
    new_rec.weight=data['space'][1]
    new_rec.character=data['space'][2]
    #new_rec.character = int(dirichlet_character_to_int(M.character(), convention='Conrey'))
    #d=sage_ambient_to_dict(M)
    print "updating!"
    return update_ModularSymbols_ambient(new_rec,data,**kwds)

def update_ModularSymbols_ambient(obj,d={},**kwds):
    r"""
    Update properties of a ModularSymbols_ambient_DB_class from dict or by key words.
    """
    basis = d.get('basis',kwds.get('basis',None))
    print "updating 0!"
#    print obj
    if basis <> None:
        obj.basis = basis
    manin = d.get('manin',kwds.get('manin',None))
    if basis <> None:
        obj.manin=d['manin']
    print "updating 1!"
    rels = d.get('rels',kwds.get('rels',None))
    if rels <> None:
        obj.rels = rels
    mod2term = d.get('mod2term',kwds.get('mod2term',None))
    if mod2term <> None:
        obj.mod2term = mod2term
    properties = [ 'dimension_modular_forms','dimension_new_cusp_forms',
                   'dimension_cusp_forms'] 
    for prop in properties:
        val = d.get(prop,kwds.get(prop))
        if val <> None:
            obj.__dict__[prop]=val
    base_field = d.get('base_field',kwds.get('base_field'))
    if base_field == None:
        raise ValuError,"Need complete dict!"
    obj.base_field = ModularSymbols_base_field(base_field)

    for N in d.get('newspace_factors',[]):
        NN = ModularSymbols_newspace_factor(N,ambient=d) #dimension = N['dimension'])
        obj.newspace_factors.append(NN)
        #NN.from_sage(N)
    return obj
            

    # def as_dict(self):
    #     d=dict()
    #     d['basis'] = self.get_basis()
    #     d['manin'] = self.get_manin()
    #     d['rels'] = self.get_rels()
    #     d['mod2term'] = self.get_mod2term()
    #     d['space'] = (self.level, self.weight, self.character)
    #     return d


        

def ModularSymbols_oldspace_factor(*args,**kwds):    
    return NotImplementedError()
    
def ModularSymbols_newspace_factor(N,**kwds):
    if isinstance(N,dict):
        ambient = N.get('ambient',kwds.get('ambient',None))
        if ambient == None:
            raise ValueError,"Need an ambient space."
        N = dict_to_factor_sage(N,ambient=ambient)
    elif not isinstance(N,sage.modular.modsym.subspace.ModularSymbolsSubspace):
        raise ValueError,"Need to be called with dict or ModularSymbolsSubspace instance!"
    return  ModularSymbols_newspace_factor_from_sage(N,**kwds)

def ModularSymbols_newspace_factor_from_sage(N,**kwds):
    names = kwds.get('names','a')
    d=factor_to_dict_sage(N)
    factor = ModularSymbols_newspace_factor_DB_class()
    factor.B  =d['B']
    factor.Bd = d['Bd']
    factor.v =d['v']
    factor.nz = d['nz']
    factor.dimension = int(N.dimension())
    #self.has_cm = has_cm(M)
    # now we set the coefficient field
    extension_field = N.eigenvalue(1,name=names).parent()
    if extension_field != N.base_ring(): # .degree() != 1 and rings.is_NumberField(extension_field):
        assert extension_field.base_field() == N.base_ring()
        minpoly = extension_field.relative_polynomial()
        degree = int(minpoly.degree())
    else:
        # QQ does not have a defining polynomial
        if extension_field.degree()==1:
            minpoly = ZZ['x'].gens()[0]
        else:                
            minpoly = extension_field.defining_polynomial()
        degree = extension_field.degree()
    coefficient_field = CoefficientField(extension_field)
    #coefficient_field.minimal_polynomial = str(minpoly)
    #coefficient_field.base_field = N.ambient().base_field()
    #coefficient_field.degree = degree
    factor.coefficient_field = coefficient_field
    if hasattr(N,'_HeckeModule_free_module__eigenvalues'):
        for n,c in N._HeckeModule_free_module__eigenvalues.iteritems():
            print n,c
            value = str(c[c.keys()[0]].list())
            print value
            a = AlgebraicNumber_DB({'value':value,
                                    'number_field':coefficient_field})
            #if not a in self.coefficient_field.algebraic_numbers:
            #    self.coefficient_field.algebraic_numbers.append(a)
            #cc.value = a
            #cc = Coefficient(index=int(n),value=a)
            factor.coefficients.append( (n,value))

#     def sage_object():
#         return NotImplementedError()

#     def sage_object():
#         return NotImplementedError()

# class CoefficientField(CoefficientField_DB_class, SageObject_DB_class):
#     def sage_object():
#         return NotImplementedError()


