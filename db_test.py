r"""
Inserting records in the database from ?

EXAMPLES::

    sage: 

"""
import inspect,os

basedir =  os.path.dirname(inspect.getabsfile(inspect.currentframe()))
os.sys.path.append("{0}/../mfdb/".format(basedir))
os.sys.path.append("{0}".format(basedir))
os.sys.path.append("{0}/mdb/".format(basedir))

import mfdb
import mdb
import schema
import schema_sage 
from schema_sage import ModularSymbols_ambient,ModularSymbols_newspace_factor,ModularSymbols_oldspace_factor,Coefficient,NumberField,ModularSymbols_base_field, CoefficientField,AlgebraicNumber
#schema.setup_all() 
#schema.create_all()
#DB=mfdb.WDB('git/mfdb/data/')
from mfdb import WDB
#print DB.known(format='web')

class WDBtoMFDB(WDB):
    r"""
    Class to pull records from database in William's format and insert in our database
    """
    def __init__(self,datadir,verbose=0):
        self._dir = datadir
        super(WDBtoMFDB,self).__init__(dir=datadir)
        self._ss = schema_sage
        self._ss.bind = schema.metadata.bind
        if verbose>0:
            self._ss.bind.echo = True
        else:
            self._ss.bind.echo = False
        self._ss.setup_all() 
        self._ss.create_all()       

    def source_db(self):
        r"""
        Return a description of the source database.
        """
        print "Source DB = Files at {0}".format(self._dir)

    def target_db(self):
        r"""
        Return a description of the target database.
        """
        print "Target DB= {0}".format(self._ss.metadata.bind)


        
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
        A = ModularSymbols_ambient(level=level,weight=weight,character=character)
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
            Anew.set_B(B)
            Anew.set_Bd(Bd)
            Anew.set_v(v)
            Anew.set_nz(nz)
            A.newspace_factors.append(Anew)
        self._ss.session.commit()

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
        for x in self._ss.session.query(level).distinct().all():
            res.append(x[0])
        return res
    
    def weights_in_DB(self,level=None):
        weight = ModularSymbols_ambient.weight
        res = []
        query = self._ss.session.query(weight)
        if level==None:
            query = query.distinct().all()
        else:
            query = query.filter_by(level=int(level))
            query = query.all()
        for x in self._ss.session.query(weight).distinct().all():
            res.append(x[0])
        return res


    def view_db(self):        
        return ModularSymbols_ambient.query.all()

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
