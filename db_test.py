r"""
Inserting records in the database from 

"""
import inspect,os

basedir =  os.path.dirname(inspect.getabsfile(inspect.currentframe()))
os.sys.path.append("{0}/../mfdb/".format(basedir))
os.sys.path.append("{0}".format(basedir))
os.sys.path.append("{0}/mdb/".format(basedir))

import mfdb
import mdb
import schema
from schema_sage import ModularSymbols_ambient,ModularSymbols_newspace_factor,ModularSymbols_oldspace_factor,Coefficient,NumberField,ModularSymbols_base_field, CoefficientField,AlgebraicNumber
#schema.setup_all() 
#schema.create_all()
#DB=mfdb.WDB('git/mfdb/data/')

#print DB.known(format='web')

class WDBtoMFDB(WDB):
    r"""
    Class to pull records from database in William's format and insert in our database
    """
    def __init__(datadir):
        super(WDB,self).__init__(datadir)
        schema.setup_all() 
        schema.create_all()       

    def insert_spaces(self,q="N=1 and k=12"):
        r"""
        Insert spaces matching query.
        """
        if q=='all':
            q = ""
        for N,k,ch,numo,nap in self.known(q):
            ## First check if record already exists.
            is_in_db = len(ModularSymbols_ambient.query.filter_by(level=N,weight=k,character=ch))
            print "Inserting {0},{1},{2},{3},{4}".format(N,k,ch,numo,nap)
            d = self.get_spaces(N,k,ch,format='data')[0]
            M = d['ambient']
            d['level']=N; d['weight']=k; d['character']=ch
            orbits = d['orbits']
            assert d['num_orbits']==numo
            d['orbits_dict'] = self.get_decomposition(N,k,ch)[(N,k,ch)]    
            insert_space_into_new_db(d)


    def insert_space_into_new_db(self,M):
        r"""
        Insert M into the database.
        M can be either a ModularSymbol_ambient or a dictionary.
        """
        if not isinstance(M,dict):
            raise NotImplementedError("Method needs to be called with dictionary")
        orbits = M['orbits']
        num_orbits = len(orbits)
        level=M['level']; weight=M['weight']; character = M['character']
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
            A.newform_orbits.append(Anew)
        schema.session.commit()

    def view_db(self):

        return ModularSymbols_ambient.query.all()
