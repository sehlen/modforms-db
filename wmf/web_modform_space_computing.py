# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <fredrik314@gmail.com>,
#  Stephan Ehlen <>
#  Distributed under the terms of the GNU General Public License (GPL)
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
r""" Class for spaces of modular forms which can be presented on the web easily


AUTHORS:

 - Fredrik Stroemberg
 - Stephan Ehlen
"""

import re
import yaml
import pymongo
from flask import url_for
import datetime

from sage.all import ZZ, QQ, DirichletGroup, CuspForms, Gamma0, ModularSymbols, Newforms, trivial_character, is_squarefree, divisors, RealField, ComplexField, prime_range, I,gcd, Cusp, Infinity, ceil, CyclotomicField, exp, pi, primes_first_n, euler_phi, RR, prime_divisors, Integer, matrix,NumberField,PowerSeriesRing,cached_function,AlphabeticStrings,Parent, SageObject, dimension_new_cusp_forms, vector, dimension_modular_forms, dimension_cusp_forms, EisensteinForms, Matrix, floor, denominator, latex, is_prime, prime_pi, next_prime, previous_prime,primes_first_n, previous_prime, factor, loads,save,dumps,deepcopy,sturm_bound
from sage.rings.power_series_poly import PowerSeries_poly
from sage.rings.number_field.number_field_base import NumberField as NumberField_class

from lmfdb.modular_forms.elliptic_modular_forms import emf_version
from lmfdb.modular_forms.elliptic_modular_forms.backend import connect_to_modularforms_db,get_files_from_gridfs
from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modform_space import WebModFormSpace
from lmfdb.modular_forms.elliptic_modular_forms.backend import WebModFormSpace

from compmf import MongoMF

from wmf import wmf_logger,WebNewForm_computing,orbit_index_from_label,orbit_label


class WebModFormSpace_computing(WebModFormSpace):
    r"""
    Class for computing properties of spaces of modular forms which will be presented on the web.

    EXAMPLES::

    sage: WS=WebModFormSpace(2,39)


    """
    def __init__(self, level=1, weight=12, character=1, cuspidal=1, prec=10, bitprec=53, update_from_db=True,host='localhost',port=37010,db='modularforms2',**kwds):
        r"""
        Init self.

        INPUT:

        
        - 'k' -- weight
        - 'N' -- level
        - 'chi' -- character
        - 'cuspidal' -- 1 if space of cuspforms, 0 if all modforms

        - 'prec' -- integer (default 10)
        - 'bitprec' -- 53
        - 'get_from_db' -- Boolean (default True)
        - 'host' -- string (default 'localhost')
        - 'port' -- integer (default 37010)
        - 'db' -- string (default modularforms2)

        
        """
        if isinstance(level,str):  ## It is probable a label
            level,weight,character = map(int,level.split("."))
        wmf_logger.debug("WebModFormSpace with k,N,chi={0}".format( (weight,level,character)))
        self._host = host; self._port=int(port); self._dbname = db
        
        super(WebModFormSpace_computing,self).__init__(level=level,weight=weight,character=character,cuspidal=cuspidal,prec=prec,bitprec=bitprec,update_from_db=update_from_db,**kwds)
        wmf_logger.debug("Super class is inited! dim of self={0}".format(self.dimension))
        self._rec = {}
        #if self.dimension == 0: ### does not work if the record has
        #not been computed...
        #    return 
        self.setup_modular_symbols_db()
        if kwds.get('recompute',True):
            self.compute_additional_properties()
        else:
            self.update_dimension_table()
    def setup_modular_symbols_db(self):
        r"""
        Connect to the mongodb with modular symbols and fetch the current record.
        """
        try: 
            self._db = MongoMF(self._host,self._port,self._dbname)
            # find the record in the database which corresponds to self
            s = {'N':int(self.level),'k':int(self.weight),
                 'character_galois_orbit':{"$in":[int(self.character.number)]}}
            self._rec = self._db._modular_symbols.find_one(s)
            if self._rec is None:
                wmf_logger.critical("Could not find the space {0} in the database! This should be computed first".format((self.level,self.weight,self.character)))
                self._rec={}
                return
        except pymongo.errors.ConnectionFailure as e:
            wmf_logger.critical("Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message))
            self._db = None  
            self._rec = {}
            
    def _repr_(self):
        r"""
        Return string representation of self.
        """
        s = 'Space of Cusp forms on ' + str(self.group()) + ' of weight ' + str(self._k)
        s += ' and dimension ' + str(self.dimension())
        return s

        
    def compute_additional_properties(self):
        r"""
        Compute additional properties. 
        """
        wmf_logger.debug("Computing addiitonal properties!")
        ### Set / Compute / fetch everything we need
        if self.group is None:
            self.group = Gamma0(self.level)
        self.set_dimensions()
        self.set_character_used_in_computation()
        self.set_character_galois_orbit()
        self.set_character_orbit_rep()
        self.set_galois_orbit_embeddings()
        wmf_logger.debug("Got all characters!")
        if self.dimension == 0:
            self.save_to_db()
            return 
        self.set_sturm_bound()
        wmf_logger.debug("Got sturm bound!")
        self.set_oldspace_decomposition()
        wmf_logger.debug("Got oldspace decomposition!")
        self.get_hecke_orbits()
        self.get_zetas()
        self.version = float(emf_version)
        self.creation_date=datetime.datetime.utcnow()
        self.save_to_db()
        self.update_dimension_table()


    def get_zetas(self):
        r"""
        Make a dictionary of all the zetas which appear in the q-expansions of the modular forms on self.
        """
        from sage.all import Infinity
        self.zeta_orders = []
        for a in self.hecke_orbits:
            f = self.hecke_orbits[a]
            K = f.q_expansion.base_ring()
            if K == QQ:
                continue
            # We have to distinguish between the cases: K = Cyclotomic, K= extension of cyclotomic and K is Q and K is non-cyclotomic extension of Q
            if K.is_relative():
                n = K.base_field().gen().multiplicative_order()
            else:
                n = K.gen().multiplicative_order()
            if n==+Infinity:
                wmf_logger.debug("Got z of infinite order! K={0}".format(K))
            elif not n in self.zeta_orders:
                self.zeta_orders.append(int(n))
                
    def set_character_galois_orbit(self):
        r"""
        Get a list of numbers of the characters in the Galois orbit of the character of self.

        """
        self._character_galois_orbit = self._rec.get('character_galois_orbit',[])
        

    def set_galois_orbit_embeddings(self):
        r"""
        Set the list of Galois orbit embeddings (according to Conrey's naming scheme) of the character of self.
        """
        from compmf.character_conversions import dirichlet_character_conrey_galois_orbit_embeddings
        self._galois_orbits_embeddings = dirichlet_character_conrey_galois_orbit_embeddings(self.character.modulus,self.character.number)
        
    def set_character_orbit_rep(self):
        r"""
        Returns canonical representative of the Galois orbit of the character of self.

        """        
        from compmf.character_conversions import conrey_character_from_number,conrey_character_number_to_conrey_galois_orbit_number
        if self._character_galois_orbit <> []:
            ci = min(self._character_galois_orbit)
        else:
            ci = 1
        on = conrey_character_number_to_conrey_galois_orbit_number(self.level,ci)[1]

        self.character_orbit_rep = conrey_character_from_number(self.level,ci)
        #self.galois_orbit_name = "{0}.{1}.{2}".format(self.level,self.weight,ci)
        self.space_orbit_label = "{0}.{1}.{2}".format(self.level,self.weight,on)
        
    def set_character_used_in_computation(self):
       r"""
       Get the character which was used in the computation of the data.
       NOTE: The character indicated by the space_label SHOULD be the character we used for the data. 
       """
       from compmf.character_conversions import conrey_character_from_number
       self.character_used_in_computation = conrey_character_from_number(self.character.modulus,self.character.number)


    def set_dimensions(self):
        r"""
        The dimension of the subspace of newforms in self. This should all be in the database.
        """

        self.dimension_modular_forms = self._rec.get('dima',-1)
        self.dimension_cusp_forms = self._rec.get('dimc',-1)
        self.dimension_new_cusp_forms = self._rec.get('dimn',-1)
        if self.dimension_modular_forms < 0 or self.dimension_cusp_forms < 0 or self.dimension_cusp_forms < 0:
            ## If we didn't get dimensions we compute them.
            if self.character.number != 1 and self.level<>1:
                x = self.character.sage_character
            else:
                x = self.level
            k = self.weight
            # Ambient modular formsspace
            self.dimension_modular_forms = int(dimension_modular_forms(x,k))
            # Cuspidal subspace
            self.dimension_cusp_forms = int(dimension_cusp_forms(x,k))
            # New cuspidal subspace 
            self.dimension_new_cusp_forms = int(dimension_new_cusp_forms(x,k))
            self.dimension_eisenstien = self.dimension_modular_forms - self.dimension_cusp_forms
        if self.cuspidal == 0:
            self.dimension = self.dimension_modular_forms
        else:
            self.dimension = self.dimension_cusp_forms
        # Old subspace of self. We only use the builtin function for cusp forms
        self.dimension_oldforms = self.dimension_cusp_forms - self.dimension_new_cusp_forms
        wmf_logger.debug("dim_mod={0}".format(self.dimension_modular_forms))
        wmf_logger.debug("dim_cusp={0}".format(self.dimension_cusp_forms))
        wmf_logger.debug("dim_new_cusp={0}".format(self.dimension_new_cusp_forms))        
        

    def update_dimension_table(self):
        if self._db is None:
            self.setup_modular_symbols_db()
        if self.version > 1.3:
            C = self._db._mongodb['dimension_table2']
        else:
            C = self._db._mongodb['dimension_table']
        in_db = True
        for x in self.hecke_orbits:
            label = self.hecke_orbits[x].hecke_orbit_label
            wmf_logger.debug("labe={0}".format(label))
            if list(WebNewForm_computing.find({'hecke_orbit_label':label,'version':self.version}))==[]:
                wmf_logger.critical("The hecke orbit {0} is not in the database".format(label))
                in_db = False
                break
        r = C.find_one({'space_label':self.space_label})
        if not r is None:
            C.update({'_id':r['_id']},{"$set":{'in_wdb':int(in_db),'in_msdb':int(1)}})
        else:
            r = {'space_orbit_label':self.space_orbit_label,
                 'space_label':self.space_label,
                 'character_orbit':map(int,self.character.character.galois_orbit()),
                 'level':int(self.level),
                 'weight':int(self.weight),
                 'cchi':int(self.character.number),
                 'd_mod':int(self.dimension_modular_forms),
                 'd_cusp':int(self.dimension_cusp_forms),
                 'd_newf':int(self.dimension_new_cusp_forms),
                 'd_eis':int(self.dimension_modular_forms - self.dimension_cusp_forms),
                 'in_wdb':int(in_db),
                 'in_msdb':int(1)}
            C.insert(r)
        # Check Gamma1 as well...??
        #norbits = M.character.character.parent()._galois_orbits()
        
    def set_sturm_bound(self):
        r""" Return the Sturm bound of S_k(N,xi), i.e. the number of coefficients necessary to determine a form uniquely in the space.
        """
        if self.sturm_bound is None or self.sturm_bound == 0:
            self.sturm_bound = sturm_bound(self.level,self.weight)


    def get_hecke_orbits(self):
        r"""
        Get the collection of WebNewforms (Hecke orbits) in self. 

        """
        from web_newforms_computing import WebNewForm_computing
        current_dim = 0; i = 0
        dim = self.dimension_new_cusp_forms
        wmf_logger.debug("Dimension={0}".format(dim))
        while current_dim < dim and i<dim: ## 
            label = orbit_label(i)
            wmf_logger.debug("WebNewForm({0},{1},{2},{3})".format(self.level,self.weight,self.character.number,label))
            F = WebNewForm_computing(self.level,self.weight,self.character.number,label,recompute=True)
            current_dim += F.dimension
            i+=1
            wmf_logger.debug("current_dim={0}".format(current_dim))
            self.hecke_orbits[label]=F
                       
    def set_oldspace_decomposition(self):
        r"""
        Get decomposition of the oldspace in S_k(N) into submodules.

        """
        if not self.cuspidal == 1:
            return 
        if not (self.oldspace_decomposition is None or self.oldspace_decomposition == []):
            return
        old_dim = self.dimension_cusp_forms - self.dimension_new_cusp_forms
        if old_dim == 0:
            return 
        wmf_logger.debug("oldspace  dimension:={0}".format(old_dim))
        N = self.level; k = self.weight
        L = []
        check_dim = 0  
        for d in divisors(N):
            if(d == 1):
                continue
            q = ZZ(N).divide_knowing_divisible_by(d)
            wmf_logger.debug("d={0} q = {1}".format(d,q))
            if self.character.is_trivial():
                Sd = dimension_new_cusp_forms(q, k)
                wmf_logger.debug("Sd={0}".format(Sd))
                if Sd > 0:
                    mult = len(divisors(ZZ(d)))
                    check_dim = check_dim + mult * Sd
                    L.append((q, 0, mult, Sd))
            else:
                xd = filter(lambda x: x.modulus() == q,self.character.character.decomposition())
                for xx in xd:
                    Sd = dimension_new_cusp_forms(xx.sage_character(), k)
                    if Sd > 0:
                        mult = len(divisors(ZZ(d)))
                        check_dim = check_dim + mult * Sd
                        L.append((q, xx.number(), mult, Sd))
                wmf_logger.debug("mult={0},N/d={1},Sd={2}".format(mult, ZZ(N / d), Sd))
                wmf_logger.debug("check_dim={0}".format(check_dim))
            if check_dim == old_dim:
                break
        if check_dim <> old_dim:
            raise ArithmeticError("Something wrong! check_dim=%s" % check_dim)
        self.oldspace_decomposition = L


