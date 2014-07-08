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
r""" Class for newforms in format which can be presented on the web easily


AUTHORS:

 - Fredrik Stroemberg
 - Stephan Ehlen
"""

import re
import yaml
import pymongo
from flask import url_for


from sage.all import ZZ, QQ, DirichletGroup, CuspForms, Gamma0, ModularSymbols, Newforms, trivial_character, is_squarefree, divisors, RealField, ComplexField, prime_range, I, join, gcd, Cusp, Infinity, ceil, CyclotomicField, exp, pi, primes_first_n, euler_phi, RR, prime_divisors, Integer, matrix,NumberField,PowerSeriesRing,cached_function,AlphabeticStrings,Parent, SageObject, dimension_new_cusp_forms, vector, dimension_modular_forms, dimension_cusp_forms, EisensteinForms, Matrix, floor, denominator, latex, is_prime, prime_pi, next_prime, previous_prime,primes_first_n, previous_prime, factor, loads,save,dumps,deepcopy,sturm_bound
from sage.rings.power_series_poly import PowerSeries_poly
from sage.rings.number_field.number_field_base import NumberField as NumberField_class

from wmf import wmf_logger,WebNewForm_computing
from lmfdb.modular_forms.elliptic_modular_forms import emf_version

from lmfdb.modular_forms.elliptic_modular_forms.backend import connect_to_modularforms_db,get_files_from_gridfs
from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modform_space import WebModFormSpace

from compmf import MongoMF
from lmfdb.modular_forms.elliptic_modular_forms.backend import WebModFormSpace

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
        wmf_logger.debug("WebModFormSpace with k,N,chi={0}".format( (weight,level,character)))      

        super(WebModFormSpace_computing,self).__init__(level,weight,character,cuspidal,prec,bitprec,update_from_db)
        try: 
            self._db = MongoMF(host,port,db)
        except pymongo.errors.ConnectionFailure as e:
            logger.critical("Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message))
            self._db = None
        self.compute_additional_properties()



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
        ### Set / Compute / fetch everything we need
        if self.group is None:
            self.group = Gamma0(self.level)
        self.set_character_used_in_computation()
        self.set_character_galois_orbit()
        self.set_character_orbit_rep()
        self.set_galois_orbit_embeddings()
        self.set_dimensions()
        self.set_sturm_bound()
        self.set_oldspace_decomposition()
        self.get_hecke_orbits()
        self.save_to_db()

    def set_character_galois_orbit(self):
        r"""
        Get a list of numbers of the characters in the Galois orbit of the character of self.

        """
        from compmf.character_conversions import dirichlet_character_conrey_galois_orbit_numbers_from_character_number
        if self._character_galois_orbit is None or self._character_galois_orbit == []:
            if self.level==1:
                self._character_galois_orbit=[int(1)]
            else:
                self._character_galois_orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(self.character.modulus,self.character.number)

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
        self.character_orbit_rep = self.character.character.galois_orbit()[0]

    def set_character_used_in_computation(self):
       r"""
       Get the character which was used in the computation of the data.
       """
       from compmf.character_conversions import dirichlet_character_conrey_used_in_computation
       if not self.character_used_in_computation is None:
           return 
       self.character_used_in_computation = dirichlet_character_conrey_used_in_computation(self.character.modulus,self.character.number)


    def set_dimensions(self):
        r"""
        The dimension of the subspace of newforms in self.
        """
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
        if self.cuspidal == 1:
            self.dimension = self.dimension_cusp_forms
        elif self.cuspidal == 0:
            self.dimension = self.dimension_modular_forms
        # Old subspace of self. We only use the builtin function for cusp forms
        self.dimension_oldforms = self.dimension_cusp_forms - self.dimension_new_cusp_forms

                
    def set_sturm_bound(self):
        r""" Return the Sturm bound of S_k(N,xi), i.e. the number of coefficients necessary to determine a form uniquely in the space.
        """
        if self.sturm_bound is None:
            self.sturm_bound = sturm_bound(self.level,self.weight)


    def get_hecke_orbits(self):
        r"""
        Get the collection of WebNewforms (Hecke orbits) in self. 

        """
        from web_newforms_computing import WebNewForm_computing
        current_dim = 0; i = 0
        while current_dim < self.dimension and i<self.dimension: ## 
            label = orbit_label(i)
            wmf_logger.debug("WebNewForm({0},{1},{2})".format(self.level,self.weight,self.character.number,label))
            F = WebNewForm_computing(self.level,self.weight,self.character.number,label)
            current_dim += F.as_factor().dimension()
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



       
@cached_function
def orbit_label(j):
    x = AlphabeticStrings().gens()
    j1 = j % 26
    label = str(x[j1]).lower()
    j  = (j - j1) / 26 -1
    while j >= 0:
        j1 = j % 26
        label = str(x[j1]).lower() + label
        j = (j - j1) / 26 - 1
    return label


@cached_function
def orbit_index_from_label(label):
    r"""
    Inverse of the above
    """
    res = 0
    A = AlphabeticStrings()
    x = A.gens()
    l = list(label)
    su = A(l.pop().upper())
    res = x.index(su)
    l.reverse()
    i = 1
    for s in l:
        su = A(s.upper())
        res+=(1+x.index(su))*26**i
        i+=1
    return res
