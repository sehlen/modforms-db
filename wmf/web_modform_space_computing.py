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

from sage.all import ZZ, QQ, DirichletGroup, CuspForms, Gamma0, ModularSymbols, Newforms, trivial_character, is_squarefree, divisors, RealField, ComplexField, prime_range, I, join, gcd, Cusp, Infinity, ceil, CyclotomicField, exp, pi, primes_first_n, euler_phi, RR, prime_divisors, Integer, matrix,NumberField,PowerSeriesRing,cached_function,AlphabeticStrings
from sage.rings.power_series_poly import PowerSeries_poly
from sage.all import Parent, SageObject, dimension_new_cusp_forms, vector, dimension_modular_forms, dimension_cusp_forms, EisensteinForms, Matrix, floor, denominator, latex, is_prime, prime_pi, next_prime, previous_prime,primes_first_n, previous_prime, factor, loads,save,dumps,deepcopy
import re
import yaml
from flask import url_for

from wmf import wmf_logger
from lmfdb.modular_forms.elliptic_modular_forms import emf_version

from sage.rings.number_field.number_field_base import NumberField as NumberField_class
from lmfdb.modular_forms.elliptic_modular_forms.backend import connect_to_modularforms_db,get_files_from_gridfs
from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modform_space import WebModFormSpace_class

def WebModFormSpace_computing(N=1, k=2, chi=1, cuspidal=1, prec=10, bitprec=53, data=None, verbose=0,**kwds):
    r"""
    Constructor for WebNewForms with added 'nicer' error message.
    """
    if data is None: data = {}
    if cuspidal <> 1:
        raise IndexError,"We are very sorry. There are only cuspidal spaces currently in the database!"
    try: 
        F = WebModFormSpace_computing_class(N=N, k=k, chi=chi, cuspidal=cuspidal, prec=prec, bitprec=bitprec,
                                        data=data, verbose=verbose,**kwds)
    except ValueError:
        raise IndexError,"We are very sorry. The requested space could not be constructed! Please run the full computations!"
    return F

from lmfdb.modular_forms.elliptic_modular_forms.backend import WebModFormSpace

class WebModFormSpace_computing_class(WebModFormSpace_class):
    r"""
    Space of cuspforms to be presented on the web.
        G  = NS.

    EXAMPLES::

    sage: WS=WebModFormSpace(2,39)


    """
    def __init__(self, N=1, k=2, chi=1, cuspidal=1, prec=10, bitprec=53, data=None, verbose=0,get_from_db=True,**kwds):
        r"""
        Init self.

        INPUT:
        - 'k' -- weight
        - 'N' -- level
        - 'chi' -- character
        - 'cuspidal' -- 1 if space of cuspforms, 0 if all modforms
        """
        wmf_logger.debug("WebModFormSpace with k,N,chi={0}".format( (k,N,chi)))
        super(WebModFormSpace_computing_class,self).__init__(N,k,chi,cuspidal,prec,bitprec,data, verbose,get_from_db=False)

        ## Properties which are here temporarily
        self._hecke_orbits_labels = None
        ## In this subclass we add properties which are not
        ## supposed to be used on the web or stored in the database
        self._dimension = None
        self._dimension_oldspace = None
        self._newforms = None
        self._modular_symbols = None
        self.compute_additional_properties()
        self.insert_into_db()
        self._check_if_all_stored()
        
    def compute_additional_properties(self):
        r"""
        Compute additional properties. 
        """
        ### Set / Compute / fetch everything we need
        if self._group is None:
            self._group = Gamma0(self._N)
        self.get_modular_symbols()
        self._newspace = self.modular_symbols().cuspidal_submodule().new_submodule()
        self.get_character_used_in_computation()
        self.get_galois_orbit_embeddings()
        self.get_newform_factors()
        if self._newforms == {} and self._newspace.dimension()>0:
            for i in self.labels():
                self._newforms[i]=None
        if len(self._ap) == 0:
            self._ap = self._get_aps(prec=self._prec)                
        self.set_dimensions()
        if self.dimension() == self.dimension_newspace():
            self._is_new = True
        else:
            self._is_new = False
        self.set_sturm_bound()
        self.character()
        self.set_oldspace_decomposition()
        self.get_character_galois_orbit()
        self.get_character_orbit_rep()

    def modular_symbols(self):
        r"""
        Return the modular symbols of self.
        """
        if self._modular_symbols is None:
            self.get_modular_symbols()
        return self._modular_symbols
        
    def newform_factors(self):
        r"""
        Return newform factors of self.
        """
        if self._newform_factors is None:
            self._newform_factors = self._get_newform_factors()
        return self._newform_factors
                            
    def get_character_orbit_rep(self,k=None):
        r"""
        Returns canonical representative of the Galois orbit nr. k acting on the ambient space of self.

        """
        from mdb.conversions import dirichlet_character_conrey_galois_orbit_rep
        self._character_orbit_rep = dirichlet_character_conrey_galois_orbit_rep(self.character().character()).number()
    ## Database fetching functions.

     
    def get_character_galois_orbit(self):
        r"""
        Get a list of numbers of the characters in the Galois orbit of the character of self.

        """
        from mdb.conversions import dirichlet_character_conrey_galois_orbit_numbers
        if self._character_galois_orbit is None or self._character_galois_orbit == []:
            if self.level()==1:
                self._character_galois_orbit=[int(1)]
            else:
                self._character_galois_orbit = dirichlet_character_conrey_galois_orbit_numbers(self.level(),self.character().character().number())
                
    def get_character_used_in_computation(self):
        r"""
        Get the character which was used in the computation of the data.
        """
        from mdb.conversions import dirichlet_character_conrey_used_in_computation
        if not self._character_used_in_computation is None:
            return 
        self._character_used_in_computation = dirichlet_character_conrey_used_in_computation(self.level(),self.character().number())

        
    def get_galois_orbit_embeddings(self):
        r"""
        Set the list of Galois orbit embeddings (according to Conrey's naming scheme) of the character of self.
        """
        from mdb.conversions import dirichlet_character_conrey_galois_orbit_embeddings
        self._galois_orbits_embeddings = dirichlet_character_conrey_galois_orbit_embeddings(self.level(),self.character().number())
        

    def _get_aps(self, prec=-1):
        r"""
        Get aps from database if they exist.
        """
        ap_files = connect_to_modularforms_db('ap.files')
        key = {'k': int(self._k), 'N': int(self._N), 'cchi': int(self._chi)}
        key['prec'] = {"$gt": int(prec - 1)}
        ap_from_db  = ap_files.find(key).sort("prec")
        wmf_logger.debug("finds={0}".format(ap_from_db))
        wmf_logger.debug("finds.count()={0}".format(ap_from_db.count()))
        fs = get_files_from_gridfs('ap')
        aplist = {}
        for i in range(len(self.labels())):
            aplist[self.labels()[i]]={}
        for rec in ap_from_db:
            wmf_logger.debug("rec={0}".format(rec))
            ni = rec.get('newform')
            if ni is None:
                for a in self.labels():
                    aplist[a][prec]=None
                return aplist
            a = self.labels()[ni]
            cur_prec = rec['prec']
            if aplist.get(a,{}).get(cur_prec,None) is None:
                aplist[a][prec]=loads(fs.get(rec['_id']).read())
            if cur_prec > prec and prec>0: # We are happy with these coefficients.
                return aplist
        return aplist

    def get_modular_symbols(self):
        r"""
        Get Modular Symbols from database they exist.
        """
        if not self._modular_symbols is None:
            return 
        modular_symbols = connect_to_modularforms_db('Modular_symbols.files')
        key = {'k': int(self._k), 'N': int(self._N), 'cchi': int(self._chi)}
        modular_symbols_from_db  = modular_symbols.find_one(key)
        wmf_logger.debug("found ms={0}".format(modular_symbols_from_db))
        if modular_symbols_from_db is None:
            raise ValueError,"Space is not found in the database!"
            #ms = None
        else:
            id = modular_symbols_from_db['_id']
            fs = get_files_from_gridfs('Modular_symbols')
            ms = loads(fs.get(id).read())
            self._id = id
        self._modular_symbols = ms

  
            
    def get_newform_factors(self):
        r"""
        Get New form factors from database they exist.
        """
        if not self._newforms is None and self._newforms == []:
            return 
        factors = connect_to_modularforms_db('Newform_factors.files')
        key = {'k': int(self._k), 'N': int(self._N), 'cchi': int(self._chi),}
        factors_from_db  = factors.find(key).sort('newform',int(1))
        wmf_logger.debug("found factors={0}".format(factors_from_db.count()))
        self._newforms = {}
        self._hecke_orbits_labels = []
        if factors_from_db.count()==0:
            raise ValueError,"Factors of space {0} is not in database!".format(key)
        else:
            facts = []
            self._labels = []
            fs = get_files_from_gridfs('Newform_factors')
            for rec in factors_from_db:
                factor = loads(fs.get(rec['_id']).read())
                label = orbit_label(rec['newform'])
                self._hecke_orbits_labels.append(label)
                self._newforms[label] = factor
                
 
    def __reduce__(self):
        r"""
        Used for pickling.
        """
        data = self.to_dict()
        return(unpickle_wmfs_v1, (self._k, self._N, self._chi, self._cuspidal, self._prec, self._bitprec, data))
            

    def _repr_(self):
        r"""
        Return string representation of self.
        """
        s = 'Space of Cusp forms on ' + str(self.group()) + ' of weight ' + str(self._k)
        s += ' and dimension ' + str(self.dimension())
        return s


    def _computation_too_hard(self,comp='decomp'):
        r"""
        See if the supplied parameters make computation too hard or if we should try to do it on the fly.
        TODO: Actually check times.
        """
        if comp=='decomp':
            if self._N > 50:
                return True
            if self._chi > 1 and self._N > 100:
                return True
            if self._k+self._N  > 100:
                return True
            return False

    # internal methods to generate properties of self
    def galois_decomposition(self):
        r"""
        We compose the new subspace into galois orbits of new cusp forms.
        """
        from sage.monoids.all import AlphabeticStrings
        if(len(self._galois_decomposition) != 0):
            return self._galois_decomposition
        if '_HeckeModule_free_module__decomposition' in self._newspace.__dict__:
            L = self._newspace.decomposition()
        else:
            decomp = self.newform_factors()
            if len(decomp)>0:
                L = filter(lambda x: x.is_new() and x.is_cuspidal(), decomp)
                wmf_logger.debug("found L:{0}".format(L))
            elif self._computation_too_hard():
                L = []
                raise IndexError,"No decomposition was found in the database!"
                wmf_logger.debug("no decomp in database!")
            else: # compute
                L = self._newspace.decomposition()
                wmf_logger.debug("newspace :".format(self._newspace))                                
                wmf_logger.debug("computed L:".format(L))
        self._galois_decomposition = L
        # we also label the compnents
        x = AlphabeticStrings().gens()
        for j in range(len(L)):
            if(j < 26):
                label = str(x[j]).lower()
            else:
                j1 = j % 26
                j2 = floor(QQ(j) / QQ(26))
                label = str(x[j1]).lower()
                label = label + str(j2)
            if label not in self._heckes_orbits_labels:
                self._hecke_orbits_labels.append(label)
        return L

    def hecke_orbit_label(self, j):
        r"""
        Return the label of the Galois orbit nr. j
        """
        if(len(self._hecke_orbits_labels) == 0):
            self.hecke_decomposition()
        return self._hecke_orbits_labels[j]

    ###  Dimension formulas, calculates dimensions of subspaces of self.
    def set_dimensions(self):
        r"""
        The dimension of the subspace of newforms in self.
        """
        if self._chi != 1:
            x = self.character().sage_character()
        else:
            x = self.level()
        k = self.weight()
        # Ambient modular formsspace
        if self._dimension_modular_forms is None:
            self._dimension_modular_forms = int(dimension_modular_forms(x,k))
        # Cuspidal subspace
        if self._dimension_cusp_forms is None:
            self._dimension_cusp_forms = int(dimension_cusp_forms(x,k))
        # New cuspidal subspace 
        if self._dimension_new_cusp_forms is None:
            self._dimension_new_cusp_forms = int(dimension_new_cusp_forms(x,k))
        # New subspace of ambient space
        if self._dimension_newspace is None:
            if self._cuspidal == 1:
                self._dimension_newspace = self.dimension_new_cusp_forms()
            else:
                self._dimension_newspace = self._newspace.dimension()

        # Old subspace of self.
        if self._dimension_oldspace is None:
            if self._cuspidal == 1:
                self._dimension_oldspace = self.dimension_cusp_forms() - self.dimension_new_cusp_forms()
            else:
                self._dimension_oldspace = self.dimension_modular_forms() - self.dimension_newforms()
                
        if self._dimension is None:
            if self._cuspidal == 1:
                self._dimension = self.dimension_cusp_forms()
            elif self._cuspidal == 0:
                self._dimension = self.dimension_modular_forms()



  
    def set_sturm_bound(self):
        r""" Return the Sturm bound of S_k(N,xi), i.e. the number of coefficients necessary to determine a form uniquely in the space.
        """
        if self._sturm_bound is None:
            self._sturm_bound = self._modular_symbols.sturm_bound()



    def set_oldspace_decomposition(self):
        r"""
        Get decomposition of the oldspace in self into submodules.

        """
        if not (self._oldspace_decomposition is None or self._oldspace_decomposition == []):
            return
        N = self._N
        k = self._k
        M = self._modular_symbols.cuspidal_submodule()
        L = list()
        L = []
        check_dim = self.dimension_newspace()
        if(check_dim == self.dimension()):
            return L
        if(self._verbose > 1):
            wmf_logger.debug("check_dim:={0}".format(check_dim))
        for d in divisors(N):
            if(d == 1):
                continue
            q = ZZ(N).divide_knowing_divisible_by(d)
            if(self._verbose > 1):
                wmf_logger.debug("d={0}".format(d))
            # since there is a bug in the current version of sage
            # we have to try this...
            try:
                O = M.old_submodule(d)
            except AttributeError:
                O = M.zero_submodule()
            Od = O.dimension()
            if(self._verbose > 1):
                wmf_logger.debug("O={0}".format(O))
                wmf_logger.debug("Od={0}".format(Od))
            if(d == N and k == 2 or Od == 0):
                continue
            if self.character().is_trivial():
                # S=ModularSymbols(ZZ(N/d),k,sign=1).cuspidal_submodule().new_submodule(); Sd=S.dimension()
                wmf_logger.debug("q={0},{1}".format(q, type(q)))
                wmf_logger.debug("k={0},{1}".format(k, type(k)))
                Sd = dimension_new_cusp_forms(q, k)
                if(self._verbose > 1):
                    wmf_logger.debug("Sd={0}".format(Sd))
                if Sd > 0:
                    mult = len(divisors(ZZ(d)))
                    check_dim = check_dim + mult * Sd
                    L.append((q, 0, mult, Sd))
            else:
                xd = self.character().character().decomposition()
                for xx in xd:
                    if xx.modulus() == q:
                        Sd = dimension_new_cusp_forms(xx.sage_character(), k)
                        if Sd > 0:
                            # identify this character for internal storage... should be optimized
                            x_k = xx.number()
                            #self.character().character()(xx).number()
                            mult = len(divisors(ZZ(d)))
                            check_dim = check_dim + mult * Sd
                            L.append((q, x_k, mult, Sd))
            if(self._verbose > 1):
                wmf_logger.debug("mult={0},N/d={1},Sd={2}".format(mult, ZZ(N / d), Sd))
                wmf_logger.debug("check_dim={0}".format(check_dim))
        check_dim = check_dim - M.dimension()
        if(check_dim != 0):
            raise ArithmeticError("Something wrong! check_dim=%s" % check_dim)
        self._oldspace_decomposition = L



       
@cached_function
def orbit_label(j):
    x = AlphabeticStrings().gens()
    if(j < 26):
        label = str(x[j]).lower()
    else:
        j1 = j % 26
        j2 = floor(QQ(j) / QQ(26))
        label = str(x[j1]).lower()
        label = label + str(j2)
    return label
