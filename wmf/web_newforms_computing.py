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


NOTE: We are now working completely with the Conrey naming scheme.
 
TODO:
Fix complex characters. I.e. embedddings and galois conjugates in a consistent way.

"""
from sage.all import ZZ, QQ, DirichletGroup, CuspForms, Gamma0, ModularSymbols, Newforms, trivial_character, is_squarefree, divisors, RealField, ComplexField, prime_range, I, join, gcd, Cusp, Infinity, ceil, CyclotomicField, exp, pi, primes_first_n, euler_phi, RR, prime_divisors, Integer, matrix,NumberField,PowerSeriesRing,cached_function,PolynomialRing
from sage.rings.power_series_poly import PowerSeries_poly
from sage.all import Parent, SageObject, dimension_new_cusp_forms, vector, dimension_modular_forms, dimension_cusp_forms, EisensteinForms, Matrix, floor, denominator, latex, is_prime, prime_pi, next_prime, previous_prime,primes_first_n, previous_prime, factor, loads,save,dumps,deepcopy,sturm_bound
import re
import yaml
from flask import url_for
import pymongo
## DB modules


from lmfdb.modular_forms.elliptic_modular_forms import emf_logger,emf_version

from sage.rings.number_field.number_field_base import NumberField as NumberField_class
from lmfdb.modular_forms.elliptic_modular_forms.backend import connect_to_modularforms_db,get_files_from_gridfs
from lmfdb.modular_forms.elliptic_modular_forms.backend.web_newforms import WebNewForm,WebEigenvalues


from lmfdb.modular_forms.elliptic_modular_forms.backend.emf_utils import newform_label, space_label

from wmf import wmf_logger
from compmf import MongoMF


class WebNewForm_computing(WebNewForm):
    r"""
    Class for representing a (cuspidal) newform on the web.
    TODO: Include the computed data in the original database so we won't have to compute here at all.
    """
    def __init__(self,level=1, weight=12, character=1, label='a', prec=10, bitprec=53, parent=None,host='localhost',port=37010,db='modularforms2',recompute=False):
        r"""
        Init self as form with given label in S_k(N,chi)

        INPUT:

        
        
        """
        
        super(WebNewForm_computing,self).__init__(level,weight,character,label,prec,bitprec,parent)
        self.hecke_orbit_label = newform_label(self.level,self.weight,self.character.number,self.label)
        try: 
            self._db = MongoMF(host,port,db)
        except pymongo.errors.ConnectionFailure as e:
            logger.critical("Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message))
            self._db = None
        wmf_logger.debug("WebNewForm_computing with N,k,chi,label={0}".format( (self.level,self.weight,self.character,self.label)))
        self._as_factor = None
        self._prec_needed_for_lfunctions = None
        self._available_precisions = []
        self._newform_number = None
        self._satake = {}
        self._twist_info = None
        self._as_polynomial_in_E4_and_E6 = None
        ## If it is in the database we don't need to compute everything unless we specify recompute=True
        if self._db._mongodb[self._collection_name].find({'hecke_orbit_label':self.hecke_orbit_label}).count()>0 and recompute is False:
            wmf_logger.debug("getting data from db")
            self.update_from_db()
            self.compute_satake_parameters_numeric()
            self.set_twist_info()
        else:
            self.compute_additional_properties()
        self.save_to_db()

    def __repr__(self):
        r"""
        """
        s = "WebNewform for computing in S_{0}({1},chi_{2}) with label {3}".format(self.weight,self.level,self.character.number,self.label)
        return s
        
        
    def compute_additional_properties(self):
        r"""
        Compute everything we need.
        """
        wmf_logger.debug("Update ap's")

        self.set_dimension()
        if self.dimension == 0:
            return 
        self.set_aps()
        self.set_q_expansion()
        self.set_q_expansion_embeddings()
        self.set_base_ring()       
        self.set_coefficient_field()

        self.set_twist_info()
        self.set_is_cm()
        self.compute_satake_parameters_numeric()

        self.set_atkin_lehner()
        self.set_absolute_polynomial()
        if self.level==1:
            self.explicit_formulas['as_polynomial_in_E4_and_E6'] = self.as_polynomial_in_E4_and_E6()

##  Internal functions
##

    def set_dimension(self):
        r"""
        The dimension of this galois orbit is not necessarily equal to the degree of the number field, when we have a character....
        We therefore need this routine to distinguish between the two cases...
        """
        if not self._properties['dimension'].has_been_set():
            try:
                self.dimension = self.as_factor().dimension()
            except ValueError:
                # check if we have a zero space
                if self.character.number == 1:
                    dim = dimension_new_cusp_forms(self.level,self.weight)
                else:
                    dim = dimension_new_cusp_forms(self.character.sage_character,self.weight)
                if dim == 0:
                    self.dimension = int(0)
                else:
                    raise ValueError,"Ambient space is not zero dimensional so this function {0} is just not computed!".format(self.hecke_orbit_label)
    def set_aps(self,reload_from_db=False):
        r"""
        We set the eigenvalues unless we already have sufficiently many and reload_from_db is False
        
        """
        wmf_logger.debug("Setting aps!")
        if self.eigenvalues.prec >= self.prec_needed_for_lfunctions() and reload_from_db is False:
            return 
        #try:
        aps = self._db.get_aps(self.level,self.weight,self.character.number,self.newform_number(),character_naming='conrey')
        wmf_logger.debug("Got ap lists:{0}".format(len(aps)))
        wmf_logger.debug("Want:{0} coefficients!".format(self.prec_needed_for_lfunctions()))
        ev_set = 0
        precs = []
        if aps<>{}:
            if isinstance(aps.values()[0],dict):
                precs = aps.values()[0].keys()
        wmf_logger.critical("precs={0}".format(precs))
        for prec in precs:
            wmf_logger.debug("Now getting prec@{0}".format(prec))
            E,v,meta = aps.values()[0][prec]
            self._available_precisions.append(prec)
            evs = WebEigenvalues(self.hecke_orbit_label,prec)
            evs.E = E
            wmf_logger.critical("E = {0}".format(E))
            evs.v = v
            evs.meta = meta
            evs.init_dynamic_properties()
            t = evs.save_to_db(update=True)
            if not t is True:
                wmf_logger.critical("Could not update webeigenvalues!")
            if prec >= self.prec and ev_set == 0:
                self.eigenvalues = evs
                ev_set = 1
            wmf_logger.debug("Got ap's with prec={0}".format(self.eigenvalues.prec))
            if self.eigenvalues.prec >= self.prec_needed_for_lfunctions():
                ## We got as many as we wanted.
                break
        if self.eigenvalues.prec < self.prec_needed_for_lfunctions():
            wmf_logger.critical("Could not find coefficients with prec:{0} Only got:{1}".format(self.prec_needed_for_lfunctions(),self.eigenvalues.prec))
        #except Exception as e:
        #wmf_logger.critical("Could not get ap's. Error:{0}".format(e.message))

    def set_q_expansion(self):
        r"""
        Set the q-expansion for self.
        """
        wmf_logger.debug("Set q-expansion")
        if not self.eigenvalues.has_eigenvalue(2):
            self.set_aps()
        QR = PowerSeriesRing(self.base_ring,name='q',order='neglex')
        q = QR.gen()
        res = 0
        m = max(self.prec,self.prec_needed_for_lfunctions())
        self.coefficients(range(1,m))
        for n in range(1,m):
            res+=self.coefficient(n)*q**n
        self.q_expansion = res


    def as_factor(self):
        r"""
        Return self as a newform factor
        """
        if self._as_factor is None:
            C = connect_to_modularforms_db('Newform_factors.files')
            s = {'N':int(self.level),'k':int(self.weight),
                 'cchi':int(self.character.number),'newform':int(self.newform_number())}
            res = C.find_one(s)
            #print res
            if not res is None:
                fid = res['_id']
                fs = get_files_from_gridfs('Newform_factors')
                self._as_factor = loads(fs.get(fid).read())
        if self._as_factor is None:
            raise ValueError,"Newform matching {0} can not be found in the database!".format(s)
        return self._as_factor

    def set_base_ring(self):
        r"""
        The base ring of self, that is, the field of values of the character of self. 
        """
        if self.base_ring is None or self.base_ring=={}:
            self.base_ring = self.as_factor().base_ring()
        if self.base_ring == QQ:
            self.is_rational = True
        
    def set_coefficient_field(self):
        r"""
        The base ring of self, that is, the field of values of the character of self. 
        """
        if not self.eigenvalues.has_eigenvalue(2):
            self.set_aps()
        try:
            self.coefficient_field = self.eigenvalues[2].parent()
        except KeyError:
            raise KeyError,"We do not have eigenvalue a(2) for this newform!"

    def set_absolute_polynomial(self):
        r"""
        Set the absolute polynomial of self.
        """
        if self.coefficient_field == QQ:
            self.absolute_polynomial = PolynomialRing(QQ,'x').gen()
        else:
            self.absolute_polynomial = self.coefficient_field.absolute_polynomial()
        

    def newform_number(self):
        r"""
        Return the index of self in the list of all Hecke orbits.
        """
        from wmf.web_modform_space_computing import orbit_index_from_label
        if self._newform_number is None:
            #print "label=",self.label,type(self.label)
            self._newform_number = orbit_index_from_label(self.label)
        return self._newform_number

    def get_character_orbit_rep(self):
        r"""
        Get the representative of the Galois orbit of the character of self.
        """
        from compmf.character_conversions import dirichlet_character_conrey_galois_orbit_rep_from_character_number
        self._character_orbit_rep = dirichlet_character_conrey_galois_orbit_rep_from_ncharacter_umber(self.character.modulus,self.character.number)
        
    def prec_needed_for_lfunctions(self):
        r"""
        Calculates the number of coefficients needed for the L-function
        main page (formula taken from their pages)
        """
        self._prec_needed_for_lfunctions = int(22 + int(RR(5)*RR(self.weight)*RR(self.level).sqrt()) + 1)
        return self._prec_needed_for_lfunctions 
    
    def set_q_expansion_embeddings(self, prec=-1, bitprec=53,format='numeric',display_bprec=26):
        r""" Compute all embeddings of self into C which are in the same space as self.
        Return 0 if we didn't compute anything new, otherwise return 1.
        """
        if prec <= 0:
            prec = self.prec_needed_for_lfunctions()
        wmf_logger.debug("computing embeddings of q-expansions : has {0} embedded coeffs. Want : {1} with bitprec={2}".format(len(self._embeddings),prec,bitprec))
        ## First check if we have sufficient data
        if self._embeddings.get('prec',0) >= prec and self._embeddings.get('bitprec',0) >= bitprec:
            return 0 ## We should already have sufficient data.
        ## Else we compute new embeddings.
        CF = ComplexField(bitprec)
        # First wee if we need higher precision, in which case we reset all coefficients:
        if self._embeddings.get('bitprec',0) < bitprec:
            self._embeddings['values']={}
            self._embeddings['prec']=int(0)
            self._embeddings['bitprec']=int(bitprec)
        # See if we have need of more coefficients
        nstart = len(self._embeddings['values'])
        wmf_logger.debug("Should have {0} embeddings".format(self._embeddings['prec']))
        wmf_logger.debug("Computing new embeddings !")
        deg = self.coefficient_field.absolute_degree()
        for n in range(self._embeddings['prec'],prec+1):
            try:
                cn = self.coefficient(n)
            except IndexError:
                break
            if hasattr(cn, 'complex_embeddings'):
                embc = cn.complex_embeddings(bitprec)
            else:
                embc = [ CF(cn) ] # this is only occuring for Q
            self._embeddings['values'][n]=embc
        self._embeddings['prec'] = prec+1
        return 1
                      
    
  
    def set_atkin_lehner(self):
        r"""
        Get the Atkin-Lehner eigenvalues from database if they exist. 
        """
        
        if not ((self.character.is_trivial() or self.character.order == 2) and not self._atkin_lehner_eigenvalues is None):
            return None
        C = connect_to_modularforms_db('Atkin_Lehner.files')
        fs = get_files_from_gridfs('Atkin_Lehner')
        s = {'N':int(self.level),'k':int(self.weight),
             'cchi':int(self.character.number),'newform':int(self.newform_number())}
        res = C.find_one(s)
        self._atkin_lehner_eigenvalues = {}
        if not res is None:
            alid = res['_id']
            al = loads(fs.get(alid).read())
            wmf_logger.debug("al = {0}".format(al))
            i = 0 
            for d in prime_divisors(self.level):
                wmf_logger.debug("d = {0}".format(d))
                self._atkin_lehner_eigenvalues[d] = 1 if al[i]=='+' else -1
                i+=1
        else: # We compute them
            A = self.as_factor()
            for p in prime_divisors(self.level):
                if self.character.is_trivial() or p==self.level:
                    self._atkin_lehner_eigenvalues[p]= int(A.atkin_lehner_operator(p).matrix()[0,0])

    def set_twist_info(self, prec=10,insert_in_db=True):
        r"""

        Try to find forms of lower level which get twisted into self.

        OUTPUT:

        -''[t,l]'' -- tuple of a Bool t and a list l. The list l contains all tuples of forms which twists to the given form.
        The actual minimal one is the first element of this list.
             t is set to True if self is minimal and False otherwise


        EXAMPLES::



        """
        from compmf.character_conversions import dirichlet_character_conrey_galois_orbits_reps,conrey_character_from_number
        from wmf.web_modform_space_computing import WebModFormSpace_computing
        N = self.level
        k = self.weight
        if(is_squarefree(ZZ(N))):
            self._twist_info = [True, None ]
            return [True, None]

        # We need to check all square factors of N
        twist_candidates = list()
        KF = self.base_ring
        # We need coefficients up to the Sturm bound so the primes we need are :
        possible_spaces = []
        ## We first find all possible parameters for spaces which can twist into self.
        character_decomp = {}
        for xx in self.character.character.decomposition():
            if not character_decomp.has_key(xx.conductor()):
                character_decomp[xx.conductor()]=xx
            else:
                character_decomp[xx.conductor()]  = xx*character_decomp[xx.conductor()]
        wmf_logger.debug("decomposition of character of self:{0}".format(character_decomp))
        for d in divisors(N):  ## We twist with character mod d
            if d == 1:
                continue
            # we look at all d such that d^2 divdes N
            wmf_logger.debug("Check character d={0}".format(d))            
            if not ZZ(d**2).divides(ZZ(N)):
                continue
            D = dirichlet_character_conrey_galois_orbits_reps(d)
            ## We know that N = lcm(M,d^2)
            M = ZZ(N).divide_knowing_divisible_by(d ** 2)
            for MM in M*d.divisors():
                wmf_logger.debug("Check space M={0}".format(MM))
                ## We have to check if self.character  = x*y^2
                ## where x is a character mod M and y mod d
                DM = dirichlet_character_conrey_galois_orbits_reps(MM)  ## These are the spaces which might be in the database
                xt = conrey_character_from_number(N,1) ## the trivial character mod N
                for x in DM:
                    for y in D:
                        ## Get the decomposition of left hand side (since multiplication of character with different modulus is not implemented in DirichletGroup_conrey)
                        lhs_decomp = {}
                        for xx in x.decomposition():
                            if not lhs_decomp.has_key(xx.conductor()):
                                lhs_decomp[xx.conductor()]=xx
                            else:
                                lhs_decomp[xx.conductor()]  = xx*lhs_decomp[xx.conductor()]
                        for yy in y.decomposition():
                            if not lhs_decomp.has_key(yy.conductor()):
                                lhs_decomp[yy.conductor()]=yy
                            else:
                                lhs_decomp[yy.conductor()]  = yy*lhs_decomp[yy.conductor()]
                        wmf_logger.debug("lhs_Decomp = {0}".format(lhs_decomp))
                        try:
                            for m in character_decomp.keys():
                                if m not in lhs_decomp.keys():
                                    continue
                                wmf_logger.debug("lhs={0}".format(lhs_decomp[m].primitive_character()))
                                wmf_logger.debug("rhs={0}".format(character_decomp[m].primitive_character()))
                                if lhs_decomp[m].primitive_character()<>character_decomp[m].primitive_character():
                                    wmf_logger.debug("characters are unequal!")
                                    raise StopIteration
                            wmf_logger.debug("Characters are equal")
                            t = (MM,x.number(),d,y.number())
                            if t not in possible_spaces:
                                possible_spaces.append(t)
                        except StopIteration:
                            pass
                    
            
        wmf_logger.debug("Possible spaces: {0}".format(possible_spaces))
        coeffs = self.coefficients(range(self.parent.sturm_bound))
        K = self.coefficient_field
        for t in possible_spaces:
            M,xi,d,yi = t
            wmf_logger.debug("Checking level {0}, CHARACTER {1}".format(M,xi))
            MF = WebModFormSpace_computing(M,k,xi)
            if MF.dimension == 0:
                continue
            y = conrey_character_from_number(d,yi).sage_character()
            for label in MF.hecke_orbits:
                F = MF.hecke_orbits[label]
                wmf_logger.debug("Checking function F={0}".format(F))
                coeffsF = F.coefficients(range(self.parent.sturm_bound))
                coeffsF_twist = []
                for j in range(self.parent.sturm_bound):
                    coeffsF_twist.append(K(coeffs[j])*K(y(j)))
                wmf_logger.debug("own coeffs={0}".format(coeffs))
                wmf_logger.debug("twisted coeffs={0}".format(coeffsF_twist))
                try:
                    for j in range(self.parent.sturm_bound):
                        if coeffsF_twist[j]<>coeffs[j]:
                            raise StopIteration
                    if F.hecke_orbit_label not in twist_candidates:
                        twist_candidates.append(F.hecke_orbit_label)
                except StopIteration: # If we are here then we don't have a twist to self
                    pass

        wmf_logger.debug("Candidates=v{0}".format(twist_candidates))
        self._twist_info = (False, twist_candidates)
        if len(twist_candidates) == 0:
            self._twist_info = [True, None]
        else:
            self._twist_info = [False, twist_candidates]
        self.twist_info = self._twist_info
        return self._twist_info

    def set_is_cm(self,insert_in_db=True):
        r"""
        Checks if f has complex multiplication and if it has then it returns the character.

        OUTPUT:

        -''[t,x]'' -- string saying whether f is CM or not and if it is, the corresponding character

        EXAMPLES::

        """
        #if(len(self._is_CM) > 0):
        #    return self._is_CM
        max_nump = self._number_of_hecke_eigenvalues_to_check()
        # E,v = self._f.compact_system_of_eigenvalues(max_nump+1)
        try:
            coeffs = self.coefficients(range(max_nump + 1))
        except IndexError: 
           return None,None
        nz = coeffs.count(0)  # number of zero coefficients
        nnz = len(coeffs) - nz  # number of non-zero coefficients
        if(nz == 0):
            self._is_CM = [False, 0]
            return self._is_CM
        # probaly checking too many
        for D in range(3, ceil(QQ(max_nump) / QQ(2))):
            try:
                for x in DirichletGroup(D):
                    if(x.order() != 2):
                        continue
                    # we know that for CM we need x(p) = -1 => c(p)=0
                    # (for p not dividing N)
                    if(x.values().count(-1) > nz):
                        raise StopIteration()  # do not have CM with this char
                    for p in prime_range(max_nump + 1):
                        if(x(p) == -1 and coeffs[p] != 0):
                            raise StopIteration()  # do not have CM with this char
                    # if we are here we have CM with x.
                    self._is_CM = [True, x]
                    return self._is_CM
            except StopIteration:
                pass
        self._is_CM = [False, 0]
        self.is_cm = self._is_CM[0]
        return self._is_CM

    def as_polynomial_in_E4_and_E6(self,insert_in_db=True):
        r"""
        If self is on the full modular group writes self as a polynomial in E_4 and E_6.
        OUTPUT:
        -''X'' -- vector (x_1,...,x_n)
        with f = Sum_{i=0}^{k/6} x_(n-i) E_6^i * E_4^{k/4-i}
        i.e. x_i is the coefficient of E_6^(k/6-i)*
        """
        if(self.level != 1):
            raise NotImplementedError("Only implemented for SL(2,Z). Need more generators in general.")
        if(self._as_polynomial_in_E4_and_E6 is not None and self._as_polynomial_in_E4_and_E6 != ''):
            return self._as_polynomial_in_E4_and_E6
        d = self.parent.dimension_modular_forms  # dimension of space of modular forms
        k = self.weight
        K = self.base_ring
        l = list()
        wmf_logger.debug("Start to get self as polynomial in E4 and E6")
        # for n in range(d+1):
        #    l.append(self._f.q_expansion(d+2)[n])
        # v=vector(l) # (self._f.coefficients(d+1))
        v = vector(self.coefficients(range(d)))
        d = dimension_modular_forms(1, k)
        lv = len(v)
        if(lv < d):
            raise ArithmeticError("not enough Fourier coeffs")
        e4 = EisensteinForms(1, 4).basis()[0].q_expansion(lv + 2)
        e6 = EisensteinForms(1, 6).basis()[0].q_expansion(lv + 2)
        m = Matrix(K, lv, d)
        lima = floor(k / 6)  # lima=k\6;
        if((lima - (k / 2)) % 2 == 1):
            lima = lima - 1
        poldeg = lima
        col = 0
        monomials = dict()
        while(lima >= 0):
            deg6 = ZZ(lima)
            deg4 = (ZZ((ZZ(k / 2) - 3 * lima) / 2))
            e6p = (e6 ** deg6)
            e4p = (e4 ** deg4)
            monomials[col] = [deg4, deg6]
            eis = e6p * e4p
            for i in range(1, lv + 1):
                m[i - 1, col] = eis.coefficients()[i - 1]
            lima = lima - 2
            col = col + 1
        if (col != d):
            raise ArithmeticError("bug dimension")
        wmf_logger.debug("m={0}".format(m, type(m)))
        wmf_logger.debug("v={0}".format(v, type(v)))
        try:
            X = m.solve_right(v)
        except:
            return ""
        #self._as_polynomial_in_E4_and_E6 = [poldeg, monomials, X]
        return [poldeg, monomials, X]

    def exact_cm_at_i_level_1(self, N=10,insert_in_db=True):
        r"""
        Use formula by Zagier (taken from pari implementation by H. Cohen) to compute the geodesic expansion of self at i
        and evaluate the constant term.

        INPUT:
        -''N'' -- integer, the length of the expansion to use.
        """
        try:
            [poldeg, monomials, X] = self.as_polynomial_in_E4_and_E6()
        except:
            return ""
        k = self.weight
        tab = dict()
        QQ['x']
        tab[0] = 0 * x ** 0
        tab[1] = X[0] * x ** poldeg
        for ix in range(1, len(X)):
            tab[1] = tab[1] + QQ(X[ix]) * x ** monomials[ix][1]
        for n in range(1, N + 1):
            tmp = -QQ(k + 2 * n - 2) / QQ(12) * x * tab[n] + (x ** 2 - QQ(1)) / QQ(2) * ((tab[
                                                                                          n]).derivative())
            tab[n + 1] = tmp - QQ((n - 1) * (n + k - 2)) / QQ(144) * tab[n - 1]
        res = 0
        for n in range(1, N + 1):
            term = (tab[n](x=0)) * 12 ** (floor(QQ(n - 1) / QQ(2))) * x ** (n - 1) / factorial(n - 1)
            res = res + term
        
        return res


    # def print_as_polynomial_in_E4_and_E6(self):
    #     r"""

    #     """
    #     if(self.level() != 1):
    #         return ""
    #     try:
    #         [poldeg, monomials, X] = self.as_polynomial_in_E4_and_E6()
    #     except ValueError:
    #         return ""
    #     s = ""
    #     e4 = "E_{4}"
    #     e6 = "E_{6}"
    #     dens = map(denominator, X)
    #     g = gcd(dens)
    #     s = "\\frac{1}{" + str(g) + "}\left("
    #     for n in range(len(X)):
    #         c = X[n] * g
    #         if(c == -1):
    #             s = s + "-"
    #         elif(c != 1):
    #             s = s + str(c)
    #         if(n > 0 and c > 0):
    #             s = s + "+"
    #         d4 = monomials[n][0]
    #         d6 = monomials[n][1]
    #         if(d6 > 0):
    #             s = s + e6 + "^{" + str(d6) + "}"
    #         if(d4 > 0):
    #             s = s + e4 + "^{" + str(d4) + "}"
    #     s = s + "\\right)"
    #     return "\(" + s + "\)"


    def compute_cm_values_numeric(self,digits=12,insert_in_db=True):
        r"""
        Compute CM-values numerically.
        """
        if isinstance(self._cm_values,dict) and self._cm_values  <> {}:
            return self._cm_values
        if self.level<>1:
           return None
         # the points we want are i and rho. More can be added later...
        bits = ceil(int(digits) * int(4))
        CF = ComplexField(bits)
        RF = ComplexField(bits)
        eps = RF(10 ** - (digits + 1))
        if(self._verbose > 1):
            wmf_logger.debug("eps={0}".format(eps))
        K = self.base_ring
        # recall that
        degree = self.degree()
        cm_vals = dict()
        rho = CyclotomicField(3).gen()
        zi = CyclotomicField(4).gen()
        points = [rho, zi]
        maxprec = 1000  # max size of q-expansion
        minprec = 10  # max size of q-expansion
        for tau in points:
            q = CF(exp(2 * pi * I * tau))
            fexp = dict()
            cm_vals[tau] = dict()
            if(tau == I and self.level == -1):
                # cv=    #"Exact(soon...)" #_cohen_exact_formula(k)
                for h in range(degree):
                    cm_vals[tau][h] = cv
                continue
            if K.absolute_degree()==1:
                v1 = CF(0)
                v2 = CF(1)
                try:
                    for prec in range(minprec, maxprec, 10):
                        if(self._verbose > 1):
                            wmf_logger.debug("prec={0}".format(prec))
                        v2 = self.as_factor().q_eigenform(prec,names='a').truncate(prec)(q)
                        err = abs(v2 - v1)
                        if(self._verbose > 1):
                            wmf_logger.debug("err={0}".format(err))
                        if(err < eps):
                            raise StopIteration()
                        v1 = v2
                    cm_vals[tau][0] = None
                except StopIteration:
                    cm_vals[tau][0] = v2
            else:
                v1 = dict()
                v2 = dict()
                err = dict()
                for h in range(degree):
                    v1[h] = 1
                    v2[h] = 0
                try:
                    for prec in range(minprec, maxprec, 10):
                        if(self._verbose > 1):
                            wmf_logger.debug("prec={0}".format(prec))
                        c = self.coefficients(range(prec),insert_in_db=insert_in_db)
                        for h in range(degree):
                            fexp[h] = list()
                            v2[h] = 0
                            for n in range(prec):
                                cn = c[n]
                                if hasattr(cn, 'complex_embeddings'):
                                    cc = cn.complex_embeddings(CF.prec())[h]
                                else:
                                    cc = CF(cn)
                                v2[h] = v2[h] + cc * q ** n
                            err[h] = abs(v2[h] - v1[h])
                            if(self._verbose > 1):
                                wmf_logger.debug("v1[{0}]={1}".format(h,v1[h]))
                                wmf_logger.debug("v2[{0}]={1}".format(h,v2[h]))
                                wmf_logger.debug("err[{0}]={2}".format(h,err[h]))
                            if(max(err.values()) < eps):
                                raise StopIteration()
                            v1[h] = v2[h]
                except StopIteration:
                    pass
                for h in range(degree):
                    if(err[h] < eps):
                        cm_vals[tau][h] = v2[h]
                    else:
                        cm_vals[tau][h] = None
        self._cm_values = cm_vals

    
    def compute_satake_parameters_numeric(self, prec=10, bits=53,insert_in_db=True):
        r""" Compute the Satake parameters and return an html-table.

        We only do satake parameters for primes p primitive to the level.
        By defintion the S. parameters are given as the roots of
         X^2 - c(p)*X*p^((k-1)/2) + chi(p) if (p,N)=1
         (in this normalization they have absolute value 1)
         we only store the Satake parameter alpha in the upper half-plane
         (the other one is then beta=1/alpha)
        INPUT:
        -''prec'' -- compute parameters for p <=prec
        -''bits'' -- do real embedings intoi field of bits precision

        """
        if self.character.order > 2:
            ## We only implement this for trival or quadratic characters.
            ## Otherwise there is difficulty to figure out what the embeddings mean... 
            return 
        K = self.coefficient_field
        degree = K.absolute_degree()
        #wmf_logger.debug("K={0}".format(K))
        #wmf_logger.debug("degree={0}".format(degree))
        RF = RealField(bits)
        CF = ComplexField(bits)
        ps = prime_range(prec)

        self._satake['ps'] = []
        alphas = dict()
        thetas = dict()
        aps = list()
        tps = list()
        k = self.weight

        for j in range(degree):
            alphas[j] = dict()
            thetas[j] = dict()
        for j in xrange(len(ps)):
            p = ps[j]
            try:
                ap = self.coefficient(p) 
            except IndexError:
                break
            # Ignore bad primes
            if p.divides(self.level):
                continue
            self._satake['ps'].append(p)
            chip = self.character.value(p)
             
            # ap=self._f.coefficients(ZZ(prec))[p]
            if K.absolute_degree()==1:
                app = RF(ap)*p**(-0.5*(k-1))
                #wmf_logger.debug("chip={0}".format(chip))
                #wmf_logger.debug("chip.parent()={0}".format(chip.parent()))
                #wmf_logger.debug("ap={0}".format(ap))
                #wmf_logger.debug("ap.parent()={0}".format(ap.parent()))
                chip = QQ(chip)
                f1 = 4 * chip - app ** 2
                #wmf_logger.debug("f1p={0}".format(f1))
                #wmf_logger.debug("f1.parent()={0}".format(f1.parent()))
                #wmf_logger.debug("f1.complex_embeddings()={0}".format(f1.complex_embeddings()))                                
                alpha_p = (app + I * f1.sqrt()) / QQ(2)
                t_p = CF(alpha_p).argument()
                thetas[0][p] = RF(t_p)
                alphas[0][p] = CF(alpha_p)
            else:
                for jj in range(degree):
                    app = ap.complex_embeddings(bits)[jj]*p**(-0.5*(k-1))
                    #wmf_logger.debug("chip={0},{1},{2}".format(chip,type(chip),chip.parent()))
                    #wmf_logger.debug("app={0}".format(app))
                    #wmf_logger.debug("jj={0}".format(jj))            
                    if not hasattr(chip,'complex_embeddings'):
                        f1 = (4 * CF(chip)  - app ** 2)
                    else:
                        wmf_logger.debug("chip.emb={0}".format(chip.complex_embeddings(bits)))                        
                        f1 = (4 * chip.complex_embeddings(bits)[0]  - app ** 2)
                    alpha_p = (app + I * abs(f1).sqrt())/RF(2)
                    #wmf_logger.debug("f1={0}".format(f1))
                    #wmf_logger.debug("alpha_p={0}".format(alpha_p))                    
                    t_p = CF(alpha_p).argument()
                    # tps.append(t_p)
                    # aps.append(alpha_p)
                    alphas[jj][p] = CF(alpha_p)
                    thetas[jj][p] = t_p
        self._satake['alphas'] = alphas
        self._satake['thetas'] = thetas
        self._satake['alphas_latex'] = dict()
        self._satake['thetas_latex'] = dict()
        for j in self._satake['alphas'].keys():
            self._satake['alphas_latex'][j] = dict()
            for p in self._satake['alphas'][j].keys():
                s = latex(self._satake['alphas'][j][p])
                self._satake['alphas_latex'][j][p] = s
        for j in self._satake['thetas'].keys():
            self._satake['thetas_latex'][j] = dict()
            for p in self._satake['thetas'][j].keys():
                s = latex(self._satake['thetas'][j][p])
                self._satake['thetas_latex'][j][p] = s

        wmf_logger.debug("satake=".format(self._satake))
        self.satake = self._satake
        return self._satake


    def _number_of_hecke_eigenvalues_to_check(self):
        r""" Compute the number of Hecke eigenvalues (at primes) we need to check to identify twists of our given form with characters of conductor dividing the level.
        """
        ## initial bound
        bd = self.as_factor().sturm_bound()
        # we do not check primes dividing the level
        bd = bd + len(divisors(self.level))
        return bd



    def twist_by(self, x):
        r"""
        twist self by a primitive Dirichlet character x
        """
        # xx = x.primitive()
        assert x.is_primitive()
        q = x.conductor()
        # what level will the twist live on?
        level = self.level
        qq = self.character.conductor()
        new_level = lcm(self.level, lcm(q * q, q * qq))
        D = DirichletGroup(new_level)
        new_x = D(self.character) * D(x) * D(x)
        ix = D.list().index(new_x)
        #  the correct space
        NS = WebModFormSpace(self._k, new_level, ix, self._prec)
        # have to find whih form wee want
        NS.galois_decomposition()
        M = NS.sturm_bound() + len(divisors(new_level))
        C = self.coefficients(range(M))
        for label in NS._hecke_orbits_labels:
            wmf_logger.debug("label={0}".format(label))
            FT = NS.f(label)
            CT = FT.f.coefficients(M)
            wmf_logger.debug("{0}".format(CT))
            K = FT.f.hecke_eigenvalue_field()
            try:
                for n in range(2, M):
                    if(new_level % n + 1 == 0):
                        continue
                    wmf_logger.debug("n={0}".format(n))
                    ct = CT[n]
                    c = K(x(n)) * K(C[n])
                    wmf_logger.debug("{0} {1}".format(ct, c))
                    if ct != c:
                        raise StopIteration()
            except StopIteration:
                pass
            else:
                wmf_logger.debug("Twist of f={0}".format(FT))
        return FT

