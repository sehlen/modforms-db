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
from sage.all import ZZ, QQ, DirichletGroup, CuspForms, Gamma0, ModularSymbols, Newforms, trivial_character, is_squarefree, divisors, RealField, ComplexField, prime_range, I, join, gcd, Cusp, Infinity, ceil, CyclotomicField, exp, pi, primes_first_n, euler_phi, RR, prime_divisors, Integer, matrix,NumberField,PowerSeriesRing,cached_function
from sage.rings.power_series_poly import PowerSeries_poly
from sage.all import Parent, SageObject, dimension_new_cusp_forms, vector, dimension_modular_forms, dimension_cusp_forms, EisensteinForms, Matrix, floor, denominator, latex, is_prime, prime_pi, next_prime, previous_prime,primes_first_n, previous_prime, factor, loads,save,dumps,deepcopy
import re
import yaml
from flask import url_for

## DB modules


from lmfdb.modular_forms.elliptic_modular_forms import emf_logger,emf_version

from sage.rings.number_field.number_field_base import NumberField as NumberField_class
from lmfdb.modular_forms.elliptic_modular_forms.backend import connect_to_modularforms_db,get_files_from_gridfs
from wmf import wmf_logger
from wmf.web_modform_space_computing import orbit_label
#from web_modforms import WebModFormSpace_computing_class,WebModFormSpace_computing
#from web_character import WebChar

    
def WebNewForm_computing(N=1, k=2, chi=1, label='', prec=10, bitprec=53, display_bprec=26, parent=None, data=None, compute=False, verbose=-1,get_from_db=True):
    r"""
    Constructor for WebNewForms with added 'nicer' error message.
    """
    ## First check
    if chi == 1:
        if k % 2 == 1:
            wmf_logger.debug("Only zero function here with N,k,chi,label={0}.".format( (N,k,chi,label)))
            return 0
    if data is None: data = {}
    wmf_logger.debug("incoming data in construction : {0}".format(data.get('N'),data.get('k'),data.get('chi')))
    try: 
        F = WebNewForm_computing_class(N=N, k=k, chi=chi, label=label, prec=prec, bitprec = bitprec, display_bprec=display_bprec, parent = parent, data = data, compute = compute, verbose = verbose,get_from_db = get_from_db)
    except ArithmeticError as e:#Exception as e:
        wmf_logger.critical("Could not construct WebNewForm with N,k,chi,label={0}. Error: {1}".format( (N,k,chi,label),e))
        raise IndexError,"We are very sorry. The sought function could not be found in the database."
    return F


from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modforms import WebNewForm_class

class WebNewForm_computing_class(WebNewForm_class):
    r"""
    Class for representing a (cuspidal) newform on the web.
    TODO: Include the computed data in the original database so we won't have to compute here at all.
    """
    def __init__(self, N=1, k=2, chi=1, label='', prec=10, bitprec=53, display_bprec=26,parent=None, data=None,compute=False, verbose=-1,get_from_db=True):
        r"""
        Init self as form with given label in S_k(N,chi)
        """
        super(WebNewForm_computing_class,self).__init__(N,k,chi,label,prec,bitprec,display_bprec,get_from_db=False)
        print "d=",self.__dict__
        wmf_logger.debug("WebNewForm with N,k,chi,label={0}".format( (N,k,chi,label)))

        self._as_factor = None
        self._prec_needed_for_lfunctions = None
        self.set_parent()
        self.set_newform_number()
        self.compute_additional_properties()
        #self.insert_into_db()

    def compute_additional_properties(self):
        r"""
        Compute everything we need.
        """
        wmf_logger.debug("Update ap's")
        self._get_aps()
        self.coefficient_field()
        self.get_base_ring()
        self.set_dimension()
        #self.get_character_galois_orbit()
        wmf_logger.debug("compute q-expansion")
        prec = self.set_prec_needed_for_lfunctions()
        self.set_q_expansion_embeddings(prec=prec)

        wmf_logger.debug("as polynomial")
        if self._N == 1:
            self.as_polynomial_in_E4_and_E6()
        wmf_logger.debug("compute twist info")
        self.twist_info()
        wmf_logger.debug("compute CM-values")

        self.compute_cm_values_numeric()
        wmf_logger.debug("Get Atkin-Lehner evs")
        self.get_atkin_lehner_eigenvalues()
        wmf_logger.debug("compute Satake parameters")
        self.set_is_CM()
        wmf_logger.debug("compute Satake parameters")
        self.compute_satake_parameters_numeric()
        

        #c = self.coefficients(self.prec(),insert_in_db=False)
        #self._check_if_all_computed()
        self.insert_into_db()

## Functions related to storing / fetching data from database
##  
    
    def insert_into_db(self):
        r"""
        Insert a dictionary of data for self into the database collection
        WebNewforms.files
        """
        wmf_logger.debug("inserting self into db! name={0}".format(self._name))
        C = connect_to_modularforms_db('WebNewforms.files')
        fs = get_files_from_gridfs('WebNewforms')
        s = {'name':self._name,'version':float(self._version)}
        rec = C.find_one(s)
        if rec:
            id = rec.get('_id')
        else:
            id = None
        if id<>None:
            wmf_logger.debug("Removing self from db with id={0}".format(id))
            fs.delete(id)
            
        fname = "webnewform-{0:0>5}-{1:0>3}-{2:0>3}-{3}".format(self._N,self._k,self._chi,self._label) 
        d = self.to_dict()
        d.pop('_ap',None)
        d.pop('_character',None)
        d.pop('_as_factor',None)
        id = fs.put(dumps(d),filename=fname,N=int(self._N),k=int(self._k),chi=int(self._chi),label=self._label,name=self._name,version=float(self._version),character_galois_orbit=map(int,self.parent().character_galois_orbit()))
        wmf_logger.debug("inserted :{0}".format(id))
    
##  Internal functions
##

    def set_parent(self):
        from wmf.web_modform_space_computing import WebModFormSpace_computing_class,WebModFormSpace_computing
        if not isinstance(self._parent,WebModFormSpace_computing_class):
            if self._verbose > 0:
                emf_logger.debug("compute parent! label={0}".format(label))
            self._parent = WebModFormSpace_computing(self._N, self._k,self._chi,get_from_db=True,get_all_newforms_from_db=False)


        
    def set_newform_number(self):
        r"""
        Find the index of self in the Galois decomposition list. 

        """
        if self._newform_number is None:
            C = connect_to_modularforms_db('Newform_factors.files')
            res = C.find({'N':int(self.level()),'k':int(self.weight()),
                          'cchi':int(self.chi())})
            num_orbits = res.count()
            print "num_orbits=",num_orbits
            for i in range(num_orbits):
                if orbit_label(i)==self.label():
                    self._newform_number = i
                    break
        if self._newform_number is None:
            raise ValueError,"Newform with this label is not in the space!"
    def as_factor(self):
        r"""
        Return self as a newform factor
        """
        if self._as_factor is None:
            C = connect_to_modularforms_db('Newform_factors.files')
            s = {'N':int(self.level()),'k':int(self.weight()),
                 'cchi':int(self.chi()),'newform':int(self.newform_number())}
            res = C.find_one(s)
            print res
            if not res is None:
                fid = res['_id']
                fs = get_files_from_gridfs('Newform_factors')
                self._as_factor = loads(fs.get(fid).read())
        if self._as_factor is None:
            raise ValueError,"Newform matching {0} can not be found".format(s)
        return self._as_factor
            
    def get_base_ring(self):
        r"""
        The base ring of self, that is, the field of values of the character of self. 
        """
        if self._base_ring is None:
            self._base_ring = self.as_factor().base_ring()
        return self._base_ring

    
    def relative_degree(self):
        r"""
        Degree of the field of coefficient relative to its base ring.
        """
        if self._relative_degree is None:
            self._relative_degree = self.coefficient_field().absolute_degree()/self.base_ring().absolute_degree()
        return self._relative_degree

    def set_dimension(self):
        r"""
        The dimension of this galois orbit is not necessarily equal to the degree of the number field, when we have a character....
        We therefore need this routine to distinguish between the two cases...
        """
        if self._dimension is None:
            self._dimension = self.as_factor().dimension()

    def set_prec_needed_for_lfunctions(self):
        r"""
        Calculates the number of coefficients needed for the L-function
        main page (formula taken from their pages)
        """
        if self._prec_needed_for_lfunctions is None:
            self._prec_needed_for_lfunctions = 22 + int(RR(5)*RR(self._k)*RR(self._N).sqrt())
        return self._prec_needed_for_lfunctions 
    
    def set_q_expansion_embeddings(self, prec=10, bitprec=53,format='numeric',display_bprec=26,insert_in_db=True):
        r""" Compute all embeddings of self into C which are in the same space as self.
        Return 0 if we didn't compute anything new, otherwise return 1.
        """
        wmf_logger.debug("computing embeddings of q-expansions : has {0} embedded coeffs. Want : {1} with bitprec={2}".format(len(self._embeddings),prec,bitprec))
        if display_bprec > bitprec:
            display_bprec = bitprec
        ## First check if we have sufficient data
        if self._embeddings['prec'] >= prec or self._embeddings['bitprec'] >= bitprec:
            return 0 ## We should already have sufficient data.
        ## Else we compute new embeddings.
        CF = ComplexField(bitprec)
        # First wee if we need higher precision, in which case we reset all coefficients:
        if self._embeddings['bitprec'] < bitprec:
            self._embeddings['values']=[]
            self._embeddings['latex']=[]
            self._embeddings['prec']=0
        # See if we have need of more coefficients
        nstart = len(self._embeddings)
        wmf_logger.debug("Should have {0} embeddings".format(self._embeddings['prec']))
        wmf_logger.debug("Computing new stuff !")
        for n in range(self._embeddings['prec'],prec):
            try:
                cn = self.coefficient(n)
            except IndexError:
                break
            if hasattr(cn, 'complex_embeddings'):
                cn_emb = cn.complex_embeddings(bitprec)
            else:
                cn_emb = [ CF(cn) for i in range(deg) ]
            self._embeddings['values'].append(cn_emb)
        self._embeddings['prec'] = len(self._embeddings['values'])
        # See if we also need to recompute the latex strings
        if display_bprec > self._embeddings['bitprec']:
            self._embeddings['latex'] = []  ## Have to redo these
        numc = len(self._embeddings['latex'])
        for n in range(numc,prec):
            cn_emb = []
            for x in self._embeddings['values'][n]:
                t = my_complex_latex(x,display_bprec)
            cn_emb_latex.append(t)
            self._embeddidngs['latex'].append(cn_emb)
        wmf_logger.debug("has embeddings_latex:{0}".format(nstart))
        return 1
                      
    
  
    def get_atkin_lehner_eigenvalues(self):
        r"""
        Get the Atkin-Lehner eigenvalues from database if they exist. 
        """
        
        if not ((self.character().is_trivial() or self.character().order() == 2) and not self._atkin_lehner_eigenvalues is None):
            return None
        C = connect_to_modularforms_db('Atkin-Lehner.files')
        fs = get_files_from_gridfs('Atkin-Lehner')
        s = {'N':int(self.level()),'k':int(self.weight()),
             'cchi':int(self.chi()),'newform':int(self.newform_number())}
        res = C.find_one(s)
        self._atkin_lehner_eigenvalues = {}
        if not res is None:
            alid = res['_id']
            al = loads(fs.get(alid).read())
            for d in prime_divisors(self.level()):
                self._atkin_lehner_eigenvalues[d] = 1 if al[d]=='+' else -1
        else: # We compute them
            A = self.as_factor()
            for p in prime_divisors(self.level()):
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
        if(len(self._twist_info) > 0):
            return self._twist_info
        N = self.level()
        k = self.weight()
        if(is_squarefree(ZZ(N))):
            self._twist_info = [True, None ]
            return [True, None]

        # We need to check all square factors of N
        twist_candidates = list()
        KF = self.base_ring()
        # check how many Hecke eigenvalues we need to check
        max_nump = self._number_of_hecke_eigenvalues_to_check()
        maxp = max(primes_first_n(max_nump))
        for d in divisors(N):
            if(d == 1):
                continue
            # we look at all d such that d^2 divdes N
            if(not ZZ(d ** 2).divides(ZZ(N))):
                continue
            D = DirichletGroup(d)
            # check possible candidates to twist into f
            # g in S_k(M,chi) wit M=N/d^2
            M = ZZ(N / d ** 2)
            if(self._verbose > 0):
                wmf_logger.debug("Checking level {0}".format(M))
            for xig in range(euler_phi(M)):
                (t, glist) = _get_newform(M,k, xig)
                if(not t):
                    return glist
                for g in glist:
                    if(self._verbose > 1):
                        wmf_logger.debug("Comparing to function {0}".format(g))
                    KG = g.base_ring()
                    # we now see if twisting of g by xi in D gives us f
                    for xi in D:
                        try:
                            for p in primes_first_n(max_nump):
                                if(ZZ(p).divides(ZZ(N))):
                                    continue
                                bf = self.as_factor().q_eigenform(maxp + 1, names='x')[p]
                                bg = g.q_expansion(maxp + 1)[p]
                                if(bf == 0 and bg == 0):
                                    continue
                                elif(bf == 0 and bg != 0 or bg == 0 and bf != 0):
                                    raise StopIteration()
                                if(ZZ(p).divides(xi.conductor())):
                                    raise ArithmeticError("")
                                xip = xi(p)
                                # make a preliminary check that the base rings match with respect to being
                                # real or not
                                try:
                                    QQ(xip)
                                    XF = QQ
                                    if(KF != QQ or KG != QQ):
                                        raise StopIteration
                                except TypeError:
                                    # we have a  non-rational (i.e. complex) value of the character
                                    XF = xip.parent()
                                    if((KF.absolute_degree() == 1 or KF.is_totally_real()) and (KG.absolute_degre() == 1 or KG.is_totally_real())):
                                        raise StopIteration
                            ## it is diffcult to compare elements from diferent rings in general but we make some checcks
                            # is it possible to see if there is a larger ring which everything can be
                            # coerced into?
                                ok = False
                                try:
                                    a = KF(bg / xip)
                                    b = KF(bf)
                                    ok = True
                                    if(a != b):
                                        raise StopIteration()
                                except TypeError:
                                    pass
                                try:
                                    a = KG(bg)
                                    b = KG(xip * bf)
                                    ok = True
                                    if(a != b):
                                        raise StopIteration()
                                except TypeError:
                                    pass
                                if(not ok):  # we could coerce and the coefficients were equal
                                    return "Could not compare against possible candidates!"
                                # otherwise if we are here we are ok and found a candidate
                            twist_candidates.append([M, g.q_expansion(prec), xi])
                        except StopIteration:
                            # they are not equal
                            pass
        wmf_logger.debug("Candidates=v{0}".format(twist_candidates))
        self._twist_info = (False, twist_candidates)
        if(len(twist_candidates) == 0):
            self._twist_info = [True, None]
        else:
            self._twist_info = [False, twist_candidates]
        return self._twist_info

    def set_is_CM(self,insert_in_db=True):
        r"""
        Checks if f has complex multiplication and if it has then it returns the character.

        OUTPUT:

        -''[t,x]'' -- string saying whether f is CM or not and if it is, the corresponding character

        EXAMPLES::

        """
        if(len(self._is_CM) > 0):
            return self._is_CM
        max_nump = self._number_of_hecke_eigenvalues_to_check()
        # E,v = self._f.compact_system_of_eigenvalues(max_nump+1)
        try:
            coeffs = self.coefficients(range(max_nump + 1),insert_in_db=insert_in_db)
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
        return self._is_CM

    def as_polynomial_in_E4_and_E6(self,insert_in_db=True):
        r"""
        If self is on the full modular group writes self as a polynomial in E_4 and E_6.
        OUTPUT:
        -''X'' -- vector (x_1,...,x_n)
        with f = Sum_{i=0}^{k/6} x_(n-i) E_6^i * E_4^{k/4-i}
        i.e. x_i is the coefficient of E_6^(k/6-i)*
        """
        if(self.level() != 1):
            raise NotImplementedError("Only implemented for SL(2,Z). Need more generators in general.")
        if(self._as_polynomial_in_E4_and_E6 is not None and self._as_polynomial_in_E4_and_E6 != ''):
            return self._as_polynomial_in_E4_and_E6
        d = self._parent.dimension_modular_forms()  # dimension of space of modular forms
        k = self.weight()
        K = self.base_ring()
        l = list()
        # for n in range(d+1):
        #    l.append(self._f.q_expansion(d+2)[n])
        # v=vector(l) # (self._f.coefficients(d+1))
        v = vector(self.coefficients(range(d),insert_in_db=insert_in_db))
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
        # return [m,v]
        if self._verbose > 0:
            wmf_logger.debug("m={0}".format(m, type(m)))
            wmf_logger.debug("v={0}".format(v, type(v)))
        try:
            X = m.solve_right(v)
        except:
            return ""
        self._as_polynomial_in_E4_and_E6 = [poldeg, monomials, X]
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
        k = self.weight()
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


    def print_as_polynomial_in_E4_and_E6(self):
        r"""

        """
        if(self.level() != 1):
            return ""
        try:
            [poldeg, monomials, X] = self.as_polynomial_in_E4_and_E6()
        except ValueError:
            return ""
        s = ""
        e4 = "E_{4}"
        e6 = "E_{6}"
        dens = map(denominator, X)
        g = gcd(dens)
        s = "\\frac{1}{" + str(g) + "}\left("
        for n in range(len(X)):
            c = X[n] * g
            if(c == -1):
                s = s + "-"
            elif(c != 1):
                s = s + str(c)
            if(n > 0 and c > 0):
                s = s + "+"
            d4 = monomials[n][0]
            d6 = monomials[n][1]
            if(d6 > 0):
                s = s + e6 + "^{" + str(d6) + "}"
            if(d4 > 0):
                s = s + e4 + "^{" + str(d4) + "}"
        s = s + "\\right)"
        return "\(" + s + "\)"


    def compute_cm_values_numeric(self,digits=12,insert_in_db=True):
        r"""
        Compute CM-values numerically.
        """
        if isinstance(self._cm_values,dict) and self._cm_values  <> {}:
            return self._cm_values
        if self.level()<>1:
           return None
         # the points we want are i and rho. More can be added later...
        bits = ceil(int(digits) * int(4))
        CF = ComplexField(bits)
        RF = ComplexField(bits)
        eps = RF(10 ** - (digits + 1))
        if(self._verbose > 1):
            wmf_logger.debug("eps={0}".format(eps))
        K = self.base_ring()
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
            if(tau == I and self.level() == -1):
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
                        v2 = self.as_factor().q_eigenform(prec,name='a').truncate(prec)(q)
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
         X^2 - c(p)X + chi(p)*p^(k-1) if (p,N)=1

        INPUT:
        -''prec'' -- compute parameters for p <=prec
        -''bits'' -- do real embedings intoi field of bits precision

        """
        if self.character().order()>2:
            ## We only implement this for trival or quadratic characters.
            ## Otherwise there is difficulty to figure out what the embeddings mean... 
            return 
        K = self.coefficient_field()
        degree = self.degree()
        RF = RealField(bits)
        CF = ComplexField(bits)
        ps = prime_range(prec)

        self._satake['ps'] = []
        alphas = dict()
        thetas = dict()
        aps = list()
        tps = list()
        k = self.weight()

        for j in range(degree):
            alphas[j] = dict()
            thetas[j] = dict()
        for j in xrange(len(ps)):
            p = ps[j]
            try:
                ap = self.coefficient(p) 
            except IndexError:
                break
            # Remove bad primes
            if p.divides(self.level()):
                continue
            self._satake['ps'].append(p)
            chip = self.character().value(p)
            wmf_logger.debug("p={0}".format(p))
            wmf_logger.debug("chip={0} of type={1}".format(chip,type(chip)))
            if hasattr(chip,'complex_embeddings'):
                wmf_logger.debug("embeddings(chip)={0}".format(chip.complex_embeddings()))
            wmf_logger.debug("ap={0}".format(ap))
            wmf_logger.debug("K={0}".format(K))                        
            
            # ap=self._f.coefficients(ZZ(prec))[p]
            if K.absolute_degree()==1:
                f1 = QQ(4 * chip * p ** (k - 1) - ap ** 2)
                alpha_p = (QQ(ap) + I * f1.sqrt()) / QQ(2)
                ab = RF(p ** ((k - 1) / 2))
                norm_alpha = alpha_p / ab
                t_p = CF(norm_alpha).argument()
                thetas[0][p] = t_p
                alphas[0][p] = (alpha_p / ab).n(bits)
            else:
                for jj in range(degree):
                    app = ap.complex_embeddings(bits)[jj]
                    wmf_logger.debug("chip={0}".format(chip))
                    wmf_logger.debug("app={0}".format(app))
                    wmf_logger.debug("jj={0}".format(jj))            
                    if not hasattr(chip,'complex_embeddings'):
                        f1 = (4 * CF(chip) * p ** (k - 1) - app ** 2)
                    else:
                        f1 = (4 * chip.complex_embeddings(bits)[jj] * p ** (k - 1) - app ** 2)
                    alpha_p = (app + I * abs(f1).sqrt())
                    # ab=RF(/RF(2)))
                    # alpha_p=alpha_p/RealField(bits)(2)
                    wmf_logger.debug("f1={0}".format(f1))
                    
                    alpha_p = alpha_p / RF(2)
                    wmf_logger.debug("alpha_p={0}".format(alpha_p))                    
                    t_p = CF(alpha_p).argument()
                    # tps.append(t_p)
                    # aps.append(alpha_p)
                    alphas[jj][p] = alpha_p
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
        return self._satake


    def _number_of_hecke_eigenvalues_to_check(self):
        r""" Compute the number of Hecke eigenvalues (at primes) we need to check to identify twists of our given form with characters of conductor dividing the level.
        """
        ## initial bound
        bd = self.as_factor().sturm_bound()
        # we do not check primes dividing the level
        bd = bd + len(divisors(self.level()))
        return bd



    def twist_by(self, x):
        r"""
        twist self by a primitive Dirichlet character x
        """
        # xx = x.primitive()
        assert x.is_primitive()
        q = x.conductor()
        # what level will the twist live on?
        level = self.level()
        qq = self.character().conductor()
        new_level = lcm(self.level(), lcm(q * q, q * qq))
        D = DirichletGroup(new_level)
        new_x = D(self.character()) * D(x) * D(x)
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

###
### Independent helper functions
###


def my_latex_from_qexp(s):
    r"""
    Make LaTeX from string. in particular from parts of q-expansions.
    """
    ss = ""
    ss += re.sub('x\d', 'x', s)
    ss = re.sub("\^(\d+)", "^{\\1}", ss)
    ss = re.sub('\*', '', ss)
    ss = re.sub('zeta(\d+)', 'zeta_{\\1}', ss)
    ss = re.sub('zeta', '\zeta', ss)
    ss += ""
    # wmf_logger.debug("ss=",ss
    return ss


def break_line_at(s, brpt=20):
    r"""
    Breaks a line containing math 'smartly' at brpt characters.
    With smartly we mean that we break at + or - but keep brackets
    together
    """
    sl = list()
    stmp = ''
    left_par = 0
    #wmf_logger.debug('Break at line, Input ={0}'.format(s))
    for i in range(len(s)):
        if s[i] == '(':  # go to the matching case
            left_par = 1
        elif s[i] == ')' and left_par == 1:
            left_par = 0
        if left_par == 0 and (s[i] == '+' or s[i] == '-'):
            sl.append(stmp)
            stmp = ''
        stmp = stmp + s[i]
        if i == len(s) - 1:
            sl.append(stmp)
    wmf_logger.debug('sl={0}'.format(sl))

    # sl now contains a split  e.g. into terms in the q-expansion
    # we now have to join as many as fits on the line
    res = list()
    stmp = ''
    for j in range(len(sl)):
        l = len_as_printed(stmp) + len_as_printed(sl[j])
        #wmf_logger.debug("l={0}".format(l))
        if l < brpt:
            stmp = join([stmp, sl[j]])
        else:
            res.append(stmp)
            stmp = sl[j]
        if j == len(sl) - 1:
            res.append(stmp)
    return res


def _get_newform(N, k, chi, fi=None):
    r"""
    Get an element of the space of newforms, incuding some error handling.

    INPUT:

     - ''k'' -- positive integer : the weight
     - ''N'' -- positive integer (default 1) : level
     - ''chi'' -- non-neg. integer (default 0) use character nr. chi
     - ''fi'' -- integer (default 0) We want to use the element nr. fi f=Newforms(N,k)[fi]. fi=-1 returns the whole list
     - ''prec'' -- integer (the number of coefficients to get)

    OUTPUT:

    -''t'' -- bool, returning True if we succesfully created the space and picked the wanted f
    -''f'' -- equals f if t=True, otherwise contains an error message.

    EXAMPLES::


        sage: _get_newform(16,10,1)
        (False, 'Could not construct space $S^{new}_{16}(10)$')
        sage: _get_newform(10,16,1)
        (True, q - 68*q^3 + 1510*q^5 + O(q^6))
        sage: _get_newform(10,16,3)
        (True, q + 156*q^3 + 870*q^5 + O(q^6))
        sage: _get_newform(10,16,4)
        (False, '')

     """
    t = False
    try:
        if(chi == 0):
            wmf_logger.debug("EXPLICITLY CALLING NEWFORMS!")
            S = Newforms(N, k, names='x')
        else:
            S = Newforms(DirichletGroup(N)[chi], k, names='x')
        if(fi >= 0 and fi < len(S)):
            f = S[fi]
            t = True
        elif(fi == -1 or fi is None):
            t = True
            return (t, S)
        else:
            f = ""
    except RuntimeError:
        if(chi == 0):
            f = "Could not construct space $S^{new}_{%s}(%s)$" % (k, N)
        else:
            f = "Could not construct space $S^{new}_{%s}(%s,\chi_{%s})$" % (k, N, chi)
    return (t, f)


def _degree(K):
    r"""
    Returns the degree of the number field K
    """
    return K.absolute_degree()


def unpickle_wnf_v1(N, k,chi, label, fi, prec, bitprec, display_bprec,parent,data):
    F = WebNewForm(N=N,k=k, chi=chi, label=label, fi=fi, prec=prec, bitprec=bitprec, display_bprec=display_bprec,parent=parent, data=data)
    return F


def unpickle_wmfs_v1(N, k,chi, cuspidal, prec, bitprec, data):
    M = WebModFormSpace(N=N, k=k, chi=chi, cuspidal=cuspidal, prec=prec, bitprec=bitprec, data=data)
    return M


def pol_to_html(p):
    r"""
    Convert polynomial p to html.
    """
    s = str(p)
    s = re.sub("\^(\d*)", "<sup>\\1</sup>", s)
    s = re.sub("\_(\d*)", "<sub>\\1</sub>", s)
    s = re.sub("\*", "", s)
    s = re.subst("x", "<i>x</i>", s)
    return s

def pol_to_latex(p):
    r"""
    Convert polynomial in string format to latex.
    """
    s = str(p)
    s = re.sub("\^(\d*)", "^{\\1}", s)
    s = re.sub("\_(\d*)", "_{\\1}", s)
    s = re.sub("\*", "", s)
    s = re.sub("zeta(\d+)", "\zeta_{\\1}", s)
    return s


def number_field_to_dict(F):

    r"""
    INPUT:
    - 'K' -- Number Field
    - 't' -- (p,gens) where p is a polynomial in the variable(s) xN with coefficients in K. (The 'x' is just a convention)

    OUTPUT:

    - 'F' -- Number field extending K with relative minimal polynomial p.
    """
    if F.base_ring().absolute_degree()==1:
        K = 'QQ'
    else:
        K = number_field_to_dict(F.base_ring())
    if F.absolute_degree() == 1:
        p = 'x'
        g = ('x',)
    else:
        p = F.relative_polynomial()
        g = str(F.gen())
        x = p.variables()[0]
        p = str(p).replace(str(x),str(g))
    return {'base':K,'relative polynomial':p,'gens':g}


def number_field_from_dict(d):
    r"""
    INPUT:

    - 'd' -- {'base':F,'p':p,'g':g } where p is a polynomial in the variable(s) xN with coefficients in K. (The 'x' is just a convention)

    OUTPUT:

    - 'F' -- Number field extending K with relative minimal polynomial p.
    """
    K = d['base']; p=d['relative polynomial']; g=d['gens']
    if K=='QQ':
        K = QQ
    elif isinstance(K,dict):
        K = number_field_from_dict(K)
    else:
        raise ValueError,"Could not construct number field!"
    F = NumberField(K[g](p),names=g)
    if F.absolute_degree()==1:
        F = QQ
    return F

def my_complex_latex(c,bitprec):
    x = c.real().n(bitprec)
    y = c.imag().n(bitprec)
    d = floor(bitprec/3.4)
    if x >= 0:
        prefx = "\\hphantom{-}"
    else:
        prefx = ""
    if y < 0:
        prefy = ""
    else:
        prefy = "+"
    xi,xf = str(x).split(".")
    xstr = "{0}.{1:0<{d}}".format(xi,xf,d=d)
    #print "xstr=",xstr
    yi,yf = str(y).split(".")
    ystr = "{0}.{1:0<{d}}".format(yi,yf,d=d)
    t = "{prefx}{x}{prefy}{y}i".format(prefx=prefx,x=xstr,prefy=prefy,y=ystr)
    return t
#     d = 
#     if y ==0:
#         return "{0:.df}
 
