r"""
File containing methods for converting data between different formats.
Sage < --- > Dictionary of primitive data, e.g. ints, text and blobs 

AUTHORS: Fredrik Stroemberg, Stephan Ehlen (code base from William Stein's mfdb module)


Note: We follow Conrey's convention for enumerating characters.

"""

from sage.all import cached_function,QQ,trivial_character,ModularSymbols
from dirichlet_conrey import *
import sage 
from sage.structure.sequence import Sequence

import collections
def extract_args_kwds(*args,**kwds):
    r"""
    Extract a list of keywords from arguments. The kwds takes priority.

    """
    res = []
    keywords = kwds.get('keywords',[])
    for key in keywords:
        val = kwds.get(key)
        if val==None and len(args)>0:
            for arg in args:
                if isinstance(arg,(dict,collections.MutableMapping)):
                    val = arg.get(key)
        res.append(val)
        print "res=",res
    if len(res)>1:
        return tuple(res)
    elif len(res)==1:
        return res[0]
    else:
        return tuple()

def sage_ambient_to_dict(M, i=None):
    """
    INPUT:
      - M -- instance of Sage class: ModularSymbols_ambient

    OUTPUT:
    
    Data structure of dictionary that is created:

        - 'space': (N, k, i), the level, weight, and and integer
          that specifies character in terms of table; all Python ints
        - 'eps': (character) list of images of gens in QQ or cyclo
          field; this redundantly specifies the character
        - 'manin': (manin_gens) list of 2-tuples or 3-tuples of
          integers (all the Manin symbols)
        - 'basis': list of integers -- index into above list of Manin symbols
        - 'rels': relation matrix (manin_gens_to_basis) -- a sparse
          matrix over a field (QQ or cyclotomic)
        - 'mod2term': list of pairs (n, c), such that the ith element
          of the list is equivalent via 2-term relations only
          to c times the n-th basis Manin symbol; these
    """
    if i is None:
        i = dirichlet_character_to_int(M.character(),convention='Conrey')
    res = {'space':(int(M.level()), int(M.weight()), int(i)),
            'manin':[(t.i,t.u,t.v) for t in M._manin_generators],
            'basis':M._manin_basis,
            'rels':M._manin_gens_to_basis,
            'mod2term':M._mod2term,
            'dimension_modular_forms'  : int(M.dimension()),
            'dimension_new_cusp_forms' : int(M.cuspidal_subspace().new_subspace().dimension()),
            'dimension_cusp_forms'     : int(M.cuspidal_subspace().dimension()),
             'base_field'  : base_ring_to_dict(M.base_ring()),
             'newspace_factors':[]}
    
    for N in M.cuspidal_subspace().new_subspace().decomposition():
        NN =  factor_to_dict_sage(N,with_ambient=False)
        #ModularSymbols_newspace_factor(dimension = N.dimension())
        res['newspace_factors'].append(NN)
    return res

def base_ring_to_dict(F):
    degree = F.degree()
    if hasattr(F,'is_cyclotomic'):
        is_cyclotomic=F.is_cyclotomic()
    elif isinstance(F,sage.rings.number_field.number_field.NumberField_cyclotomic):
        is_cyclotomic=True
    else:
        is_cyclotomic=False
    if degree > 1:
        minimal_polynomial = str(F.polynomial())
    else:
        minimal_polynomial = 'x'           
    return {'minimal_polynomial' : minimal_polynomial,
            'degree' : degree,
            'is_cyclotomic': is_cyclotomic }
        
def dict_to_ambient_sage(modsym,convention='Conrey'):
    r"""
    Convert a dictionary to a Sage object of type ModularSymbols_ambient

    INPUT:

    - modsym -- dictionary
    - convention -- String. 'Conrey' or 'Sage' depending on the numbering convention used for Dirichlet characters in the dict.
    """
    N,k,i = modsym['space']
    manin = modsym['manin']
    basis = modsym['basis']
    rels  = modsym['rels']
    mod2term  = modsym['mod2term']
    F = rels.base_ring()
    if i == 0:
        eps = trivial_character(N)
    else:
        if F<>QQ or convention=='Sage':
            print "i=",i
            eps = DirichletGroup(N, F)(i)
        else:
            eps = DirichletGroup_conrey(N)[i].sage_character()
    from sage.modular.modsym.manin_symbols import ManinSymbolList, ManinSymbol
    manin_symbol_list = ManinSymbolList(k, manin)

    def custom_init(M):
        # reinitialize the list of Manin symbols with ours, which may be
        # ordered differently someday:
        syms = M.manin_symbols()
        ManinSymbolList.__init__(syms, k, manin_symbol_list)

        M._manin_generators = [ManinSymbol(syms, x) for x in manin]
        M._manin_basis = basis
        M._manin_gens_to_basis = rels
        M._mod2term = mod2term
        return M

    return ModularSymbols(eps, k, sign=1, custom_init=custom_init, use_cache=False)

@cached_function
def dirichlet_character_sage_galois_orbits_reps(N):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    return [X[0] for X in DirichletGroup(N).galois_orbits()]

@cached_function
def dirichlet_character_sage_galois_orbit_rep(x):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    N = x.modulus()
    reps = sage_dirichlet_character_galois_orbits_reps(N)
    if x in reps:
        return reps.index(x)
    else: 
        for i in len(reps):
            if x in reps[i].galois_orbit():
                return i
    raise ArithmeticError("Could not find representative of Galois orbit of {0}".format(x))

@cached_function
def dirichlet_character_conrey_galois_orbit_numbers(x):
    return [y.number() for y in x.galois_orbit()]

@cached_function
def dirichlet_character_conrey_galois_orbits_reps(N):
    """
    Return list of representatives for the Galois orbits of Conrey Dirichlet characters of level N.
    We always take the one that has the smallest index.
    """
    D = DirichletGroup_conrey(N)
    Ds = dirichlet_character_sage_galois_orbit_reps(N)
    Dl = list(D)
    reps=[]
    for x in D:
        if x not in Dl:
            continue
        orbit_of_x = sorted(x.galois_orbit())
        reps.append(orbit_of_x[0])
        for xx in orbit_of_x:
            Dl.remove(xx)
    return reps

@cached_function
def dirichlet_character_sage_to_conrey(x):
    for y in DirichletGroup_conrey(x.modulus()):
        if y.sage_character()==x:
            return y

@cached_function
def dirichlet_character_conrey_used_in_computation(x):
    r"""
      INPUTS:
       - ```x```: A Conrey Dirichlet Character
    
      Returns the number of Conrey Dirichlet Character ```c```,
      such that ```c.sage_character()``` is the representative that
      was used to compute the spaces of modular forms with character ```x```.

      OUTPUT:
       - int: the number of the corresponding Conrey Dirichlet Character.
    """
    
    reps_sage = dirichlet_character_sage_galois_orbit_reps(x.modulus())

    for c in x.galois_orbit():
        if c.sage_character() in reps_sage:
            return c.number()
    

@cached_function
def dirichlet_character_conrey_galois_orbit_embeddings(x):
    r"""
       Returns a dictionary that maps the Conrey numbers
       of the Dirichlet characters in the Galois orbit of x
       to the powers of $\zeta_{\phi(N)}$ so that the corresponding
       embeddings map the labels.

       Let $\zeta_{\phi(N)}$ be the generator of the cyclotomic field
       of $N$-th roots of unity which is the base field
       for the coefficients of a modular form contained in the database.
       Considering the space $S_k(N,\chi)$, where $\chi = \chi_N(m, \cdot)$,
       if embeddings()[m] = n, then $\zeta_{\phi(N)}$ is mapped to
       $\mathrm{exp}(2\pi i n /\phi(N))$.
    """    
    embeddings = {}
    base_number = 0
    N = x.modulus()
    base_number = dirichlet_character_conrey_used_in_computation(x).number()
    embeddings[base_number] = 1
    for n in range(2,N):
        if gcd(n,N) == 1:
            embeddings[Mod(base_number,N)**n] = n
    return embeddings

def dirichlet_character_conrey_galois_orbit_rep(x):
    """
    Return a representative of the Galois orbit of the Dirichlet character x
    """
    N = x.modulus()
    reps = dirichlet_character_conrey_galois_orbits_reps(N)
    ## Change this
    xx = dirichlet_character_sage_to_conrey(x)
    if xx in reps:
        return xx
        #return reps.index(xx)
    else: 
        for i in range(len(reps)):
            if xx in reps[i].galois_orbit():
                return reps[i]
    raise ArithmeticError('Did not find representative of {0}'.format(xx))

@cached_function
def dirichlet_character_to_int(chi,convention='Conrey'):
    r"""
    Returns integer representing the character 
    """
    N = chi.modulus()
    if convention=='Sage':
        x = dirichlet_character_sage_galois_orbit_rep(chi)
        return dirichlet_character_sage_galois_orbit_reps(chi.modulus()).index(x)
    elif convention=='Conrey':
        x = dirichlet_character_conrey_galois_orbit_rep(chi)
        return x.number()
    else:
        raise ValueError("convention must be one of 'Sage' or 'Conrey' ")


def dict_to_factor_sage(factor,**kwds):
    r"""
    Returns a newform factor as in instance of the Sage class ModularSymbolsSubspace
    """
    print "factor=",factor
    M = factor.get('ambient',kwds.get('ambient',None))
    if M == None:
        raise ValueError('Dictionary does not contain an ambient space!')
    B = factor.get('B',None)
    Bd = factor.get('Bd',None)
    v = factor.get('v',None)
    nz = factor.get('nz',None)
    if M is None or B is None or Bd is None or v is None or nz is None:
        raise ValueError('Dictionary does not contain correct data!')
    M = dict_to_ambient_sage(M)
    B._cache['in_echelon_form'] = True
    Bd._cache['in_echelon_form'] = True
    # These two lines are scary, obviously, since they depend on
    # private attributes of modular symbols.
    A = sage.modular.modsym.subspace.ModularSymbolsSubspace(M, B.row_module(), Bd.row_module(), check=False)
    A._HeckeModule_free_module__dual_eigenvector = {('a',nz):(v,False)}
    A._is_simple = True
    A._HeckeModule_free_module__decomposition = {(None,True):Sequence([A], check=False)}
    return A

def factor_to_dict_sage(factor,with_ambient=False):
    r"""
    Take a newform factor of the type ModularSymbolSubspace and return a dictionary.
    THe parameter with_ambient is there to avoid recursion.
    """
#    import sage.modular.modsym.subspace
    if with_ambient:
        M  = sage_ambient_to_dict(factor.ambient())
    B  = factor.free_module().basis_matrix()
    Bd = factor.dual_free_module().basis_matrix()
    v  = factor.dual_eigenvector(names='a', lift=False) 
    nz = factor._eigen_nonzero()
    print "factor=",factor,type(factor)
    print "B=",B,type(B)
    print "B._cache=",B._cache
    B._cache = {}
    B._cache['in_echelon_form'] = True
    Bd._cache['in_echelon_form'] = True   
    if with_ambient:
        res = {'B':B,'Bd':Bd,'v':v,'nz':nz,'ambient':M}
    else:
        res = {'B':B,'Bd':Bd,'v':v,'nz':nz}
    return res




def an_from_aps(n,aplist):
    pass
    
