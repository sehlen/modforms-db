r"""
File containing methods for converting data between different formats.
Sage < --- > Dictionary of primitive data, e.g. ints, text and blobs 

AUTHORS: Fredrik Stroemberg, Stephan Ehlen (code base from William Stein's mfdb module)


Note: We follow Conrey's convention for enumerating characters.

"""

from sage.all import cached_function,QQ,trivial_character,ModularSymbols
from dirichlet_conrey import *
 

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
    return {'space':(int(M.level()), int(M.weight()), int(i)),
            'manin':[(t.i,t.u,t.v) for t in M._manin_generators],
            'basis':M._manin_basis,
            'rels':M._manin_gens_to_basis,
            'mod2term':M._mod2term}

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
            eps = DirichletGroup(N, F)(eps)
        else:
            eps = DirichletGroup_conrey(N)[eps].sage_character()
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
def dirichlet_character_conrey_galois_orbits_reps(N):
    """
    Return list of representatives for the Galois orbits of Conrey Dirichlet characters of level N.
    """
    D = DirichletGroup_conrey(N)
    Dl = list(D)
    reps=[]
    for x in D:
        if x not in Dl:
            continue
        orbit_of_x=x.galois_orbit()        
        reps.append(orbit_of_x[0])
        for xx in orbit_of_x:
            Dl.remove(xx)
    return reps

@cached_function
def dirichlet_character_sage_to_conrey(x):
    for y in DirichletGroup_conrey(x.modulus()):
        if y.sage_character()==x:
            return y

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


def dict_to_factor_sage(factor):
    r"""
    Returns a newform factor as in instance of the Sage class ModularSymbolsSubspace
    """
    M = factor.get('ambient',None)
    B = factor.get('B',None)
    Bd = factor.get('Bd',None)
    d = factor.get('d',None)
    v = factor.get('v',None)
    nz = factor.get('nz',None)
    if M is None or B is None or Bd is None or d is None or v is None or nz is None:
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

def factor_to_dict_sage(factor):
    r"""
    Take a newform factor of the type ModularSymbolSubspace and return a dictionary.
    """
    import sage.modular.modsym.subspace
    M  = sage_ambient_to_dict(factor.ambient())
    B  = factor.free_module().basis_matrix()
    Bd = factor.dual_free_module().basis_matrix()
    v  = factor.dual_eigenvector(names='a', lift=False) 
    nz = factor._eigen_nonzero()
    B._cache['in_echelon_form'] = True
    Bd._cache['in_echelon_form'] = True   
    res = {'ambient':M,'B':B,'Bd':Bd,'v':v,'nz':nz}
    return res




def an_from_aps(n,aplist):
    pass
    
