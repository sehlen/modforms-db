r"""

Utilities for modform computing classes.

"""

from sage.all import cached_function,AlphabeticStrings,QQ,Matrix
from compmf import clogger

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
    label = str(label)
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

@cached_function
def newform_label(N,k,i,d):
    if isinstance(d,(int,Integer)):
        return "{0}.{1}.{2}.{3}".format(N,k,i,orbit_label(d))
    elif isinstance(d,basestring):
        return "{0}.{1}.{2}.{3}".format(N,k,i,d)
    else:
        raise ValueError,"Could not find label from {0}".format(d)

@cached_function
def are_compatible(modulus,weight,si,si_format='orbit_no'):
    r"""
    Return true if x(-1) == (-1)**k
    where x is a representative of the si-th galois orbit of charaters of modulus N in Sage numbering scheme.
    """
    import character_conversions
    if si_format == 'orbit_no':
        i = character_conversions.dirichlet_character_conrey_galois_orbits_reps(modulus)[si]
        x = character_conversions.conrey_character_from_number(modulus,i)
    else:
        x = character_conversions.conrey_character_from_number(modulus,si)
    return (x.is_even() and weight % 2 == 0) or (x.is_odd() and weight % 2 == 1)

def label_from_param(N,k,i,d=None):
    if not d is None:
        return "{0}.{1}.{2}{3}".format(N,k,i,orbit_label(d))
    else:
        return "{0}.{1}.{2}{3}".format(N,k,i)

def param_from_label(lab):
    t = lab.split(".")
    if len(t)<>3:
        raise ValueError
    N=int(t[0]); k = int(t[1])
    s = t[2]
    if s.isdigit():
        i = int(s)
        return (N,k,i)
    else:
        l = map(lambda x:x.isalpha(),s)
        delim = l.index(True)
        i = int(s[0:delim])
        label =s[delim:]
        return (N,k,i,label)        





def multiply_mat_vec(E,v):
    EE = convert_matrix_to_extension_fld(E,v.base_ring())
    return EE*v

def convert_matrix_to_extension_fld(E,K):
    if E.base_ring() == K:
        return E
    EE=Matrix(K,E.nrows(), E.ncols())
    KE = E.base_ring()
    if KE == QQ:
        return E
    # we need to make the names agree since sometimes a field is given by a1 and sometimes by a or a2... 
    clogger.debug("KE={0}".format(KE))
    if len(K._names) == len(KE._names):
        KE = KE.change_names(K._names)
        clogger.debug("New KE={0}".format(KE))
    if KE.is_relative(): # KE = NF / Cyclotomic
        gen = KE.base_ring().gen()
    else:
        gen = KE.gen()
    z = K(gen)
    x = E[0,0].polynomial().parent().gen()
    EE=E.apply_map(lambda y: y.polynomial().substitute({x:z}))
#    for a in range(E.nrows()):
#        for b in range(E.ncols()):
#            EE[a,b]=E[a,b].polynomial().substitute({x:z})
    return EE
    
