r"""

Utilities for modform computing classes.

"""

from sage.all import cached_function,AlphabeticStrings
       
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
def are_compatible(modulus,weight,ci):
    r"""
    Return true if x(-1) == (-1)**k
    where x is the ci-th character of modulus N in Sage numbering scheme.
    """
    import character_conversions
    x = character_conversions.dirichlet_character_conrey_from_sage_character_number(modulus,ci)
    return (x.is_even() and weight % 2 == 0) or (x.is_odd() and weight % 2 == 1)
