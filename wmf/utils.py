r"""

Utilities for modform computing classes.

"""

from sage.all import cached_function
       
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
