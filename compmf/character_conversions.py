# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2014
#  Fredrik Strömberg <fredrik314@gmail.com>,
#  Stephan Ehlen <stephan.j.ehlen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPLv2)
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
r"""

Routines for converting characters between the Sage characters and the Conrey characters and their Galois orbits.


AUTHORS: Fredrik Stromberg, Stephan Ehlen (code base from William Stein's mfdb module)

TODO: Do we want to store characters in a database instead of a simple cache? 

Available functions:

Groups of Dirichlet Characters:
 - dirichlet_group_conrey(n)
 - dirichlet_group_sage(n)

Convert between characters and numbers (in same system) 
 - conrey_character_from_number(n,i) -- return a Conrey character from a conrey number
 - conrey_character_to_number(x) -- return the number of character 
 - sage_character_from_number(n,i)
 - sage_character_to_number(x)
 - sage_character_to_conrey_galois_orbit_number(x)
 - sage_character_to_conrey_character(x)
 - conrey_character_number_to_conrey_galois_orbit_number(n,i)
 

Galois orbits:
 - dirichlet_group_conrey_galois_orbits(n)
 - dirichlet_character_sage_galois_orbits_reps(N)
 - dirichlet_character_conrey_galois_orbits_reps(N)

- sage_character_to_sage_galois_orbit_number(x) -- take a Sage character and return its index in the Galois orbit given by DirichletGroup(n).galois_orbits(reps_only=True)

 - dirichlet_character_conrey_galois_orbit_numbers_from_character_number(n,xi)
 - dirichlet_character_conrey_galois_orbit_rep_from_character_number(n,xi)
 - dirichlet_character_conrey_from_sage_character_number(n,i)
 - dirichlet_character_conrey_from_sage_character(x):
 - dirichlet_character_conrey_from_sage_galois_orbit_number(n,i):
 - dirichlet_character_sage_from_conrey_character(x):
 - dirichlet_character_sage_from_conrey_character_number(n,i):
 - dirichlet_character_conrey_galois_orbit_rep_from_sage_character(x):
 - dirichlet_character_sage_galois_orbit_rep_from_number(N,i):
 - dirichlet_character_sage_galois_orbit_rep_from_sage_character(N,xi):
 - conrey_character_number_from_sage_galois_orbit_number(n,i):
 - sage_character_number_from_conrey_number(n,i):
 - dirichlet_character_conrey_used_in_computation(N,xi):
 - dirichlet_character_conrey_galois_orbit_embeddings(N,xi):
 - dirichlet_character_to_int(chi,convention='Conrey'):

 - conrey_orbit_number_
"""

from sage.all import cached_function,QQ,trivial_character,ModularSymbols,Mod
from dirichlet_conrey import *
import sage 
from sage.structure.sequence import Sequence
from compmf import clogger
## The Dirichlet group

@cached_function
def dirichlet_group_conrey(n):
    r"""
    Return DirichletGroup_conrey(n) 
    """
    if n > 0:
        return DirichletGroup_conrey(n)
    raise ValueError,"No Dirichlet Group of modulus 0!"

@cached_function
def dirichlet_group_sage(n):
    r"""
    Return DirichletGroup(n) 
    """
    if n>0:
        return DirichletGroup(n)
    raise ValueError,"No Dirichlet Group of modulus 0!"        

# Characters to and from number
# Note that Conrey character correspond to their actual number, i.e.
# x <---> x.number()
# And Sage characters correspond to their index in their associated DirichletGroup
# x <---> DirichletGroup(N).index(x)
@cached_function    
def conrey_character_from_number(n,i):
    r"""
    Return character nr. i in DirichletGroup_conrey(n) 
    """
    return dirichlet_group_conrey(n)[i]
    
def conrey_character_to_number(x):
    r"""
    Return the odulus and number of the character nr. x
    """
    if x.modulus() == 1:
        return 1,1
    return x.modulus(),x.number()

@cached_function    
def sage_character_from_number(n,i):
    r"""
    Return character nr. i in DirichletGroup(n) 
    """
    return dirichlet_group_sage(n)[i]

def sage_character_to_number(x):
    r"""
    Return the modulus and number of the character nr. x
    """
    try:
        n = x.modulus()
        i = dirichlet_group_sage(n).list().index(x)
    except Exception as e:
        raise ValueError,"{0} is not a sage character! Error:{1}".format(x,e.message)
    return n,i

def sage_character_to_sage_galois_orbit_number(x):
    r"""
    Return the number of the (Sage) galois orbit which contains x
    """
    if x.is_trivial():
        return int(0)
    orbit = x.galois_orbit()
    i = 0
    for c in dirichlet_character_sage_galois_orbits_reps(x.modulus()):
        if c in orbit:
            return i
        i+=1
    raise ValueError,"Could not find Galois orbit of {0}".format(x)

def sage_character_to_conrey_character(x):
    r"""
    Return the Corney character corresponding to the sage character x
    """
    for y in dirichlet_group_conrey(x.modulus()):
        if x == y.sage_character():
            return y

def sage_character_to_conrey_galois_orbit_number(x):
    r"""
    Return the number of the Galois orbit in DirichletGroup_conrey whjich contains x
    """
    i = 0
    for o in dirichlet_group_conrey_galois_orbits(x.modulus()):
        for y in o:
            if y.sage_character() == x:
                return i
        i+=1
    raise ArithmeticError,"Could not find an orbit for x={0}".format(x)


def conrey_character_number_to_conrey_galois_orbit_number(n,i):
    r"""
    Return the number of the galois orbit in DirichletGroup_conrey(n).galois_orbits()
    which contains the Dirichlet character DirichletGroup_conrey(n)[i]
    """
    if n == 1:
        if i == 1 or i == 0:
            return 1
        else:
            raise ArithmeticError,"Could not find an orbit for (n,i)={0}".format((n,i))    
    j = 0
    for o in dirichlet_group_conrey_galois_orbits(n):
        for x in o:
            if x.number() == i:
                return j
        j+=1
    raise ArithmeticError,"Could not find an orbit for (n,i)={0}".format((n,i))    
# Galois orbits

@cached_function
def dirichlet_group_conrey_galois_orbits(n):
    if n > 0:
        return dirichlet_group_conrey(n).galois_orbits()
    raise ValueError,"No Dirichlet Group of modulus 0!"

@cached_function
def dirichlet_group_sage_galois_orbits(n):
    r"""
    Returns the Galois orbits of the (Sage) DirichletGroup(n)
    in the same ordering as given by the reps_only=True option.
    """
    if n == 0:
        raise ValueError,"No Dirichlet Group of modulus 0!"        
    l = dirichlet_group_sage(n).galois_orbits(reps_only=True)
    res = []
    for x in l:
        res.append(x.galois_orbit())
    return res
    


## Galois orbit representatives
        
@cached_function
def dirichlet_character_sage_galois_orbits_reps(N):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    return DirichletGroup(N).galois_orbits(reps_only=True)

@cached_function
def dirichlet_character_conrey_galois_orbits_reps(N):
    """
    Return list of representatives for the Galois orbits of Conrey Dirichlet characters of level N.
    We always take the one that has the smallest index.
    """
    D = DirichletGroup_conrey(N)
    if N == 1:
        return [D[1]]
    Ds = dirichlet_character_sage_galois_orbits_reps(N)
    Dl = list(D)
    reps=[]
    for x in D:
        if x not in Dl:
            continue
        orbit_of_x = sorted(x.galois_orbit())
        reps.append(orbit_of_x[0])
        for xx in orbit_of_x:
            if xx not in Dl:
                continue
            Dl.remove(xx)
    return reps
    

## Representatives as numbers

@cached_function
def dirichlet_character_conrey_galois_orbit_numbers_from_character_number(n,xi):
    r"""
    Return the numbers of the characters in the galois orbit of the conrey character chi(n,xi)
    """
    if n==1:
        return [int(1)]
    x = conrey_character_from_number(n,xi)
    orbit = [y.number() for y in x.galois_orbit()]
    orbit.sort()
    return orbit

@cached_function    
def dirichlet_character_conrey_galois_orbit_rep_from_character_number(n,xi):
    r"""
    Return the representative of the Galois orbit nr. i modulo N
    """    
    D = dirichlet_character_conrey_galois_orbits_reps(n)
    if n == 1:
        if xi == 0 or xi == 1:
            return D[0]
        else:
            raise ValueError,"Character nr {0} mod {1} does not have a rep!".format(xi,n)
    orbit = dirichlet_character_conrey_galois_orbit_numbers_from_character_number(n,xi)
    for x in D:
        if x.number() in orbit:
            return x
    else:
        raise ValueError,"Character nr {0} mod {1} does not have a rep!".format(xi,n)


### Conversions between the two types 



@cached_function    
def dirichlet_character_conrey_from_sage_character_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    x = sage_character_from_number(n,i)
    for c in dirichlet_group_conrey(n):
        if c.sage_character() == x:
            return c

@cached_function
def dirichlet_character_conrey_from_sage_character(x):
    r"""
    Return the Conrey character corresponding to DirichletGroup(n)[xi]
    """
    n,xi = sage_character_to_number(x)
    return dirichlet_character_conrey_from_sage_character_number(n,xi)

@cached_function    
def dirichlet_character_conrey_from_sage_galois_orbit_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    D = dirichlet_group_sage(n)
    # We need to find the correct orbit as given by the ordering in .galois_orbits(reps_only=True)
    x = D.galois_orbits(reps_only=True)[i]
    for c in dirichlet_group_conrey(n):
        if c.sage_character() in x.galois_orbit():
            return c
    raise ValueError,"No Conrey character for {0},{1}".format(n,i)
    
def dirichlet_character_sage_from_conrey_character(x):
    r"""
    Return the Conrey character corresponding to DirichletGroup(n)[xi]
    """
    n,xi = conrey_character_to_number(x)
    return dirichlet_character_sage_from_conrey_character_number(n,xi)

@cached_function    
def dirichlet_character_sage_from_conrey_character_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    x = conrey_character_from_number(n,i)
    for c in dirichlet_group_sage(n):
        if c==x.sage_character():
            return c
    raise ValueError,"No rep for Conrey number {0},{1}".format(n,i)
    
def dirichlet_character_conrey_galois_orbit_rep_from_sage_character(x):
    r"""
    Return the representative of the Galois orbit nr. i modulo N
    """    
    n,i = sage_character_to_number(x)
    y = dirichlet_character_conrey_from_sage_character_number(n,i)
    rep = dirichlet_character_conrey_galois_orbit_rep_from_character_number(n,x.number())
    return rep
    
@cached_function    
def dirichlet_character_sage_galois_orbit_rep_from_number(N,i):
    r"""
    Return the representative of the Galois orbit nr. i modulo N
    """    
    D = dirichlet_character_sage_galois_orbits_reps(N)
    if i>=0 and i < len(D):
        return D[i]
    else:
        raise ValueError,"Galois orbit nr. {0} does not exist in DirichletGroup({1})!".format(i,N)

@cached_function
def dirichlet_character_sage_galois_orbit_rep_from_sage_character(N,xi):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N (using the sage representation of characters).
    """
    if N == 1:
        return 1
    x = conrey_character_from_number(N,xi).sage_character()
    reps = dirichlet_character_sage_galois_orbits_reps(N)
    #clogger.debug("reps={0}".format(reps))
    clogger.debug("reps={0}".format(reps))
    if x in reps:
        return x #reps.index(x)
    else: 
        for i in range(len(reps)):
            if x in reps[i].galois_orbit():
                return reps[i]
    raise ArithmeticError("Could not find representative of Galois orbit of {0}".format(x))

@cached_function
def conrey_character_number_from_sage_galois_orbit_number(n,i):
    r"""
    Get the number of the character x corresponding to the (Sage) Galois orbit nr. i in DirichletGroup(n)
    """
    if n == 1:
        return 1
    return dirichlet_character_conrey_from_sage_galois_orbit_number(n,i).number()

@cached_function
def conrey_galois_orbit_number_from_sage_galois_orbit_number(n,i):
    r"""
    Get the number of the character x corresponding to the (Sage) Galois orbit nr. i in DirichletGroup(n)
    """
    D = dirichlet_group_conrey(n)
    cchar = D[conrey_character_number_from_sage_galois_orbit_number(n,i)]
    i = 0
    for orbit in D.galois_orbits():
        if cchar in orbit:
            return i
        i+=1
    raise ValueError,"Could not find galois orbit number!"
    #return dirichlet_character_conrey_from_sage_galois_orbit_number(n,i).number()



@cached_function
def sage_character_number_from_conrey_number(n,i):
    r"""
    Get the number of the orbit of the character x corresponding to the Conrey character nr. i.
    """
    x = dirichlet_character_sage_from_conrey_character_number(n,i)
    orbit = x.galois_orbit()
    i = 0
    for c in dirichlet_character_sage_galois_orbits_reps(n):
        if c in orbit:
            return i
        i+=1
    raise ArithmeticError("Could not find sage character number of COnrey character {0},{1}".format(n,i))
    

@cached_function
def dirichlet_character_conrey_used_in_computation(N,xi):
    r"""
      INPUTS:
       - ```x```: A Conrey Dirichlet Character
    
      Returns the number of Conrey Dirichlet Character ```c```,
      such that ```c.sage_character()``` is the representative that
      was used to compute the spaces of modular forms with character ```x```.

      OUTPUT:
       - int: the number of the corresponding Conrey Dirichlet Character.
    """
    if N == 1:
        return int(1)
    reps_sage = dirichlet_character_sage_galois_orbits_reps(N)
    x = conrey_character_from_number(N,xi)
    for c in x.galois_orbit():
        if c.sage_character() in reps_sage:
            return c.number()
    

@cached_function
def dirichlet_character_conrey_galois_orbit_embeddings(N,xi):
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
    x = conrey_character_from_number(N,xi)
    #N = x.modulus()
    base_number = dirichlet_character_conrey_used_in_computation(N,xi)
    embeddings[base_number] = 1
    for n in range(2,N):
        if gcd(n,N) == 1:
            embeddings[Mod(base_number,N)**n] = n
    return embeddings


def dirichlet_character_to_int(chi,convention='Conrey'):
    r"""
    Returns integer representing the character 
    """
    N = chi.modulus()
    if convention=='Sage':
        x = dirichlet_character_sage_galois_orbit_rep(chi)
        return dirichlet_character_sage_galois_orbits_reps(chi.modulus()).index(x)
    elif convention=='Conrey':
        x = dirichlet_character_conrey_galois_orbit_rep(chi)
        if x.modulus() == 1:
            return 1
        return x.number()
    else:
        raise ValueError("convention must be one of 'Sage' or 'Conrey' ")


def test_conversion(nmax=10):
    r"""

    """
    for n in range(1,nmax):
        D = dirichlet_group_sage(n)
        for x in D:
            y = dirichlet_character_conrey_from_sage_character(x)
            z = dirichlet_character_sage_from_conrey_character(y)
            if x<>z:
                raise ArithmetiError,"Problem in n={0}, x={1}".format(n,x)
    return True      
