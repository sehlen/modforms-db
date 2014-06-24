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

"""

from sage.all import cached_function,QQ,trivial_character,ModularSymbols,Mod
from dirichlet_conrey import *
import sage 
from sage.structure.sequence import Sequence


@cached_function
def dirichlet_group_conrey(n):
    r"""
    Return DirichletGroup_conrey(n) 
    """
    return DirichletGroup_conrey(n)
@cached_function
def dirichlet_group_sage(n):
    r"""
    Return DirichletGroup(n) 
    """
    return DirichletGroup(n)
@cached_function    
def conrey_character_from_number(n,i):
    r"""
    Return character nr. i in DirichletGroup_conrey(n) 
    """
    return dirichlet_group_conrey(n)[i]
    
@cached_function    
def conrey_from_sage_character_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    D = dirichlet_group_sage(n)
    x = D.galois_orbits()[i][0]
    for c in dirichlet_group_conrey(n):
        if c.sage_character() == x:
            return c
    
@cached_function
def dirichlet_character_sage_galois_orbits_reps(N):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    return [X[0] for X in DirichletGroup(N).galois_orbits()]

@cached_function    
def sage_character_galois_orbit_rep_from_number(N,i):
    r"""
    Return the representative of the Galois orbit nr. i
    """    
    D = dirichlet_character_sage_galois_orbits_reps(N)
    if i>=0 and i < len(D):
        return D[i]
    else:
        raise ValueError,"Galois orbit nr. {0} does not exist in DirichletGroup({1})!".format(i,N)

@cached_function
def dirichlet_character_sage_galois_orbit_rep(N,xi):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    if n == 1:
        return 1
    x = conrey_character_from_number(N,xi)
    reps = sage_dirichlet_character_galois_orbits_reps(N)
    if x in reps:
        return reps.index(x)
    else: 
        for i in len(reps):
            if x in reps[i].galois_orbit():
                return i
    raise ArithmeticError("Could not find representative of Galois orbit of {0}".format(x))

            
@cached_function
def dirichlet_character_conrey_galois_orbit_numbers(n,xi):
    r"""
    Return the numbers of the characters in the galois orbit of the conrey character chi(n,xi)
    """
    if n==1:
        return [int(1)]
    x = conrey_character_from_number(n,xi)
    return [y.number() for y in x.galois_orbit()]

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
            Dl.remove(xx)
    return reps

@cached_function
def dirichlet_character_sage_to_conrey(n,xi):
    r"""
    Return the Conrey character corresponding to DirichletGroup(n)[xi]
    """
    x = conrey_from_sage_character_number(n,xi)
    for y in dirichlet_group_conrey(n):
        if y.sage_character()==x:
            return y

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

def dirichlet_character_conrey_galois_orbit_rep(x):
    """
    Return a representative of the Galois orbit of the Dirichlet character x
    """
    N = x.modulus()
    if N == 1:
        return x
    reps = dirichlet_character_conrey_galois_orbits_reps(N)
    for i in range(len(reps)):
        if x in reps[i].galois_orbit():
            return reps[i]
    raise ArithmeticError('Did not find representative of {0}'.format(x))


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
        return x.number()
    else:
        raise ValueError("convention must be one of 'Sage' or 'Conrey' ")


