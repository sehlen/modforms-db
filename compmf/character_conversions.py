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
 - dirichlet_group_conrey_galois_orbits(n)
 - dirichlet_group_sage(n)
 - dirichlet_group_sage_galois_orbits(n)

 - dirichlet_character_sage_galois_orbits_reps(N)
 - dirichlet_character_conrey_galois_orbits_reps(N)

Convert between characters and numbers (in same system) 
 - conrey_character_from_number(n,i) -- return a Conrey character from a conrey number
 - conrey_character_to_number(x) -- return the number of character 
 - conrey_character_from_galois_orbit_number(n,i)
 - conrey_character_to_galois_orbit_number(x):

 - sage_character_from_number(n,i)
 - sage_character_to_number(x)
 - sage_character_from_sage_galois_orbit_number(n,i)
 - sage_character_to_sage_galois_orbit_number(x)

##
## Primitive conversion functions
##
 - sage_to_sage_(x,number_format='galois_orbit_number',output='character')
 - conrey_to_conrey_(x,number_format='character_number',output='character')


 - sage_character_to_conrey_character(x)
 - sage_character_to_conrey_galois_orbit_number(x)
 - conrey_character_number_to_conrey_galois_orbit_number(n,i)


 - dirichlet_character_conrey_galois_orbit_numbers_from_character_number(n,xi)
 - dirichlet_character_conrey_galois_orbit_rep_from_character_number(n,xi)
 - dirichlet_character_conrey_from_sage_character_number(n,i)
 - dirichlet_character_conrey_from_sage_galois_orbit_number(n,i):

 - dirichlet_character_sage_from_conrey_character_number(n,i):
 - dirichlet_character_sage_from_conrey_galois_orbit_number(n,i):

 - dirichlet_character_conrey_galois_orbit_rep_from_sage_character(x):
 - dirichlet_character_sagey_galois_orbit_rep_from_sage_character(x):
 - conrey_character_number_from_sage_galois_orbit_number(n,i):
 - conrey_galois_orbit_number_from_sage_galois_orbit_number(n,i)


## Should not be used since it is not always correct
## - dirichlet_character_conrey_used_in_computation(N,xi):
 - dirichlet_character_conrey_galois_orbit_embeddings(N,xi):
 - dirichlet_character_to_int(chi,convention='Conrey'):

 - conrey_orbit_number_

### Basic conversion function:
- _sage_to_sage_(x,input='character',output='character')
"""

from sage.all import cached_function,QQ,trivial_character,ModularSymbols,Mod
from dirichlet_conrey import *
import sage 
from sage.structure.sequence import Sequence
from compmf import clogger
## The Dirichlet group
from dirichlet_conrey import DirichletCharacter_conrey

@cached_function
def dirichlet_group_conrey(n):
    r"""
    Return DirichletGroup_conrey(n) 
    """
    if n > 0:
        return DirichletGroup_conrey(n)
    raise ValueError,"No Dirichlet Group of modulus 0!"

@cached_function
def dirichlet_group_conrey_galois_orbits(n):
    if n > 0:
        return dirichlet_group_conrey(n).galois_orbits()
    raise ValueError,"No Dirichlet Group of modulus 0!"


@cached_function
def dirichlet_group_sage(n):
    r"""
    Return DirichletGroup(n) 
    """
    if n>0:
        return DirichletGroup(n)
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
    D = dirichlet_group_conrey(N)
    if N == 1:
        return [D[1]]
    reps = []
    for o in D._galois_orbits():
        reps.append(min(o))
    return reps
    


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
    return conrey_to_conrey_((n,i),number_format='character_number',output='character')
    
def conrey_character_to_number(x):
    r"""
    Return the modulus and number of the character nr. x
    """
    return conrey_to_conrey_(x,output='character_number')

@cached_function    
def conrey_character_from_galois_orbit_number(n,i):
    r"""
    Return character nr. i in DirichletGroup_conrey(n) 
    """
    return conrey_to_conrey_((n,i),number_format='galois_orbit_number',output='character')
    
def conrey_character_to_galois_orbit_number(x):
    r"""
    Return the modulus and number of the character nr. x
    """
    return conrey_to_conrey_(x,output='galois_orbit_number')

@cached_function    
def sage_character_from_number(n,i):
    r"""
    Return character nr. i in DirichletGroup(n) 
    """
    return sage_to_sage_((n,i),number_format='character_number',output='character')

def sage_character_to_number(x):
    r"""
    Return the modulus and number of the character nr. x
    """
    return sage_to_sage_(x,output='character_number')

def sage_character_from_sage_galois_orbit_number(n,i):
    r"""
    Return the number of the (Sage) galois orbit which contains x
    """
    return sage_to_sage_((n,i),number_format='galois_orbit_number',output='character')

def sage_character_to_sage_galois_orbit_number(x):
    r"""
    Return the number of the (Sage) galois orbit which contains x
    """
    return sage_to_sage_(x,output='galois_orbit_number')

## The two main conversion functions
def sage_to_sage_(x,number_format='galois_orbit_number',output='character'):
    r"""
    Return the number of the (Sage) galois orbit which contains x
    INPUT:
     - x -- tuple (N,i) or DirichletCharacter
     - 'output' -- string (in ['character_number','galois_orbit','galois_orbit_number']
    """
    
    if output not in ['character','galois_orbit','galois_orbit_number']:
        raise ValueError,"Conversion from Sage character to {0} is not implemented!".format(output)
    if not isinstance(x,(tuple,list)) and not isinstance(x,type(trivial_character(1))):
        raise ValueError,"Conversion from Sage character in format {0} to {1} is not implemented!".format(x,output)           
    if isinstance(x,tuple):
        if len(x) <> 2 or not isinstance(x[0],(int,Integer)) or not isinstance(x[1],(int,Integer)):
            raise ValueError,"Input: {0} is not a tuple of the form (n,i) !".format(x)
        N,i = x
        if N <= 0:
            raise ValueError,"No Dirichlet character of modulus {0}".format(N)
        if number_format == 'galois_orbit_number':
            num_orbits = len(dirichlet_character_sage_galois_orbits_reps(N))
            if i<0 or i > num_orbits:
                raise ValueError,"Input: {0} is not an index of a galois orbit of characters modulo {1}. We only have {2} orbits.".format(i,N,num_orbits)
            x = dirichlet_character_sage_galois_orbits_reps(N)[i]            
        if number_format == 'character_number':
            D = dirichlet_group_sage(N)
            if i < 0 or i > len(D):
                raise ValueError,"Input: {0} is not an index of a character modulo {1}. We only have {2} orbits.".format(i,N,len(D))                
            x = D[i]
    if output == 'character':
        return x
    if output == 'character_number':
        return x.parent().list().index()
    if output == 'galois_orbit':
        return x.galois_orbit()
    i = 0
    for c in dirichlet_group_sage_galois_orbits(x.modulus()):
        if x in c:
            if output == 'galois_orbit_rep':
                return dirichlet_character_sage_galois_orbits_reps(N)[i] 
            else:
                return x.modulus(),i
        i+=1

    raise ValueError,"Could not find Galois orbit of {0}".format(x)

def conrey_to_conrey_(x,number_format='character_number',output='character'):
    r"""
    Convert a Conrey character form one format to another.
    (as an instance of DirichletCharacter_conrey or a tuple (n,i)
    INPUT:
     - x -- tuple (N,i) or DirichletCharacter_conrey. In the first case N is the modulus and
    i is either a character number or a Galois orbit number as indicated by .
     - 'number_format' -- string (in ['character_number','galois_orbit_number','galois_orbit_rep']
     - 'output' -- string (in ['character_number','galois_orbit','galois_orbit_number']
    """
    if output not in ['character','character_number','galois_orbit','galois_orbit_number','galois_orbit_rep','galois_orbit_numbers_list']:
        raise ValueError,"Conversion from Sage character to {0} is not implemented!".format(output)
    if not isinstance(x,(tuple,list)) and not isinstance(x, DirichletCharacter_conrey):
        raise ValueError,"Conversion from Conrey character in format {0} to {1} is not implemented!".format(x,output)           
    if isinstance(x,tuple):
        if len(x) <> 2 or not isinstance(x[0],(int,Integer)) or not isinstance(x[1],(int,Integer)):
            raise ValueError,"Input: {0} is not a tuple of the form (n,i) !".format(x)
        N,i  = x
        if N<1:
            raise ValueError,"There is no Dirichlet group of modulus {0}!".format(N)
        D =  dirichlet_group_conrey(N)
        if number_format == 'character_number':
            numbers = map(lambda y : y.number(),D)
            if not i in numbers:
                raise ValueError,"There is no Character of number {0} in {1}!".format(i,D)
            x = DirichletCharacter_conrey(D,i)
        elif number_format == 'galois_orbit_number':
            if i> len(D._galois_orbits()) or i <0:
                    raise ValueError,"There is no Character of number {0} in {2}!".format(i,D)
            x = dirichlet_group_conrey_galois_orbits_reps(D,i)
        else:
            raise ValueError,"There is no format {0} implemented!".format(number_format)
    if output == 'character':
        return x
    if output == 'character_number':
        if x.modulus()==1:
            return int(1)
        return x.number()
    o = x.galois_orbit()
    o.sort()
    if output == 'galois_orbit':
        return o
    if output == 'galois_orbit_numbers_list':
        return map(lambda x:x.number(),o)
    if output == 'galois_orbit_rep':
        return o[0]
    try:
        i = dirichlet_character_conrey_galois_orbits_reps(N).index(o[0].number())
        return x.modulus(),i
    except IndexError:
        raise ValueError,"Could not find Galois orbit of {0}".format(x)


def sage_character_to_conrey_character(x):
    r"""
    Return the Corney character corresponding to the sage character x
    """
    return dirichlet_group_conrey(x.modulus()).from_sage_character(x)

def sage_character_to_conrey_galois_orbit_number(x):
    r"""
    Return the number of the Galois orbit in DirichletGroup_conrey whjich contains x
    """
    y = sage_character_to_conrey_character(x)
    return conrey_to_conrey_(y,output='galois_orbit_number')


def conrey_character_number_to_conrey_galois_orbit_number(n,i):
    r"""
    Return the number of the galois orbit in DirichletGroup_conrey(n).galois_orbits()
    which contains the Dirichlet character DirichletGroup_conrey(n)[i]
    """
    return conrey_to_conrey_((n,i),number_format='character_number',output='galois_orbit_number')

    


## Representatives as numbers

@cached_function
def dirichlet_character_conrey_galois_orbit_numbers_from_character_number(n,xi):
    r"""
    Return the numbers of the characters in the galois orbit of the conrey character chi(n,xi)
    """
    conrey_to_conrey_((n,xi),number_format='character_number',output='galois_orbit_number')


@cached_function    
def dirichlet_character_conrey_galois_orbit_rep_from_character_number(n,xi):
    r"""
    Return the representative of the Galois orbit nr. i modulo N
    """
    return conrey_to_conrey_((n,xi),number_format='character_number',output='galois_orbit_rep')

### Conversions between the two types 

@cached_function    
def dirichlet_character_conrey_from_sage_character_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    x = sage_to_sage_((n,i),number_format='character_number',output='character')
    return dirichlet_character_conrey_from_sage_character(x)


@cached_function    
def dirichlet_character_conrey_from_sage_galois_orbit_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    x = sage_to_sage_((n,i),number_format='galois_orbit_number',output='character')
    return dirichlet_character_conrey_from_sage_character(x)
    

@cached_function    
def dirichlet_character_sage_from_conrey_character_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    x = conrey_to_conrey_((n,i),number_format='character_number',output='character')
    return x.sage_character()

@cached_function    
def dirichlet_character_sage_from_conrey_galois_orbit_number(n,i):
    r"""
    Return x in DirichletGroup_conrey(n) corresponding to y=DirichletGroup(n).galois_orbits()[i][0]
    """
    x = conrey_to_conrey_((n,i),number_format='galois_orbit_number',output='character')
    return x.sage_character()
  

def dirichlet_character_conrey_galois_orbit_rep_from_sage_character(x):
    r"""
    Return the representative of the Galois orbit nr. i modulo N
    """    
    y = sage_charcter_to_conrey_character(x)
    return conrey_to_conrey_(y,output='galois_orbit_rep')
    
@cached_function
def dirichlet_character_sage_galois_orbit_rep_from_sage_character(x):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N (using the sage representation of characters).
    """
    o = x.galois_orbit()
    for r in  dirichlet_character_sage_galois_orbits_reps(N):
        if r in o:
            return r
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
    

# @cached_function
# def dirichlet_character_conrey_used_in_computation(N,xi):
#     r"""
#       INPUTS:
#        - ```x```: A Conrey Dirichlet Character
    
#       Returns the number of Conrey Dirichlet Character ```c```,
#       such that ```c.sage_character()``` is the representative that
#       was used to compute the spaces of modular forms with character ```x```.

#       OUTPUT:
#        - int: the number of the corresponding Conrey Dirichlet Character.
#     """
#     if N == 1:
#         return int(1)
#     reps_sage = dirichlet_character_sage_galois_orbits_reps(N)
#     x = conrey_character_from_number(N,xi)
#     for c in x.galois_orbit():
#         if c.sage_character() in reps_sage:
#             return c.number()
    

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
    #base_number = dirichlet_character_conrey_used_in_computation(N,xi)
    base_number = conrey_to_conrey((N,xi),number_format='character_number',
                                   output='galois_orbit_rep')
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
