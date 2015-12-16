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

import pymongo
import lmfdb
from lmfdb.modular_forms.elliptic_modular_forms import *
from lmfdb.number_fields.number_field import make_disc_key
from sage.all import ZZ,QQ

def Modf_changevar(f,NF,Bfacto=10^6):
    r"""
     Usage : f a hecke_orbit, NF=lmfdb.base.getDBConnection()['numberfields']['fields']
     Returns : [v2,label], where v2 is v expressed on a nice model of the coeff field, and label is the lmfdb label of the field (or '' if not in the database)
    """
    ZZx=ZZ['x']
    QQx=QQ['x']
    P=f.absolute_polynomial
    # If f is rational, nothing to do :)
    if f.is_rational:
        return [f.eigenvalues.v,u'1.1.1.1']
    # Is the coefficient field already identified ?
    Klabel=f.coefficient_field.lmfdb_label
    if Klabel:
        # It is, let us make the isomorphism explicit
        K=NF.find_one({'label':Klabel})
        Q=ZZx([ZZ(a) for a in K['coeffs'].split(',')])
        # Compute max order in gp
        maxram=max([ZZ(a) for a in K['ramps']])
        pkQ=gp.nfinit([Q,maxram+1])
        iso=QQx(str(gp.nfisisom(pkQ,P)[1]))
    else:
        # It is not, let us see if it exists in the DB
        plist=[p for p in range(200) if is_prime(p)]
        query={}
        # Lazy order
        pKP=gp.nfinit([P,Bfacto])
        # Sign
        [r1,r2]=pKP[2]
        query['signature']=str(r1)+','+str(r2)
        DpKP=pKP[3]
        # Is the lazy order maximal ?
        if len(gp.nfcertify(pKP)):
            # The lazy order is not maximal
            ur=[]
            ram=[]
            # Primes<200 known to be unramified
            for p in plist:
                if Mod(DpKP,p):
                    ur.append(str(p))
            # Lazy factorisation of the disc of the lazy order
            faD=gp.factor(DpKP,Bfacto)
            faD=str(faD).replace('[','').replace(']','').replace('Mat(','').replace(')','').split(';')
            # Primes known to be ramified
            for s in faD:
                p=s.split(',')[0].replace(' ','')
                if ZZ(p)<Bfacto:
                    ram.append(p)
            query['$nor'] = [{'ramps': x} for x in ur]
            query['ramps'] = {'$all': ram}
        else:
            # The lazy order is maximal :)
            # Query on disc
            s,D=make_disc_key(ZZ(DpKP))
            query['disc_sign']=s
            query['disc_abs_key']=D
            LK=NF.find(query)
            Klabel=''
            for K in LK:
                # Found a candidate in the nf DB, here is its defining polynomial
                Q=ZZx([ZZ(a) for a in K['coeffs'].split(',')])
                # Compute max order in gp
                maxram=max([ZZ(a) for a in K['ramps']])
                pkQ=gp.nfinit([Q,maxram+1])
                # Check for isomorphism
                iso=gp.nfisisom(pkQ,P)
                if iso:
                    iso=QQx(str(iso[1]))
                    Klabel=K['label']
                    break

    if Klabel=='':
        # Field nothing found, so we reduce the initial polynomial as we can
        [Q,iso]=gp.polredbest(P,1)
        Q=ZZx(str(Q))
        iso=QQx(str(gp.lift(iso)))

    # Finally, apply isomorphism
    KQ.<a>=NumberField(Q)
    iso=KQ(iso)
    v=f.eigenvalues.v
    newv=[l.lift()(iso) for l in v]
    return [newv,Klabel]
