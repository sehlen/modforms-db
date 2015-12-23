r"""

Utilities for modform computing classes.

"""

from sage.all import cached_function,AlphabeticStrings
from wmf import wmf_logger

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
from sage.all import ZZ,QQ,is_prime,gp,Mod,NumberField,vector

def Modf_changevar2(f,NF,Bfacto=10^6):
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
    KQ=NumberField(Q,name='a')
    iso=KQ(iso)
    v=f.eigenvalues.v
    newv=[l.lift()(iso) for l in v]
    return [newv,Klabel]



def Modf_changevar_f(f,NF=None,Bfacto=10^6):
    try:
        label = f.coefficiewnt_field.lmfdb_label
    except:
        label = ''
    E = f.eigenvalues.E
    v = f.eigenvalues.v
    return Modf_changevar_Ev(E,v,NF,Bfacto=Bfaco,Klabel=label)


def get_lmfdb_label(v,NF=None,Bfacto=10^6):
    return Modf_changevar_Ev(E=None,v=v,NF=NF,Bfacto=Bfacto,Klabel='',label_only=True)


def Modf_changevar_Ev(E=None,v=None,NF=None,Bfacto=10^6,Klabel='',label_only=False):
    r"""
    Usage : f a hecke_orbit, NF=lmfdb.base.getDBConnection()['numberfields']['fields']
    Returns : [v2,E2,Q,emb,label], where v2 and E2 are v and E
    expressed on a nice model of the coeff field, Q is the absolute
    defining polynomial of this model, emb is the embeddding of the
    generator of the cycltomic subfield (for Gamma1), and label is the
    lmfdb label of the field (or '' if not in the database)

    if label_only = True we only want to see if we can find the label,
    not apply the isomorphism.
    """
    import lmfdb
    if v is None:
        return [None,None,None,None,None]
    if NF is None:
        NF=lmfdb.base.getDBConnection()['numberfields']['fields']
    coefficient_field = v[0].parent()
    wmf_logger.debug("K={0}".format(coefficient_field))
    # If f is rational, nothing to do :)
    if coefficient_field.absolute_degree() == 1:
        if label_only:
            return u'1.1.1.1'
        return [E,v,QQ,'x',-1,u'1.1.1.1']
    P=coefficient_field.absolute_polynomial()
    ZZx=ZZ['x']
    QQx=QQ['x']
    wmf_logger.debug("P={0}".format(P))
    # Is the coefficient field already identified ?
    #Klabel=coefficient_field.lmfdb_label
    if Klabel != '':
        # It is, let us make the isomorphism explicit
        K=NF.find_one({'label':Klabel})
        Q=ZZx([ZZ(a) for a in K['coeffs'].split(',')])
        # Compute max order in gp
        pkQ=gp.nfinit([Q,[ZZ(p) for p in K['ramps']]])
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
        #print "query=",query
        Klabel=''
        for K in LK:
            #print "K=",K
            # Found a candidate in the nf DB, here is its defining polynomial
            Q=ZZx([ZZ(a) for a in K['coeffs'].split(',')])
            # Compute max order in gp
            pkQ=gp.nfinit([Q,[ZZ(p) for p in K['ramps']]])
            # Check for isomorphism
            iso=gp.nfisisom(pkQ,P)
            if iso:
                iso=QQx(str(iso[1]))
                Klabel=K['label']
                break

    if label_only:
        return Klabel
            
    if Klabel=='':
        # Field not found, so we reduce the initial polynomial as we can
        [Q,iso]=gp.polredbest(P,1)
        Q=ZZx(str(Q))
        pkQ=gp.nfinit([Q,Bfacto])
        iso=QQx(str(gp.lift(iso)))
    # Now we have the model we want for the absolute field.
    # We now want the explicit embedding of the cyclotomic field, the
    # relative polynomial for thi new field, and the relative version of the isomorphism
    #E=f.eigenvalues.E
    #v=f.eigenvalues.v
    KQ=NumberField(Q,name='a')
    a = KQ.gen()
    Kcyc=v[0].parent().base_ring()
    if Kcyc.degree()>1:
        polcyc=Kcyc.defining_polynomial()
        relP=v[0].parent().defining_polynomial()
        emb=QQx(str(gp.nfisincl(polcyc,pkQ)[1]))(a)
        Krel=Kcyc.extension(relP,name='a')
        a = Krel.gen()
        osi=gp.lift(gp.modreverse(gp.Mod(iso,Q)))
        osi=QQx(str(osi))
        relQ=osi(a).charpoly()
        R=Kcyc.extension(relQ,name='a')
        a = R.gen()
        relIso=iso(a)
        newv=vector([l.lift()(relIso) for l in v])
        if E.base_ring() != Kcyc:
            E=E.apply_map(lambda x:x[0])
        return [E,newv,Q,emb,Klabel]

    # Trivial cycltomic field case
    KQ = NumberField(Q,name='a')
    iso=KQ(iso)
    newv=vector([l.lift()(iso) for l in v])
    Enew=E.apply_map(lambda x: x.lift()(iso))
    return [Enew,newv,Q,-1,Klabel]
