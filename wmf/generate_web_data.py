# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <fredrik314@gmail.com>,
#  Stephan Ehlen <>
#  Distributed under the terms of the GNU General Public License (GPL)
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
r""" Class for newforms in format which can be presented on the web easily


AUTHORS:

 - Fredrik Stroemberg
 - Stephan Ehlen
"""

import pymongo

from sage.all import parallel
from wmf import wmf_logger,WebNewForm_computing,WebModFormSpace_computing
from compmf import MongoMF

def generate_web_modform_spaces(level_range=[],weight_range=[],chi_range=[],ncpus=1,recompute=False,host='localhost',port=int(37010)):
    r"""
    Compute and insert objects of type WebModFormSpace with levels in the given range.

    NOTE: We only compute forms which have an entry in the mongodb so you need to use the CompMF class first to generate these.

    """
    try: 
        D  = MongoMF(host,port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    args = []; s={}
    if level_range <> []:
        if len(level_range)==1:
            s['N']=int(level_range[0])
        else:
            s['N']={"$gt":int(level_range[0]-1),"$lt":int(level_range[-1]+1)}
    if weight_range <> []:
        if len(weight_range)==1:
            s['k']=int(weight_range[0])
        else:
            s['k']={"$gt":int(weight_range[0]-1),"$lt":int(weight_range[-1]+1)}
    if chi_range <>[]:
        if len(chi_range)==1:
            s['chi']=int(chi_range[0])
        else:
            s['chi']={"$gt":int(chi_range[0]-1),"$lt":int(chi_range[-1]+1)}
    s['complete']=int(3)
    q = D._modular_symbols.find(s).sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING)])
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace_test'
    for r in q:
        N = r['N']; k=r['k']; chi=r['cchi']
        if recompute is False:
            if D._mongodb['webmodformspace'].find({'level':int(N),'weight':int(k),'character':int(chi)}).count()>0:
                continue
        args.append((N,k,chi))
    print "s=",s
    print "args=",args
    print "ncpus=",ncpus
    if ncpus>=32:
        l = generate_one_webmodform_space32(args)
    elif ncpus>=16:
        l = generate_one_webmodform_space16(args)        
    elif ncpus>=8:
        l = generate_one_webmodform_space8(args)
    elif ncpus>=4:
        l = generate_one_webmodform_space4(args)
    else:
        l =  generate_one_webmodform_space1_par(args)
    return list(l)
    
@parallel(ncpus=32)
def generate_one_webmodform_space32(level,weight,chi,**kwds):
    return generate_one_webmodform_space1(level,weight,chi,**kwds)
@parallel(ncpus=16)
def generate_one_webmodform_space16(level,weight,chi,**kwds):
    return generate_one_webmodform_space1(level,weight,chi,**kwds)
@parallel(ncpus=8)
def generate_one_webmodform_space8(level,weight,chi,**kwds):
    return generate_one_webmodform_space1(level,weight,chi,**kwds)
@parallel(ncpus=4)
def generate_one_webmodform_space4(level,weight,chi,**kwds):
    return generate_one_webmodform_space1(level,weight,chi,**kwds)
@parallel(ncpus=1)
def generate_one_webmodform_space1_par(level,weight,chi,**kwds):
    return generate_one_webmodform_space1(level,weight,chi,**kwds)

    
def generate_one_webmodform_space1(level,weight,chi):    
    r"""
    Generates one modform space.

    """
#    print "generate:",level,weight,chi
    M = WebModFormSpace_computing(level,weight,chi)


def generate_table(level_range=[1,500],weight_range=[2,12],chi_range=[],ncpus=1,host='localhost',only_new=True,port=int(37010)):
    r"""
    Generates a table of the available (computed) WebModFormSpaces.
    In addition we also add data for (level,weight,chi) with
    level in level_range, weight in weight_range and chi in chi_range
    
    """
    try: 
        D  = MongoMF(host,port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
#    M = WebNewForm(1,12,1)
#    q = D.['
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace_test'
    ## Do Gamma0 data first
    tbl0 = {}
    q = D._mongodb[webmodformspace].find({'character':int(1)})
    for n in level_range:
        tbl0[n]={}
        for k in weight_range:
            tbl0[n][k]=0,dimension_new_cusp_forms(n,k)
    for r in q:
        n = r['level']; k = r['weight']
        if n > level_range(1):
            tbl0['level'][n] = {}
        if tbl0[n].has_key(k):
            t,d = tbl0[n][k]
        else:
            d = r['dimension_new_cusp_forms']
        tbl0[n][k] = 1,d
        

    
