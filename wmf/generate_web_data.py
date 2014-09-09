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

from sage.all import parallel,dumps,Gamma1
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

from sage.all import dimension_new_cusp_forms
from compmf import character_conversions
import json 
#import bson

def my_dumps(s):
    r"""
    Use my own dump funciton to make it esier to switch between json and bson dumps.
    """
    return json.dumps(s)

def my_loads(s):
    return json.loads(s)
    
def generate_table(level_range=[1,500],weight_range=[2,12],chi_range=[],ncpus=1,host='localhost',only_new=True,port=int(37010)):
    r"""
    Generates a table of the available (computed) WebModFormSpaces.
    In addition we also add data for (level,weight,chi) with
    level in level_range, weight in weight_range and chi in chi_range
    

    We use the field level_max and weight_max in the database to indicate the endpoints of the range of systematically computed dimensions.
    """
    try: 
        D  = MongoMF(host,port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace_test'
    ## If we have an old table we load it.
    ## Do Gamma0 data first
    r0 = D._mongodb['webmodformspace_dimension'].find_one({'group':'gamma0'})
    id0 = None; id1 = None
    tbl0 = {}; tbl1={}
    level_max_in_db = 0;weight_max_in_db = 0
    if r0:
        if r0.get('data'):
            #print r0.get('data')
            tbl0 = my_loads(r0.get('data'))
            id0 = r0.get('_id')
            level_max_in_db = r0.get('level_max',0)
            weight_max_in_db = r0.get('weight_max',0)
    if level_max_in_db < level_range[1] or weight_max_in_db < weight_range[1]:
        for n in range(level_range[0],level_range[1]+1):
            tbl0[n]={}
            for k in range(weight_range[0],weight_range[1]+1):
                if tbl0.has_key(n):
                    if tbl0[n].has_key(k):
                        continue
                tbl0[n][k]=int(dimension_new_cusp_forms(n,k)),int(0)
    q = D._mongodb[webmodformspace].find({'character':int(1)}).sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)])
    for r in q:
        n = r['level']; k = r['weight']
        if not tbl0.has_key(n):
            tbl0[n] = {}
        if tbl0[n].has_key(k):
            d,t = tbl0[n][k]
        else:
            d = r['dimension_new_cusp_forms']
        tbl0[n][k] = (int(d),int(1))
    rec0 = {'group':'gamma0','data':my_dumps(tbl0),
            'level_max':int(level_range[1]),
            'weight_max':int(level_range[1])}
    # remove the old record
    if id0:
        D._mongodb['webmodformspace_dimension'].remove(id0)
    D._mongodb['webmodformspace_dimension'].insert(rec0)
    # Now compute the gamma1 data
    r1 = D._mongodb['webmodformspace_dimension'].find_one({'group':'gamma0'})
    level_max_in_db = 0; weight_max_in_db = 0
    if r1:
        if r1.get('_id'):
            tbl1 = my_loads(r1.get('data'))
            level_max_in_db = r1.get('level_max',0)
            weight_max_in_db = r1.get('weight_max',0)
    if level_max_in_db < level_range[1] or weight_max_in_db < weight_range[1]:
        for n in range(level_range[0],level_range[1]+1):
            tbl1[n]={}
            for k in range(weight_range[0],weight_range[1]+1):
                tbl1[n][k]={}
                ds = 0
                for x in character_conversions.dirichlet_character_conrey_galois_orbits_reps(n):
                    xi = x.number()
                    if (k % 2)==0 and not x.is_even():
                        d = 0
                    elif (k % 2)==1 and not x.is_odd():
                        d = 0
                    else:
                        d = dimension_new_cusp_forms(x.sage_character(),k)
                    tbl1[n][k][xi]=(int(d),int(0))
                    ds+=d
                tbl1[n][k][-1]=(int(ds),int(0))
    q = D._mongodb[webmodformspace].find({'character':int(1)})
    for r in q:
        n = r['level']; k = r['weight']
        i = r['character_orbit_rep']
        if not tbl1.has_key(n):
            tbl1[n] = {}
        if not tbl1[n].has_key(k):
            tbl1[n][k] = {}
        if  tbl1[n][k].has_key(i):
            d,t = tbl1[n][k][i]
        else:
            d = r['dimension_new_cusp_forms']
        tbl1[n][k][i] = (int(d),int(1))
        if not tbl1[n][k].has_key(-1):
            tot_dim = dimension_new_cusp_forms(Gamma1(n))
            d = sum(filter(lambda x:x[0],tbl1[n][k].values()))
            if d <> tot_dim:
                emf_log.warning("The sum of the computed dimensions does not ad up at N,k={0}".format((N,k)))
            tbl1[n][k][-1] = (int(tot_dim),int(1))
    ### Also add the total dimensions...
    rec1 = {'group':'gamma1','data':my_dumps(tbl1),
            'level_max':int(level_range[1]),
            'weight_max':int(level_range[1])}
    if id1:
        D._mongodb['webmodformspace_dimension'].remove(id1)
    D._mongodb['webmodformspace_dimension'].insert(rec1)
        
    return tbl0,tbl1

    
