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
import gridfs
import bson
from sage.all import parallel,dumps,Gamma1,QQ,prime_pi,RR,deepcopy,nth_prime
from wmf import wmf_logger,WebNewForm_computing,WebModFormSpace_computing
from compmf import MongoMF,MongoMF,data_record_checked_and_complete,CompMF,CheckingDB
from compmf.utils import multiply_mat_vec,convert_matrix_to_extension_fld
from sage.misc.cachefunc import cached_function
from lmfdb.modular_forms.elliptic_modular_forms import emf_version

def generate_web_modform_spaces(level_range=[],weight_range=[],chi_range=[],ncpus=1,recompute=False,host='localhost',port=int(37010),user=None,password=None):
    r"""
    Compute and insert objects of type WebModFormSpace with levels in the given range.

    NOTE: We only compute forms which have an entry in the mongodb so you need to use the MongoMF class first to generate these.

    """
    try:
        D  = MongoMF(host=host,port=port,user=user,password=password)
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
            s['cchi']=int(chi_range[0])
        else:
            s['cchi']={"$gt":int(chi_range[0]-1),"$lt":int(chi_range[-1]+1)}
    s['complete']={"$gt":int(data_record_checked_and_complete-1)}
    wmf_logger.debug("checking to level:{0}".format(data_record_checked_and_complete))
    q = D._modular_symbols.find(s).sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING)])
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace'
    for r in q:
        N = r['N']; k=r['k']; cchi=r['cchi']
        if recompute is False:
            if D._mongodb['webmodformspace'].find({'level':int(N),'weight':int(k),'character':int(cchi),'version':emf_version}).count()>0:
                continue
        args.append((N,k,cchi))
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
        l = []
        for N,k,cchi in args:
            l.append(generate_one_webmodform_space1_par(N,k,cchi))
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
#@parallel(ncpus=1)
def generate_one_webmodform_space1_par(level,weight,chi,**kwds):
    return generate_one_webmodform_space1(level,weight,chi,**kwds)


def generate_one_webmodform_space1(level,weight,cchi,host='localhost',port=int(37010),recompute=True):
    r"""
    Generates one modform space.

    """
    #    print "generate:",level,weight,chi
    D = MongoMF(host=host,port=port)
    cid = D.register_computation(level=level,weight=weight,cchi=cchi,typec='wmf')
    M = WebModFormSpace_computing(level,weight,cchi,recompute=recompute)
    M.save_to_db()
    D.register_computation_closed(cid)


def generate_web_eigenvalues(level_range=[],weight_range=[],chi_range=[],ncpus=1,recompute=False,host='localhost',port=int(37010),user=None,password=None):
    r"""
    Compute and insert objects of type WebModFormSpace with levels in the given range.

    NOTE: We only compute forms which have an entry in the mongodb so you need to use the MongoMF class first to generate these.

    """
    try:
        D  = MongoMF(host=host,port=port,user=user,password=password)
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
            s['cchi']=int(chi_range[0])
        else:
            s['cchi']={"$gt":int(chi_range[0]-1),"$lt":int(chi_range[-1]+1)}
    for r in D._aps.find(s):
        label = r['hecke_orbit_label']
        prec = r['prec']
        if D._mongodb['webeigenvalues.files'].find({'hecke_orbit_label':label,'prec':prec}):
            continue
        args.append((label,prec))
    print "s=",s
    print "args=",args
    print "ncpus=",ncpus
    l = generate_one_webeigenvalue32(args)
    return list(l)

@parallel(32)
def generate_one_webeigenvalue32(label,prec):
    E = WebEigenvalues(label,prec=prec)
    return E.save_to_db()
    
def web_modformspace_collection(host='localhost',port=int(37010)):
    try:
        D = MongoMF(host=host,port=port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    try:
        col = WebModFormSpace_computing._collection_name
    except AttributeError:
        col = 'webmodformspace'
    return D._mongodb[col]

def web_newform_collection(host='localhost',port=int(37010)):
    try:
        D = MongoMF(host=host,port=port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    try:
        col = WebNewForm_computing._collection_name
    except AttributeError:
        col = 'webnewforms'
    return D._mongodb[col]

## Create the indices we want on the collections

from pymongo import ASCENDING, DESCENDING

def create_index(host='localhost',port=int(37010),only=None):
    #c = web_newform_collection(host,port)
    C = MongoMF(host=host,port=port)
    D = C._mongodb
    collections_indices = [
        {'name':'webnewforms', 'index':[
                ('level',pymongo.ASCENDING),
                ('weight',pymongo.ASCENDING),
                ('chi',pymongo.ASCENDING)],
         'unique':False},
        {'name':'webnewforms', 'index':[
            ('hecke_orbit_label',pymongo.ASCENDING),
            ('version',pymongo.ASCENDING)],
         'unique':True},
        {'name': 'webnewforms.files', 'index': [
            ('hecke_orbit_label',pymongo.ASCENDING),
            ('version',pymongo.ASCENDING)],
         'unique':True},
        {'name': 'webmodformspace' ,'index':[
                ('level',pymongo.ASCENDING),
                ('weight',pymongo.ASCENDING),
                ('chi',pymongo.ASCENDING)],
         'unique':False},
        {'name': 'webmodformspace', 'index':[
                ('space_label',pymongo.ASCENDING),
                ('version',pymongo.ASCENDING)],
         'unique':True},
        {'name':  'webmodformspace.files', 'index':[
            ('space_label',pymongo.ASCENDING),
            ('version',pymongo.ASCENDING)],
         'unique':True},
        {'name':  'webchar', 'index': [
                ('modulus',pymongo.ASCENDING),
                ('version',pymongo.ASCENDING),
                ('number',ASCENDING)],
         'unique':True},
        #       {'name':  'webchar', 'index': [
        #            ('label',ASCENDING)],
        #        'unique':True},
        {'name' : 'webeigenvalues', 'index':[
                ('hecke_orbit_label',ASCENDING),
            ('version',pymongo.ASCENDING)],
         'unique':True}
        ]
    for r in collections_indices:
        if not only is None:
            if r['name'] <> only:
                continue
        if r['name'] in D.collection_names():
            D[r['name']].create_index(r['index'],unique=r['unique'])
            print "Created indeex r={0}".format(r)





from sage.all import dimension_new_cusp_forms
from compmf import character_conversions
from compmf.character_conversions import  dirichlet_group_conrey_galois_orbits,conrey_character_from_number,dirichlet_group_conrey_galois_orbits_numbers,dirichlet_character_sage_from_conrey_character_number,sage_character_to_conrey_character
import json
#import bson

def my_dumps(s):
    r"""
    Use my own dump funciton to make it esier to switch between json and bson dumps.
    """
    return json.dumps(s)
#    return dumps(s)

def my_loads(s):
    #    return json.loads(s)
    return json.loads(s)

def generate_dimension_tables(level_range=[1,500],weight_range=[2,12],chi_range=[],ncpus=1,only_new=True,host='localhost',port=int(37010)):
    r"""
    Generates a table of the available (computed) WebModFormSpaces.
    In addition we also add data for (level,weight,chi) with
    level in level_range, weight in weight_range and chi in chi_range


    We use the field level_max and weight_max in the database to indicate the endpoints of the range of systematically computed dimensions.
    """
    try:
        D = MongoMF(host=host,port=port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace'
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
        #    if level_max_in_db < level_range[1] or weight_max_in_db < weight_range[1]:
        #print tbl0.keys()
        for n in range(level_range[0],level_range[1]+1):
            ns = str(n)
            if not tbl0.has_key(ns):
                tbl0[ns]={}
                for k in range(weight_range[0],weight_range[1]+1):
                    ks = str(k)
                    if tbl0[ns].has_key(ks):
                        continue
                    tbl0[ns][ks]=int(dimension_new_cusp_forms(n,k)),int(0)
    q = D._mongodb[webmodformspace].find({'character':int(1)}).sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)])
    for r in q:
        n = str(r['level']); k = str(r['weight'])
        if not tbl0.has_key(n):
            tbl0[n] = {}
        if tbl0[n].has_key(k):
            d,t = tbl0[n][k]
        else:
            d = r['dimension_new_cusp_forms']
        tbl0[n][k] = (int(d),int(1))
    rec0 = {'group':'gamma0','data':my_dumps(tbl0),
            'level_max':int(level_range[1]),
            'weight_max':int(level_range[1]),
            'date':bson.datetime.datetime.now()}
    # remove the old record
    if id0:
        D._mongodb['webmodformspace_dimension'].remove(id0)
    D._mongodb['webmodformspace_dimension'].insert(rec0)
    print "inserted new table!"
    # Now compute the gamma1 data
    r1 = D._mongodb['webmodformspace_dimension'].find_one({'group':'gamma1'})
    level_max_in_db = 0; weight_max_in_db = 0

    if r1:
        if r1.get('_id'):
            id1 = r1.get('_id')
            tbl1 = my_loads(r1.get('data'))
            level_max_in_db = r1.get('level_max',0)
            weight_max_in_db = r1.get('weight_max',0)
    wmf_logger.debug("level_max_in_db={0}".format(level_max_in_db))
    wmf_logger.debug("level_range={0}".format(level_range))
    #    if level_max_in_db < level_range[1] or weight_max_in_db < weight_range[1]:
    for n in range(level_range[0],level_range[1]+1):
        ns = str(n)
        if not tbl1.has_key(ns):
            tbl1[ns]={}
        #wmf_logger.debug("n={0}".format(n))
        for k in range(weight_range[0],weight_range[1]+1):
            ks = str(k)
            if tbl1[ns].has_key(ks):
                continue
            tbl1[ns][ks]={}
            ds = 0
            for x in character_conversions.dirichlet_character_conrey_galois_orbits_reps(n,format='character'):
                xi = x.number()
                if (k % 2)==0 and not x.is_even():
                    d = 0
                elif (k % 2)==1 and not x.is_odd():
                    d = 0
                else:
                    d = dimension_new_cusp_forms(x.sage_character(),k)
                tbl1[ns][ks][str(xi)]=(int(d),int(0))
                ds+=d
            tbl1[ns][ks][str(-1)]=(int(ds),int(0))
#            if n==190:
#                wmf_logger.warning("tbl[190][{0}]={1}".format(k,tbl1[n][k]))

    q = D._mongodb[webmodformspace].find().sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)])
    for r in q:
        n = r['level']; k = r['weight']
        i = r['character_orbit_rep']
        ns = str(n); ks=str(k); ist=str(i)
        if not tbl1.has_key(ns):
            tbl1[ns] = {}
        if not tbl1[ns].has_key(ks):
            tbl1[ns][ks] = {}
        if  tbl1[ns][ks].has_key(ist):
            d,t = tbl1[ns][ks][ist]
        else:
            d = r['dimension_new_cusp_forms']
        tbl1[ns][ks][ist] = (int(d),int(1))
        if not tbl1[ns][ks].has_key(str(-1)):
            tot_dim = dimension_new_cusp_forms(Gamma1(n),k)
            t = 0
            if n <= level_range[1] and k<=weight_range[1]:
                d = sum(map(lambda x:x[0],tbl1[ns][ks].values()))
                if d <> tot_dim:
                    wmf_logger.warning("The sum of the computed dimensions does not add up at N,k={0} d={1}, tot_dim={2}".format((n,k),d,tot_dim))
                    #                    wmf_logger.warning("tbl1[n][k]={0}".format(tbl1[n][k]))

                else:
                    t = 1
            tbl1[ns][ks][str(-1)] = (int(tot_dim),int(t))
    ### Also add the total dimensions...
    rec1 = {'group':'gamma1','data':my_dumps(tbl1),
            'level_max':int(level_range[1]),
            'date':bson.datetime.datetime.now(),
            'weight_max':int(level_range[1])}
    if id1:
        D._mongodb['webmodformspace_dimension'].remove(id1)
    D._mongodb['webmodformspace_dimension'].insert(rec1)

    return tbl0,tbl1


def update_dimension_tables(host='localhost',port=int(37010)):
    r"""
    Update the tabel by adding all known spaces.


    We use the field level_max and weight_max in the database to indicate the endpoints of the range of systematically computed dimensions.
    """
    try:
        D  = MongoMF(host=host,port=port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace'
    ## If we have an old table we load it.
    ## Do Gamma0 data first
    r0 = D._mongodb['webmodformspace_dimension'].find_one({'group':'gamma0'})
    r1 = D._mongodb['webmodformspace_dimension'].find_one({'group':'gamma1'})
    id0 = None; id1 = None
    tbl0 = {}; tbl1={}
    level_max_in_db = 0;weight_max_in_db = 0
    if r0:
        if r0.get('data'):
            tbl0 = my_loads(r0.get('data'))
            id0 = r0.get('_id')
    if r1:
        if r1.get('data'):
            tbl1 = my_loads(r1.get('data'))
            id1 = r1.get('_id')
    #q = D._mongodb[webmodformspace].find().sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)])
    q = D._modular_symbols.find().sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING)])
    wmf_logger.debug("Update dimension tables!")
    for r in q:
        #n = r['level']; k = r['weight']; i = r['character_orbit_rep']
        n = r['N']; k = r['k']; i = min(r['character_galois_orbit'])
        #if n==1 and k==76:
        #    wmf_logger.debug("r[1][76]={0}\n i={1}".format(r,i))
        n = str(n); k=str(k); i=str(i)
        if i == '1':
            if not tbl0.has_key(n):
                tbl0[n] = {}
            d = r['dimn']
            if D._mongodb['webmodformspace'].find({'space_orbit_label':r['space_orbit_label']}).count()>0: #if  D._mongodb[webmodformspace].find({'level':int(n),'weight':int(k)}):
                tbl0[n][k] = (int(d),int(1))
            else:
                if n=='8':
                    wmf_logger.debug("name {0} not in db!".format(r['space_orbit_label']))
                tbl0[n][k] = (int(d),int(0))
            #if n=='1' and k=='76':
            #    print tbl0[n][k]
        if not tbl1.has_key(n):
            tbl1[n] = {}
        if not tbl1[n].has_key(k):
            tbl1[n][k] = {}
        if tbl1[n][k].has_key(i):
            d,t = tbl1[n][k][i]
        else:
            d = r['dimn']
            #d = r['dimension_new_cusp_forms']
        if D._mongodb[webmodformspace].find({'galois_orbit_name':r['space_orbit_label']}).count()>0:
            tbl1[n][k][i] = (int(d),int(1))
        else:
            x = conrey_character_from_number(int(n),int(i))
            d = dimension_new_cusp_forms(x.sage_character(),int(k))
            tbl1[n][k][i] = (int(d),int(0))
        #if not tbl1[n][k].has_key("-1"):
        tbl1[n][k]["-1"] = (int(dimension_new_cusp_forms(Gamma1(int(n)),int(k))),int(1))
    for n in tbl1.keys():
        orbits = dirichlet_group_conrey_galois_orbits_numbers(int(n))
        if int(n)==17:
            wmf_logger.debug("orbits={0} ".format(orbits))
        for k in tbl1[n].keys():
            dtot = 0
            for i in tbl1[n][k].keys():
                mul = 0
                if int(n)==17 and int(k)==10:
                    wmf_logger.debug(" i={0}".format(int(i)))
                if int(i) > 0:

                    # Remember that these are for Galois conjugacy classes
                    for o in orbits:
                        if int(n)==17 and int(k)==10:
                            wmf_logger.debug("o={0} i={1}".format(o,int(i)))
                        if int(i) in o:
                            mul = len(o)
                            break
                    dtot+=int(tbl1[n][k][i][0])*mul
            d = tbl1[n][k].get("-1",(0,0))[0]
            if int(n)==17 and int(k)==10:
                wmf_logger.debug("dtot={0} d={1}".format(dtot,d))
            if str(d)==str(dtot):
                tbl1[n][k]["-1"]=(int(dtot),int(1))
            else:
                tbl1[n][k]["-1"]=(int(dtot),int(0))
                if n == "10" and k=="4":
                    print "sum of dims = {0} and true dim = {1}".format(dtot,d)
    t0 = my_dumps(tbl0)
    t1 = my_dumps(tbl1)
    d0 = D._mongodb['webmodformspace_dimension'].update({'_id':id0},{"$set": {'data':t0,'date':bson.datetime.datetime.now()}},upsert=True)
    d1 = D._mongodb['webmodformspace_dimension'].update({'_id':id1},{"$set": {'data':t1,'date':bson.datetime.datetime.now()}},upsert=True)
    print "update ",id0,id1
    print "res=",d0,d1
    return tbl0,tbl1,t0,t1


def drop_webmodform_data(host='localhost',port=int(37010)):
    r"""
    Drop all collections related to webnewforms and webmodforms
    """
    try:
        D  = MongoMF(host=host,port=port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    try:
        webmodformspace = WebModFormSpace_computing._collection_name
        webnewforms = WebNewForm_computing._collection_name
    except AttributeError:
        webmodformspace = 'webmodformspace'
        webnewforms = 'webnewforms'

    for col in [webmodformspace,webnewforms]:
        D._mongodb.drop_collection(col)
        D._mongodb.drop_collection(col+'.files')
        D._mongodb.drop_collection(col+'.chunks')

    D._mongodb.drop_collection('webmodformspace_dimension')

@cached_function
def dimension_from_db(level,weight,chi=None,group='gamma0'):
    import json
    try:
        D  = MongoMF(host=host,port=port)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    db = D._mongodb['webmodformspace_dimension']
    q = db.find_one({'group':group})
    dim_table = {}
    if q:
        dim_table = q.get('data',{})
        dim_table = json.loads(dim_table)
    if group=='gamma0' and chi<>None:
        d,t = dim_table.get(str(level),{}).get(str(weight),{}).get(str(chi),[-1,0])
        return  d,t
    elif chi is None:
        d,t = dim_table.get(str(level),{}).get(str(weight),[-1,0])
        return  d,t
    elif chi == 'all':
        res = {level: {weight:{}}}
        dtable = dim_table.get(str(level),{}).get(str(weight),{})
        for i in dtable.keys():
            res[level][weight][int(i)] = dtable[i]
        return res

def web_modformspace_in_db(host='localhost',port=int(37010)):
    res = []
    col = web_newform_collection(host,port)
    q = col.find().sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)])
    for r in q:
        res.append((r['level'],r['weight'],r['character']))
    return res



### Fix the addiiton of names.


from wmf.web_modform_space_computing import orbit_label
def add_orbit_labels_to_aps(host='localhost',port=int(37010)):
    D  = MongoMF(host=host,port=port)
    for r in D._aps.find({"hecke_orbit_label":{"$exists":False}}).sort('N',int(1)):
    #{'name':{"$exists":False}}):
        N=r['N']
        k=r['k']
        fid=r['_id']
        chi=r['chi']
        cchi=r.get('cchi')
        #if cchi==None:
        #    cchi = dirichlet_character_conrey_from_sage_character_number(N,chi)
        #    DB._aps.update({'_id':fid},{"$set":{'cchi':cchi.number()}})
        d=r['newform']
        label = orbit_label(d)
        name = '{0}.{1}.{2}{3}'.format(N,k,cchi,label)
        D._aps.update({'_id':fid},{"$set":{'hecke_orbit_label':name}})
        D._aps.update({'_id':fid},{"$unset":{'name':""}})


def add_hecke_orbits(host='localhost',port=int(37010)):
    import compmf
    from utils import orbit_label
    D  = MongoMF(host=host,port=port)
    spaces = D._mongodb.webmodformspace.distinct('space_label')
    for label in spaces:
        N,k,i = map(int,label.split("."))
        M = WebModFormSpace_computing(N,k,i)
        if M.hecke_orbits == {}:
            for d in range(M.dimension_new_cusp_forms):
                flabel = orbit_label(d)
                F = WebNewForm(N,k,i,flabel)
                M.hecke_orbits[flabel]=F
                M.save_to_db()
            print "Fixed {0}".format((N,k,i))

def remove_duplicates(D,label=None):
    col = D._mongodb['webchar']
    if label is None:
        labels = col.distinct("label")
    else:
        labels = [label]
    for lab in labels:
        q = col.find({'label':lab})
        if q.count()>1:
            print "Duplicates at {0}: num={1}".format(lab,q.count())
            l = list(q)
            nl = len(l)
            kept=0
            for i in range(nl):
                remove = True
                r = l[i]
                if r['label']<>[] and kept==0:
                    # keep
                    kept  = 1
                    remove = False
                if i == nl-1 and kept == 0: # at the end we keep at least one...
                    remove = False
                if remove:
                    col.remove({'_id':r['_id']})
                    print "removed {0}:{1}".format(lab,r['_id'])


def fix_galois_orbit_labels(D):
    import compmf
    from lmfdb.modular_forms.elliptic_modular_forms.backend.emf_utils import space_label
    import pymongo
    col = D._mongodb['webmodformspace.files']
    for q in col.find(): #{"version":float(1.2)}):
        gal_old = q.get('galois_orbit_name')
        space_label_old = q.get('space_label')
        if q.has_key('weight'):
            k = q['weight']; N=q['level']; i = q['character']
        else:
            N,k,i = space_label_old.split(".")
            N=int(N); k=int(k); i=int(i)
        space_label_new = space_label(level=N, weight=k, character=i)
        N,o=compmf.character_conversions.conrey_character_number_to_conrey_galois_orbit_number(N,i)
        gal_new = '{0}.{1}.{2}'.format(N,k,o)
        if gal_old <> gal_new:
            fid = q['_id']
            if pymongo.version_tuple[0]>=3:
                col.update_one({'_id':fid},{"$set":{'galois_orbit_name':gal_new}})
            else:
                col.update({'_id':fid},{"$set":{'galois_orbit_name':gal_new}})
            print "Want to fix label {0} -> {1}".format(gal_old,gal_new)
        if space_label_old <> space_label_new:
            print "Want to fix space label {0} -> {1}".format(space_label_old,space_label_new)
            if pymongo.version_tuple[0]>=3:
                col.update_one({'_id':fid},{"$set":{'space_label':space_label_new}})
            else:
                col.update({'_id':fid},{"$set":{'space_label':space_label_new}})

def fix_galois_orbit_labels_files(D,verbose=0):
    import compmf
    from lmfdb.modular_forms.elliptic_modular_forms.backend.emf_utils import space_label
    import pymongo
    import gridfs
    from sage.all import loads
    col = D._mongodb['webmodformspace.files']
    fs =  gridfs.GridFS(D._mongodb,'webmodformspace')
    for q in col.find():
        if verbose>0:
            print "Checking: ",q['space_label'],q['galois_orbit_name']
        label = q['space_label']
        N,k,i = label.split(".")
        N=int(N); k=int(k); i=int(i)

        name_new = q.get('galois_orbit_name')
        fid = q['_id']
        fd = loads(fs.get(fid).read())
        name_old = fd['galois_orbit_name']
        if name_new <> name_old:
            print "Want to fix label {0} -> {1}".format(name_old,name_new)
            fd['galois_orbit_name']=name_new


def remove_gridfs_duplicates(D,label_in=None):
    import gridfs
    fs = gridfs.GridFS(D._mongodb,collection='webmodformspace')
    col = D._mongodb['webmodformspace.files']
    if label_in is None:
        labs = col.distinct('space_label')
    else:
        labs = [label_in]
    for label in labs:
        q = col.find({'space_label':label}).sort("uploadDate", 1)
        n = q.count()
        if n>1:
            print "duplicate:",label
            i = 0
            for f in q:
                fid = f['_id']
                if i < n-1:
                    fs.delete(fid)
                    print "removed f=",f
                else:
                    print "keep =",f
                i+=1
#            fsq = fs.find({'hecke_orbit_label':label})
#            if fsq.count()>1:


def recompute_existing(D,ncpus=1,llim=10):
    llim = int(llim)
    args = []
    #    for r in D._mongodb['webnewforms'].find({'version':{"$lt":float(1.3)}}).limit(llim):
    for r in D._mongodb['webmodformspace'].find({'creation_date':{"$exists":False}}):
        level = r['level']; weight=r['weight']; character = r['character']
        args.append((level,weight,character))
    print "Recomputing {0} spaces!".format(len(args))
    if ncpus>=32:
        l = generate_one_webmodform_space32(args,recompute=True)
    elif ncpus>=16:
        l = generate_one_webmodform_space16(args,recompute=True)
    elif ncpus>=8:
        l = generate_one_webmodform_space8(args,recompute=True)
    elif ncpus>=4:
        l = generate_one_webmodform_space4(args,recompute=True)
    else:
        l =  generate_one_webmodform_space1_par(args,recompute=True)
    return list(l)

def remove_newform_with_label(D,hecke_orbit_label):
    import gridfs
    coll = D._mongodb['webnewforms']
    key = {'hecke_orbit_label':hecke_orbit_label}
    if all:
        r = coll.delete_many(key) # delete meta records
    else:
        r = coll.delete_one(key) # delete meta records
    if r.deleted_count == 0:
        wmf_logger.debug("There was no meta record present matching {0}".format(key))
    file_collection = D._mongodb['webnewforms.files']
    r = file_collection.find_one(key)
    if r is None:
        raise IndexError("Record does not exist: {0}".format(key))
    fid = r['_id']
    fs = gridfs.GridFS(D._mongodb,'webnewforms')
    fs.delete(fid)


def recompute_newforms(D):
    q = D._mongodb['webnewforms'].find({"$where": "this.q_expansion.length < 10"}).distinct('hecke_orbit_label')
    for label in q:
        remove_newform_with_label(D,label)
    args = []
    for label in q:
        # remove the alphabetic label at the end
        la = "".join([x for x in label if x.isalpha()])
        label = label.replace(la,"")
        N,k,i = label.split(".")
        N = int(N); k = int(k); i = int(i)
        args.append((N,k,i))
    l =  generate_one_webmodform_space32(args,recompute=True)
    return list(l)


def recheck_and_compute_1(D):
    from sage.all import loads,RR
    import gridfs
    coll = D._mongodb['webnewforms']
    file_collection = D._mongodb['webnewforms.files']
    labels_to_remove = []
    for r in coll.find():
        label = r['hecke_orbit_label']
        rr = file_collection.find_one({'hecke_orbit_label':label})
        fid = rr['_id']
        fs = gridfs.GridFS(D._mongodb,'webnewforms')
        frec = loads(fs.get(fid).read())
        emb = frec.get('_embeddings')
        if emb is None:
            continue
        k =RR(r['weight'])
        if not emb.has_key('values'):
            wmf_logger.debug("Record {0} embedding {1} has no values".format(label,emb))
        c2 = emb['values'][2][0]/2.0**((k-1)/2.0)
        if abs(c2)>2:
            wmf_logger.debug("Record {0} has too large c(2):{1}".format(label,c2))
            labels_to_remove.append(label)
    args = []
    for label in labels_to_remove:
        remove_newform_with_label(D,label)
        la = "".join([x for x in label if x.isalpha()])
        label = label.replace(la,"")
        N,k,i = label.split(".")
        N = int(N); k = int(k); i = int(i)
        args.append((N,k,i))
    l =  generate_one_webmodform_space32(args,recompute=True)
    return len(list(l))

def add_zetas(D):
    from sage.all import loads,RR
    import gridfs
    coll = D._mongodb['webmodformspace']
    #file_collection = D._mongodb['webmodformspace']
    args = []
    for r in coll.find({"zeta_orders":{"$exists":False}}): #,'character':{"$ne":int(1)}}):
        args.append((r['level'],r['weight'],r['character']))
    l =  add_zeta_parallel(args)
    return len(list(l))
    #M=webmodformspace_computing(r['level'],r['weight'],r['character'],compute=False)
    #    M.get_zetas()
    #    M.save_to_db()
@parallel(ncpus=32)
def add_zeta_parallel(level,weight,cchi,host='localhost',port=int(37010)):
    r"""
    Generates one modform space.

    """
    #    print "generate:",level,weight,chi
    D = MongoMF(host=host,port=port)
    M = WebModFormSpace_computing(level,weight,cchi,recompute=False)
    M.get_zetas()
    M.save_to_db()


def fix_orbit_labels(D):
    from compmf.character_conversions import conrey_character_number_to_conrey_galois_orbit_number
    for r in D._modular_symbols.find():
        space_orbit_label = r.get('space_orbit_label')
        print space_orbit_label
        #if not space_orbit_label is None:
        #    l =space_orbit_label.split(".")
        #    if  l[2].isdigit():
        #        continue
        #    new_label = "{0}.{1}.{2}".format(l[0],l[1],l[2][1])
        #else:
        N = r.get('N'); k=r.get('k'); ci=r.get('cchi')
        if N is None or k is None or ci is None:
            wmf_logger.critical("r={0}".format(r))
            D._modular_symbols.remove({'_id':r['_id']})
            continue
        on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)[1]
        new_label="{0}.{1}.{2}".format(N,k,on)
        D._modular_symbols.update({'_id':r['_id']},{"$set":{'space_orbit_label':new_label}})

def fix_orbit_labels_2(D):
    from compmf.character_conversions import conrey_character_number_to_conrey_galois_orbit_number
    for r in D._mongodb['webmodformspace'].find(): #{'space_orbit_label':{"$exists":False}}):
        N = r.get('level'); k=r.get('weight'); ci=r.get('character')
        if N is None or k is None or ci is None:
            wmf_logger.critical("r={0}".format(r))
            continue
        on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)[1]
        new_label="{0}.{1}.{2}".format(N,k,on)
        D._mongodb['webmodformspace'].update({'_id':r['_id']},{"$set":{'space_orbit_label':new_label}})

def fix_orbit_labels_ap(D):
    from compmf.character_conversions import conrey_character_number_to_conrey_galois_orbit_number
    for r in D._aps.find({'conrey_galois_orbit_number':{"$exists":False}}):
        N = r.get('N'); k=r.get('k'); ci=r.get('cchi')
        on = conrey_character_number_to_conrey_galois_orbit_number(N,ci)[1]
        D._aps.update({'_id':r['_id']},{"$set":{'conrey_galois_orbit_number':int(on)}})
        label = r.get('hecke_orbit_label')
        wmf_logger.debug("{0}".format(label))


def fix_spaces(D):
    for r in D._mongodb['webmodformspace'].find():
        label = r['space_orbit_label']
        q = D._modular_symbols.find_one({'space_orbit_label':label})
        if q is None:
            wmf_logger.critical("No space {0} in D._modular_symbols!".format(label))
        else:
            if q['dimn']<>r['dimension']:
                D._mongodb['webmodformspace'].update({'_id':r['_id']},{"$set":{'dimension':q['dimn']}})

#### LAtest dimension table routines
#### These should be used.
def update_database_of_dimensions(D,nrange=[1,500],krange=[2,20]):
    r"""
    Update the dimension table in the collection 'dimension_table'
    """
    from compmf.character_conversions import dirichlet_character_conrey_galois_orbits_reps
    from sage.all import  Gamma1
    C = D._mongodb['dimension_table']
    for n in range(nrange[0],nrange[1]+1):
        orbits = dirichlet_group_conrey_galois_orbits_numbers(n)
        G = Gamma1(n)
        for xi in range(len(orbits)):
            orbit = orbits[xi]
            orbit.sort()
            orbit = map(int,orbit)
            x = orbit[0]
            xc = conrey_character_from_number(n,x).sage_character()
            if xc.is_even():
                parity = int(1)
            else:
                parity = int(-1)
            for k in range(krange[0],krange[1]+1):
                label = '{0}.{1}.{2}'.format(n,k,x)
                space_orbit_label = '{0}.{1}.{2}'.format(n,k,xi)
                if C.find({'space_orbit_label':space_orbit_label}).count()==0:
                    r = C.find_one({'space_label':label}); fid = None
                    if not r is None:
                        d_new = r['d_newf']
                        d_mod = r['d_mod']
                        d_cusp= r['d_cusp']
                        d_eisen = r['d_eis']
                        fid = r['_id']
                    elif n <= 2:
                        d_new = G.dimension_new_cusp_forms(k)
                        d_mod = G.dimension_modular_forms(k)
                        d_eisen = G.dimension_eis(k)
                        d_cusp = G.dimension_cusp_forms(k)
                    else:
                        d_new = G.dimension_new_cusp_forms(k,eps=xc)
                        d_mod = G.dimension_modular_forms(k,eps=xc)
                        d_eisen = G.dimension_eis(k,eps=xc)
                        d_cusp = G.dimension_cusp_forms(k,eps=xc)
                    cw= D._mongodb['webmodformspace'].find({'space_orbit_label':space_orbit_label}).count()
                    cm= D._modular_symbols.find({'space_orbit_label':space_orbit_label,'complete':{"$gt":int(data_record_checked_and_complete-1)}}).count()
                    r = {'space_orbit_label':space_orbit_label,
                         'space_label':label,
                         'character_orbit':orbit,
                         'character_parity':parity,
                         'level':int(n),
                         'weight':int(k),
                         'cchi':int(x),
                         'd_mod':int(d_mod),
                         'd_cusp':int(d_cusp),
                         'd_newf':int(d_new),
                         'd_eis':int(d_eisen),
                         'in_wdb':int(cw),
                         'in_msdb':int(cm)}
                    if not fid is None:
                        C.update({'_id':fid},{"$set":r})
                    else:
                        C.insert(r)
        # For Gamma1 -- total of the above
        num_orbits = len(orbits)
        for k in range(krange[0],krange[1]+1):
            label = '{0}.{1}'.format(n,k)
            r = C.find_one({'gamma1_label':label}); fid = None
            if r is None:
                d_new = G.dimension_new_cusp_forms(k)
                d_mod = G.dimension_modular_forms(k)
                d_eisen = G.dimension_eis(k)
                d_cusp = G.dimension_cusp_forms(k)
            else:
                d_mod = r['d_mod']
                d_new = r['d_newf']
                d_cusp = r['d_cusp']
                d_eisen = r['d_eis']
                fid = r['_id']
            num_in_db = len(D._mongodb['webmodformspace'].find({'level':int(n),'weight':int(k)}).distinct('character'))
            r = {'gamma1_label':label,
                 'one_in_wdb': int(num_in_db)>0,
                 'level':int(n),
                 'weight':int(k),
                 'd_mod':int(d_mod),
                 'd_cusp':int(d_cusp),
                 'd_newf':int(d_new),
                 'd_eis':int(d_eisen),
                 'all_in_db': int(num_in_db) >= (num_orbits)
                 }
            if fid is None:
                C.insert(r)
            else:
                C.update({'_id':fid},{"$set":r})

    print "Updated table!"

def update_existing_database_of_dimensions(D,nrange=[1,500],krange=[2,20],
        only_Gamma0=True,verbose=0,check_db=False):
    r"""
    Update the dimension table in the collection 'dimension_table'
    with information about what exists in our webmodform 
    """
    from compmf.character_conversions import dirichlet_character_conrey_galois_orbits_reps
    from sage.all import  Gamma1
    C = D._mongodb['dimension_table']
    #cw= D._mongodb['webmodformspace'].find({'space_orbit_label':space_orbit_label}).count()
    #cm= D._modular_symbols.find({'space_orbit_label':space_orbit_label,'complete':{"$gt":int(data_record_checked_and_complete-1)}}).count()
    for n in C.find().distinct('level'):
        if nrange != []:
            if n < nrange[0] or n > nrange[1]:
                continue
        for xi in C.find({'level':n}).distinct('cchi'):
            if verbose>0:
                wmf_logger.debug("xi={0}".format(xi))
            x = conrey_character_from_number(n,xi)
            if x.is_even():
                parity = int(1)
            else:
                parity = int(-1)
            for r in C.find({'level':int(n),'cchi':xi}).sort([('weight',int(1))]):
                if krange != []:
                    if r['weight'] < krange[0] or r['weight'] > krange[1]:
                        continue
                q = {
                    'character_parity':parity
                }
                if check_db:
                    q['in_wdb'] = D._mongodb['webmodformspace'].find({'version':float(1.3),'space_label':r['space_label']}).count()
                    q['in_msdb']= D._modular_symbols.find({'space_label':r['space_label'],'complete':{"$gt":int(data_record_checked_and_complete-1)}}).count()
                dr = C.find_one({'space_label':r['space_label']})
                if not dr is None:
                    fid = dr['_id']
                    u = C.update({'_id':fid},{"$set":q})
                    if verbose > 0:
                        wmf_logger.debug("updated {0} for label={1}".format(u,r['space_label']))
                else:
                    wmf_logger.debug("dimension record of {0} not in database!".format(r['space_label']))
                        
        # For Gamma1 -- total of the above
    wmf_logger.info("Updated the Gamma0 table!")
    if only_Gamma0:
        return 
    for n in D._mongodb['webmodformspace'].distinct('level'):
        orbits = dirichlet_group_conrey_galois_orbits_numbers(n)
        num_orbits = len(orbits)
        for k in D._mongodb['webmodformspace'].find({'level':n}).distinct('weight'):
            label = '{0}.{1}'.format(n,k)
            r = C.find_one({'gamma1_label':label}); fid = None
            if r is None:
                continue
            fid = r['_id']
            num_in_db = len(D._mongodb['webmodformspace'].find({'level':int(n),'weight':int(k)}).distinct('character'))
            r = {
                 'one_in_wdb': int(num_in_db)>0,
                 'all_in_db': int(num_in_db) >= (num_orbits)
                 }
            C.update({'_id':fid},{"$set":r})

    print "Updated table!"



def check_all_in_db(D):
    C = D._mongodb['dimension_table']
    for r in C.find({'gamma1_label':{"$exists":True}}):
        fid = r['_id']
        label = r.get('gamma1_label')
        d_newf = r['d_newf']
        indb = 0
        for q in C.find({'space_label':{"$exists":True},'level':r['level'],'weight':r['weight']}):
            indb += q['d_newf']*len(q['character_orbit'])
        if indb == d_newf:
            C.update({'_id':fid},{"$set":{'all_in_db':int(1)}})
        else:
            C.update({'_id':fid},{"$set":{'all_in_db':int(0)}})
        r1 = C.find_one({'level':r['level'],'weight':r['weight'],'space_label':{"$exists":True}})
        if not r1 is None:
            indb = r1['in_wdb']
            if indb == 1:
                C.update({'_id':fid},{"$set":{'one_in_wdb':int(1)}}) ## This simply means that at least one is in the db

def add_character(D):
    C = D._mongodb['dimension_table']
    for r in C.find({'space_label':{"$exists":True}}):
        fid = r['_id']
        label = r.get('space_label')
        N,k,i = label.split(".")
        C.update({'_id':fid},{"$set":{'cchi':int(i)}})

def fixing_spaces(D,N,k,i):
    r"""
    Check the space with this label and make sure that the coefficients are associated with the same space...
    PROBLEM: Data in files are sometimes at wrong place...
    """
    if isinstance(N,basestring):
        N,k,i=N.split(".")
    M = D.get_ambient(N,k,i)
    x1 = M.character()
    x2 = dirichlet_character_sage_from_conrey_character_number(N,i)
    if x1 != x2:
        return False
    factors = D.get_factors(N,k,i,'all')
    d = 0
    for F in factors:
        wmf_logger.debug("Checking {0}".format(F))
        if not F.is_submodule(M):
            wmf_logger.debug("Newform {0} is not a subspace of ambient space!".format(F))
            return False
        if not F.is_new() or not F.is_cuspidal():
            wmf_logger.debug("Newform {0} is not new or not cuspidal!".format(F))
            return False
        # check coefficients fields
        K1 = F.eigenvalue(1).parent()
        a = D.get_aps(N,k,i,d)
        for prec in a.keys():
            v = a[prec][1]
            K2 = v[0].parent()
            if K1 != QQ:
                if not K1.is_isomorphic(K2):
                    wmf_logger.debug("Parents are not isomorphic: K1={0} and K2={1}".format(K1,K2))
            else:
                if K2 != QQ:
                    wmf_logger.debug("Parents are not isomorphic: K1={0} and K2={1}".format(K1,K2))
                    # so the coefficients have the correct field...
            if len(v) <> F.dimension():
                wmf_logger.debug("Length of eigenvector is incorrect!")
                return False
            E = a[prec][0]
            if E.nrows()<>prime_pi(prec):
                wmf_logger.debug("Number of rows in E is incorrect!")
                return False
        d+=1
        ## We can also check the coefficients explicitly
        c1 = [F.eigenvalue(n) for n in range(1,F.sturm_bound()+1)]
        c2 = (E*v)[1:F.sturm_bound()+1]
        if c1 <> c2:
            wmf_logger.debug("Coefficients are not equal!")
    return True


def add_name_to_AL(D):
    C = D._mongodb['Atkin_Lehner.files']
    for r in C.find({'space_label':{"$exists":True}}):
        fid = r['_id']
        label = r.get('space_label')
        N,k,i = label.split(".")
        C.update({'_id':fid},{"$set":{'cchi':int(i)}})


def check_files_of_coefficients(D,s="",ncpus=32):
    args=[]
    C = D._mongodb['file_checked']
    l = []
    wmf_logger.debug("s={0}".format(s))
    if s <> "":
        s += " AND newforms > 0 AND maxp>0"
    else:
        s = "newforms > 0 AND maxp>0"
    for N,k,ci,nd,maxn in D._db.known(s):
        l.append((N,k,ci,nd,maxn))
    wmf_logger.debug("l has {0} elements".format(len(l)))
    for N,k,ci,nd,maxn in l:
        #s = {'N':N,'k':k,'ci':ci,'d':{"$in": [int(nd),int(nd-1)]},'maxn':maxn}
        #wmf_logger.debug("search for {0}".format(s))
        #    continue
        if maxn == 0 or nd==0:
            continue
        for d in range(nd):
            s = {'N':N,'k':k,'ci':ci,'d':int(d),'maxn':maxn,'checked':False}
            if C.find(s).count()==0:
                continue
            args.append((N,k,ci,d,maxn))
    wmf_logger.debug("Will check {0} records!".format(len(args)))
    if ncpus >= 32:
        l = check_coefficients_32_record(args)
    elif ncpus >= 16:
        l = check_coefficients_16_record(args)
    elif ncpus >= 8:
        l = check_coefficients_8_record(args)
    else:
        l = []
        for r in args:
            l.append(check_coefficients_one_record(**r))
    return list(l)

@parallel(ncpus=8)
def check_coefficients_8_record(N,k,ci,d,maxn,datadir='/home/stromberg/data/modforms-db/',host='localhost',port=int(37010),dryrun=False):
    return check_coefficients_one_record(N,k,ci,d,maxn,datadir=datadir,host=host,port=port,dryrun=dryrun)

@parallel(ncpus=16)
def check_coefficients_16_record(N,k,ci,d,maxn,datadir='/home/stromberg/data/modforms-db/',host='localhost',port=int(37010),dryrun=False):
    return check_coefficients_one_record(N,k,ci,d,maxn,datadir=datadir,host=host,port=port,dryrun=dryrun)

@parallel(ncpus=32)
def check_coefficients_32_record(N,k,ci,d,maxn,datadir='/home/stromberg/data/modforms-db/',host='localhost',port=int(37010),dryrun=False):
    return check_coefficients_one_record(N,k,ci,d,maxn,datadir=datadir,host=host,port=port,dryrun=dryrun)

def check_coefficients_one_record(N,k,ci,d,maxn,datadir='/home/stromberg/data/modforms-db/',host='localhost',port=int(37010),dryrun=False):
    r"""
    Check coefficients in the file for one record

    """
    D = CompMF(datadir=datadir,host=host,port=port)
    C = D._mongodb['file_checked']
    a = D.get_aps(N,k,ci,d,sources=['files'],prec_needed='all')
    try:
        pprecs = deepcopy(a.keys())
        wmf_logger.debug("Got pprecs: {0}".format(pprecs))        
    except AttributeError as e:
        wmf_logger.debug("Can not get aps: {0}".format(e))
        pprecs=[]
    if pprecs == []:
        C.update({'N':int(N),'k':int(k),'ci':int(ci),'d':int(d),'maxn':int(maxn)},{"$set":{'checked':True,'pprec':[]}})
        return
    for pprec in pprecs:
        if C.find({'N':int(N),'k':int(k),'ci':int(ci),'d':int(d),'maxn':int(maxn),'pprec':[int(pprec[0]),int(pprec[1])],'checked':True}).count()>0:
            continue
        wmf_logger.debug("Checking {0}".format((N,k,ci,d,pprec)))
        E,v = a[pprec][0:2]
        dim = dimension_new_cusp_forms(conrey_character_from_number(N,ci).sage_character(),k)
        if E.ncols() <> len(v) or dim <> len(v):
            wmf_logger.debug("E.ncols={0} len(v)={1} dim={2}".format(E.ncols(),len(v),dim))
            if pprec[0]>0:
                fname = D._db.factor_aplist(N,k,ci,d,False,pprec[0],pprec[1])
                wmf_logger.critical("Removing file for {0} : {1}".format((N,k,ci,d,pprec),fname))
                if dryrun:
                    continue
                D._db.delete_file(fname)
                C.remove({'N':int(N),'k':int(k),'ci':int(ci),'d':int(d),'maxn':int(maxn),'pprec':[int(pprec[0]),int(pprec[1])]})
                continue
            # better check MongoDB as well:
            s = {'N':int(N),'k':int(k),'cchi':int(ci),'newform':int(d)}
            for x in D._aps.find(s):
                E,v = D.load_from_mongo('ap',x['_id'])
                if len(v)<>dim or E.ncols()<> dim or len(v)<>E.ncols():
                    D.delete_from_mongo('ap',x['_id'])
                    wmf_logger.debug("Remove {0}".format(s))
        try:
            c = multiply_mat_vec(E,v)
            ok = True            
        except TypeError:
            # in this case we might have problems with MongoDB as well
            s = {'N':int(N),'k':int(k),'cchi':int(ci),'newform':int(d),'prec':int(maxn)}
            r = D._aps.find_one(s)
            if not r is None:
                D._mongodb['mongo_problem'].insert(s)
                D.delete_from_mongo('ap',r['_id'])
                wmf_logger.debug("Removed record {0} from mongo!".format(r))
            ok = False
            c = []
        if len(c) <> prime_pi(pprec[1])-prime_pi(pprec[0]):
            ok = False
        if ok:
            p = nth_prime(pprec[0]+1)
            if c[0].parent() is QQ:
                ap = abs(c[0])/RR(p)**(RR(k-1)/RR(2))
            else:
                ap = abs(c[0].complex_embedding())/RR(p)**(RR(k-1)/RR(2))
            if ap > 2.0:
                ok = False
        if not ok:
            fname = D._db.factor_aplist(N,k,ci,d,False,pprec[0],pprec[1])
            wmf_logger.critical("Removing file for {0} : {1}".format((N,k,ci,d,pprec),fname))
            if dryrun:
                continue
            D._db.delete_file(fname)
            C.remove({'N':int(N),'k':int(k),'ci':int(ci),'d':int(d),'maxn':int(maxn),'pprec':[int(pprec[0]),int(pprec[1])]})
            pprecs.remove(pprec)
            if pprec==maxn:
                wmf_logger.critical("Removing from known db also!")
                db = sqlite3.connect(D._known_db_file)
                cursor = db.cursor()
                # remove old record
                cursor.execute("DELETE FROM known WHERE N={0} AND k={1} AND i={2} AND newforms={3} and maxp={4}".format(N,k,ci,d,maxn))
                if pprecs != []:
                    prec_max = max(pprecs)
                    # insert new with possible smaller precision given
                    cursor.execute("INSERT INTO known VALUES(?,?,?,?,?)", (N,k,ci,d,prec_max))
                db.commit()
        s = {'N':int(N),'k':int(k),'ci':int(ci),'d':int(d),'maxn':int(pprec[1])}
        q = C.find(s)
        wmf_logger.debug("Updating record for s={0}. count={1}".format(s,q.count()))
        if q.count()>0:
            new_rec = {"$set":{'checked':True,'pprec':[int(pprec[0]),int(pprec[1])]}}
            for r in q:
                C.update({'_id':r['_id']},new_rec,multi=True)
        else:
            new_rec = {'N':int(N),'k':int(k),'ci':int(ci),'d':int(d),'maxn':int(pprec[1]),'pprec':[int(pprec[0]),int(pprec[1])],'checked':True} 
            fid = C.insert(new_rec)
            wmf_logger.debug("inserted new record: {0}".format(fid))            
            #D._db.delete_file(apfile)



def fix_aps_nmax(D,nmax=10,ncpus=1,verbose=0):
    args = []
    for r in D._aps.find({'is_converted':False,'N':{"$lt":int(nmax)}}).sort([('N',int(1)),('k',int(1))]):
        args.append(r['_id'])
    wmf_logger.debug("num of args={0}".format(len(args)))
    if ncpus>=32:
        return list(fix_aps_parallel_32(args))
    elif ncpus>=8:
        return list(fix_aps_parallel_8(args))
    else:
        l = []
        for fid in args:
            l.append(fix_aps_parallel_one(fid,verbose=verbose))
        return list(l)

@parallel(ncpus=8)
def fix_aps_parallel_8(fid):
    return fix_aps_parallel_one(fid)


@parallel(ncpus=32)
def fix_aps_parallel_32(fid):
    return fix_aps_parallel_one(fid)

def fix_aps_parallel_one(fid,verbose=0):
    r"""
    Now use the new routines...
    """
    from sage.all import prime_pi,nth_prime
    D = MongoMF(host='localhost',port=int(37010))
    s = {'_id':fid,'is_converted':False}
    r = D._aps.find_one(s)
    aps_new = D._mongodb['ap2'] 
    if r is None:
        return
    wmf_logger.debug("checking: {0}, id={1}".format(r['hecke_orbit_label'],r['_id']))
    if verbose > 0:
        print "record = ",r
    try:
        E,v = D.load_from_mongo('ap',fid)
        #if verbose > 0:
        wmf_logger.debug("Checking E and v for {0}".format(r['hecke_orbit_label']))
        if E.base_ring() <> v.base_ring():
            wmf_logger.debug("converting E")
            #l = Modf_changevar_Ev(E,v,NF=None,Bfacto=10^6)
            #E=l[0]; v=l[1]
            E = convert_matrix_to_extension_fld(E,v.base_ring())
            #gen = str(l[2]); emb=str(l[3]); label=str(l[4])
        a2 = (E*v)[0] # sum(E[0,x]*v[x] for x in range(len(v)))
        k = r['k']
        if a2.parent() != QQ:
            a2 = max(map(abs,a2.complex_embeddings()))
        else:
            a2 = abs(a2)
        a2 = a2/RR(2.0)**(RR(k-1)/RR(2))
        if a2 > 2.0:
            wmf_logger.critical("a(2)={0} does not satisfy the Ramanujan bound".format(a2))
            D._aps.update({'_id':r['_id']},{"$set":{'recheck':True}})
            return
        rr = deepcopy(r)
        rr.pop('_id')
        if not r.has_key('nmax'):
            nmax = int(nth_prime(E.nrows()+1)-1)
            nmin = int(0)
            rr['nmax'] = nmax
            rr['nmin'] = nmin
        if E.nrows()<>prime_pi(rr['nmax'])-prime_pi(rr['nmin']) :
            clogger.critical("Record {0} does not have correct  number of coefficients".format(r['hecke_orbit_label']))
            D._aps.update({'_id':r['_id']},{"$set":{'recheck':True}})
            return
        fs_ap = gridfs.GridFS(D._mongodb, 'ap')
        rr['is_converted']=True
        #rr['field_gen']=gen
        #rr['field_emb']=emb
        #rr['field_label']=label
        rold = D._aps.find_one({'hecke_orbit_label':r['hecke_orbit_label'],'prec':r['prec']})
        if not rold is None:
            if float(pymongo.version)==2.8:
                D._aps.remove(rold['_id'])
            else:
                D._aps.delete_one(rold['_id'])
        t = fs_ap.put(dumps( (E,v)),**rr)
        return t
    #res = D._aps.update({'_id':r['_id']},{"$set":{'nmax':nmax,'nmin':nmin,'pmax':int(nn)}})            
    except Exception as e:
        wmf_logger.critical("Wrongly formatted record!. Error{0}".format(e))
        D._aps.update({'_id':r['_id']},{"$set":{'recheck':True}})
        #wmf_logger.debug("Removing record {0} which has old class number field elements!".format(r['hecke_orbit_label']))
        #if not 'out of memory' in str(e):
        #    D.delete_from_mongo('ap',r['_id'])
        raise ValueError,"Could not load E,v for {0}. Error:{1}".format(r['hecke_orbit_label'],e)




def change_base_ring(D,nmax=None,nmin=0,nlimit=None,ncpus=1,verbose=0):
    from sage.all import Integer
    args = []
    C = D._mongodb['converted_E']
    C1 = D._mongodb['not_converted_E']
    s = {}
    if isinstance(nmax,(int,Integer)):
        s['N']={"$lt":int(nmax),"$gt":int(nmin)}
    if nlimit is None:
        q = D._aps.find(s).sort([('N',int(1)),('k',int(1))])
    else:
        q = D._aps.find(s).sort([('N',int(1)),('k',int(1))]).limit(int(nlimit))
    for r in q:
        label = r.get('hecke_orbit_label')
        if label is None:
            wmf_logger.debug("No label for r={0}".format(r))
        elif C.find({'hecke_orbit_label':label}).count()>0:
            continue
        if C1.find({'hecke_orbit_label':label}).count()==0:
            C1.insert({'hecke_orbit_label':label})
        args.append(r['_id'])
    wmf_logger.debug("Will change ring for {0} records!".format(len(args)))
    if ncpus>=32:
        return list(change_base_ring_32(args))
    elif ncpus>=16:
        return list(change_base_ring_16(args))
    elif ncpus>=8:
        return list(change_base_ring_8(args))
    else:
        l = []
        for fid in args:
            l.append(change_base_ring_one(fid))
        return list(l)

@parallel(ncpus=4)
def change_base_ring_4(fid):
    return change_base_ring_one(fid)

@parallel(ncpus=8)
def change_base_ring_8(fid):
    return change_base_ring_one(fid)

@parallel(ncpus=16)
def change_base_ring_16(fid):
    return change_base_ring_one(fid)

@parallel(ncpus=32)
def change_base_ring_32(fid):
    return change_base_ring_one(fid)


def change_base_ring_one(fid):
    D = MongoMF(host='localhost',port=int(37010))
    s = {'_id':fid}
    r = D._aps.find_one(s)
    C = D._mongodb['converted_E']
    C1 = D._mongodb['not_converted_E']
    if r is None:
        return
    if not C.find_one({'hecke_orbit_label':r['hecke_orbit_label'],'aid':r['_id']}) is None:
        return
    wmf_logger.debug("want to change base ring for {0}".format(r['hecke_orbit_label']))
    try:
        E,v = D.load_from_mongo('ap',fid)
        if E.base_ring() == v.base_ring():
            C.insert({'hecke_orbit_label':r['hecke_orbit_label'],'aid':r['_id']})
        else:
            EE = convert_matrix_to_extension_fld(E,v.base_ring())
            #D.delete_from_mongo('ap',r['_id'])
            # Insert an updated version of E with changed base ring
            fs_ap = gridfs.GridFS(D._mongodb, 'ap')
            rr = deepcopy(r)
            rr.pop('_id')
            t = fs_ap.put(dumps( (EE,v)),**rr)
            # N=r['N'],k=r['k'],chi=r['chi'],cchi=r['cchi'],
            # character_galois_orbit=r['character_galois_orbit'],
            # conrey_galois_orbit_number=r['conrey_galois_orbit_number'],
            # newform=r['newform'],
            # hecke_orbit_label=r['hecke_orbit_label'],

            # nmin=r['nmin'],nmax=r['nmax'],
            # cputime = r['cputime'],
            # sage_version = r['sage_version'],
            # ambient_id = r['ambient_id'])
            if t is not None:
                # delete old record
                D.delete_from_mongo('ap',r['_id'])
                C.insert({'hecke_orbit_label':r['hecke_orbit_label'],'aid':r['_id']})
                C1.remove({'hecke_orbit_label':r['hecke_orbit_label']})
                wmf_logger.debug("did change base ring for {0}".format(r['hecke_orbit_label']))
            else:
                wmf_logger.debug("Could not change base ring  for {0}".format(r['hecke_orbit_label']))
    except Exception as e:
        if not 'out of memory' in str(e):
            wmf_logger.debug("Removing record {0} which has old class number field elements!".format(r['hecke_orbit_label']))
            #    D.delete_from_mongo('ap',r['_id'])
        raise ValueError,"Could not load E,v for {0}. Error:{1}".format(r['hecke_orbit_label'],e)
    return r['hecke_orbit_label']

def clear_checked(D):
    for r in D._mongodb['file_checked'].find({'checked':False}):
        s = " N={0} AND k={1} AND i={2} and maxp={3}".format(r['N'],r['k'],r['ci'],r['maxn'])
        q = D._db.known(s)
        if len(list(q)) == 0:
            fid = r['_id']
            D._mongodb['file_checked'].remove({'_id':fid})
    # for N,k,ci,nd,maxn in D._db.known(s):
    #     l.append((N,k,ci,nd,maxn))
    # for r in D._mongodb['file_checked'].find({'checked':True}):
    #     N = int(r['N']); k= int(r['k']); ci=int(r['ci']); d=int(r['d'])
    #     maxn=int(r['maxn']); pprec=r['pprec']
    #     s = deepcopy(r)
    #     s.pop('_id');
    #     s['checked']=False
    #     D._mongodb['file_checked'].remove(s)


def check_aps_in_mongo(D,nmin=1,nmax=10,nlim=10):
    i = 0
    C=D._mongodb['aps_mongo_checked']
    args = []
    if nlim > 0:
        for q in D._aps.find({'N':{"$lt":int(nmax)+1,"$gt":int(nmin-1)}}).sort([('N',int(1)),('k',int(1))]).limit(nlim):
            N=q['N']; k=q['k']; ci=q['cchi']; fid=q['_id']
            if C.find({'record_id':fid}).count()>0:
                continue
            args.append(fid)
    else:
        for q in D._aps.find({'N':{"$lt":int(nmax)+1,"$gt":int(nmin-1)}}).sort([('N',int(1)),('k',int(1))]):
            N=q['N']; k=q['k']; ci=q['cchi']; fid=q['_id']
            if C.find({'record_id':fid}).count()>0:
                continue
            args.append(fid)
    return list(check_aps_in_mongo32(args))

@parallel(ncpus=16)
def check_aps_in_mongo32(fid):
    import sage
    i = 0
    D = MongoMF(host='localhost',port=int(37010))
    C=D._mongodb['aps_mongo_checked']
    if C.find({'record_id':fid}).count()>0:
        return
    sage.modular.modsym.modsym.ModularSymbols_clear_cache()
    for q in D._aps.find({'_id':fid}):
        N=q['N']; k=q['k']; ci=q['cchi']; fid=q['_id']
        ambient_id = q['ambient_id']; label = q['hecke_orbit_label']
        prec = q.get('prec',int(0))
        wmf_logger.debug("Checking:{0}".format(label))
        try: 
            M = D.load_from_mongo('Modular_symbols',ambient_id)
            if M is None:
                M = D.get_ambient(N,k,ci,sources=['mongo'])            
        except Exception as e:
            wmf_logger.critical("Error with import for {0}. Probably incompatible Sage versions! ERROR:{1}".format(label,e))
            continue
        ok = False
        if not M is None:
            E,v=D.load_from_mongo('ap',fid)
            c = E*v
            K=v.base_ring()
            DEC = M.new_subspace().cuspidal_subspace().decomposition()
            d = q.get('newform',-1)
            S = None
            if d <0 or len(DEC)<=d:
                wmf_logger.critical("Incorrect newform number for {0}. d={1} and len(decomp)={2}\n DEC={3}".format(label,d,len(DEC),DEC))
            else:
                try:
                    S = DEC[d]
                    E1,v1=S.compact_system_of_eigenvalues([2])
                    K1 = v1.base_ring()
                except IndexError as e:
                    wmf_logger.critical("Incorrect decomposition for {0}! d={1} and len(decomp)={2}\n DEC={3}".format(label,d,len(DEC),DEC))
            if S is None:
                ok = False
            elif len(v)<>len(v1):
                wmf_logger.critical("length are different! len(v(mongo))={0}, len(v(dec))={1}. Need to remove!".format(len(v),len(v1)))
                ok = False
            elif not (K1==QQ and K == QQ):
                if K1 == QQ and K != QQ:
                    wmf_logger.critical("K1=Q, K<>Q")
                    ok = False
                elif K1 != QQ and K == QQ:
                    wmf_logger.critical("K1<>Q, K==Q")                    
                    ok = False
                elif not K1.is_isomorphic(K):
                    wmf_logger.critical("K1 !~ K2")
                    ok = False
                else:
                    ok = True
                if not ok:
                    wmf_logger.critical("parents of v are different! need to remove!")
            else:
                if map(lambda x:x.norm(),v)==map(lambda x:x.norm(),v1):
                    ok = True
        else:
            wmf_logger.critical("No ambient space for coefficients with {0}".format(q['hecke_orbit_label']))
        if not ok:
            #D.delete_from_mongo('ap',fid)
            C.update({'record_id':fid},{'hecke_orbit_label':label,'prec':prec,'record_id':fid,'ok':False},upsert=True)            
            wmf_logger.critical("Removing record for {0}".format(q['hecke_orbit_label']))
        else:
            prec = q.get('prec',int(0))
            if prec == 0:
                wmf_logger.critical("Record without prec: {0}".format(label))
            C.update({'record_id':fid},{'hecke_orbit_label':label,'prec':prec,'record_id':fid,'ok':True},upsert=True)


def check_ambient_in_mongo(D,nmin=1,nmax=10,nlim=10):
    import sage
    from sage.all import ModularSymbols
    args = []
    C=D._mongodb['ambient_mongo_checked']
    for q in D._modular_symbols.find({'N':{"$lt":int(nmax)+1,"$gt":int(nmin-1)}}).sort([('N',int(1)),('cchi',int(1)),('k',int(1))]):
        fid=q['_id']
        if C.find({'record_id':fid,'ok':True}).count()>0:
            continue
        args.append(fid)
    return list(check_ambient_in_mongo16(args))

@parallel(ncpus=16)
def check_ambient_in_mongo16(fid):
    import sage
    from sage.all import ModularSymbols
    D = MongoMF(host='localhost',port=int(37010))
    C=D._mongodb['ambient_mongo_checked']
    if C.find({'record_id':fid,'ok':True}).count()>0:
        return 
    for q in D._modular_symbols.find({'_id':fid}):
        N=q['N']; k=q['k']; ci=q['cchi']; fid=q['_id']
        label = q['space_label']
        x = conrey_character_from_number(N,ci)
        sage.modular.modsym.modsym.ModularSymbols_clear_cache()
        #M = ModularSymbols(x.sage_character(),k,sign=1)
        try:
            M1 = D.load_from_mongo('Modular_symbols',fid)
        except Exception as e:
            wmf_logger.critical("Error with import for {0}. Probably incompatible Sage versions! ERROR:{1}".format(label,e))
            continue
        ci_true = sage_character_to_conrey_character(M1.character()).number()
        if ci <> ci_true:
            wmf_logger.critical("x1<>x!: {0} \n x={1}\n x1={2}".format(label,x.sage_character(),M1.character()))
            # change the space label to the correct one...
            true_label = '{0}.{1}.{2}'.format(N,k,ci_true)
            wmf_logger.debug("changing labels from {0} to {1}".format(label,true_label))
            D._modular_symbols.update({'_id':fid},{"$set":{'space_label':true_label,'cchi':ci_true}})
            ## files are ok
            for r in D._newform_factors.find({'N':N,'k':k,'cchi':ci}):
                F = D.load_from_mongo('Newform_factors',r['_id'])
                cii_true = int(sage_character_to_conrey_character(F.character()).number())
                if cii_true <> ci:
                    nn,kk,lab = r['hecke_orbit_label'].split('.')
                    lab = lab.replace(str(ci),str(cii_true))
                    wmf_logger.critical("Updating character nr. in newform from {0} to {1}".format(r['hecke_orbit_label'],lab))
                    
                    D._newform_factors.update({'_id':r['_id']},{"$set":{'cchi':cii_true,'hecke_orbit_label':'{0}.{1}.{2}'.format(N,k,lab)}})
                    for ap_r in D._aps.find({'N':N,'k':k,'cchi':ci,'newform':int(r['newform'])}):
                        wmf_logger.critical("ALso updating ap's")
                        D._newform_factors.update({'_id':ap_r['_id']},{"$set":{'cchi':cii_true,'hecke_orbit_label':'{0}.{1}.{2}'.format(N,k,lab)}})
                
            #C.update({'record_id':fid},{'space_label':label,'record_id':fid,'ok':False},upsert=True)
            # check if this space is not the representative
            #
            #if ci <> min(q['character_galois_orbit']):
            #    wmf_logger.debug("Character not orbit representative!")
            #s = {'cchi':min(q['character_galois_orbit']),'N':N,'k':k}
            #r = D._modular_symbols.find_one(s)
            #if not r is None:
            #    wmf_logger.debug("Space with orbit representative exists Label={0}! Deleting this record! {1} orbit={2}".format(r['space_label'],q['space_label'],q['character_galois_orbit']))
            #    #D.delete_from_mongo('Modular_symbols',fid)
            #    #C.remove({'record_id':fid})
        else:
            t = C.update({'record_id':fid},{'space_label':label,'record_id':fid,'ok':True},upsert=True)        
  
# @parallel(ncpus=16)
# def check_ambient_in_mongo16(fid):
#     import sage
#     from sage.all import ModularSymbols
#     D = MongoMF(host='localhost',port=int(37010))
#     C=D._mongodb['ambient_mongo_checked']
#     if C.find({'record_id':fid}).count()>0:
#         return 
#     for q in D._modular_symbols.find({'_id':fid}):
#         N=q['N']; k=q['k']; ci=q['cchi']; fid=q['_id']
#         label = q['space_label']
#         x = conrey_character_from_number(N,ci)
#         sage.modular.modsym.modsym.ModularSymbols_clear_cache()
#         #M = ModularSymbols(x.sage_character(),k,sign=1)
#         M1 = D.load_from_mongo('Modular_symbols',fid)
#         if M1.character() <> x.sage_character():
#             #wmf_logger.critical("M1<>M!: {0} \nM={1}\nM1={2}".format(label,M,M1))
#             wmf_logger.critical("x1<>x!: {0} \n x={1}\nM1={2}".format(label,M,M1))            
#             C.update({'record_id':fid},{'space_label':label,'record_id':fid,'ok':False},upsert=True)
#             # check if this space is not the representative
#             if ci <> min(q['character_galois_orbit']):
#                 wmf_logger.debug("Character not orbit representative!")
#             s = {'cchi':min(q['character_galois_orbit']),'N':N,'k':k}
#             if D._modular_symbols.find(s).count()>0:
#                 wmf_logger.debug("Space with orbit representative exists! Deleting this wrong instance!")
#                 #D.delete_from_mongo('Modular_symbols',fid)
#                 M2 = D._db.load_ambient_space(N,k,ci)
                
#                 if M2 <> M:
#                      wmf_logger.critical("M2<>M!: {0} \nM={1}\nM2={2}".format(label,M,M2))
#                 else:
#                     pass
#             t = C.update({'record_id':fid},{'space_label':label,'record_id':fid,'ok':True},upsert=True)        
  
#def check_duplicates_in_orbits(fid):
    
    


def check_twist_info(D,nmax=10,nmin=1):
    res = []
    from lmfdb.modular_forms.elliptic_modular_forms import emf_logger
    emf_logger.setLevel(100)
    for r in D._mongodb['webnewforms.files'].find():
        level = int(r['hecke_orbit_label'].split(".")[0])
        if level < nmin or level > nmax:
            continue
        fid = r['_id']
        f = D.load_from_mongo('webnewforms',fid)
        if f is None:
            res.append(r['hecke_orbit_label'])
            continue
        t = f.get('twist_info',None)
        if t is not None:
            try:
                if  len(t)==2:
                    if not isinstance(t[1][0],basestring):
                        print type(t[1][0])
                        res.append(r['hecke_orbit_label'])
            except:
                res.append(r['hecke_orbit_label'])
    wmf_logger.debug("Redoing {0} newforms!".format(len(res)))
    return list(recompute_one(res))


@parallel(ncpus=32)
def recompute_one(label):
    F=WebNewForm_computing(label,recompute=False)
    return True

def get_duplicate_keys(D):
    #C=D._mongodb['webmodformspace.files']
    C=D._newform_factors
    for r in C.find({'character_galois_orbit':{"$exists":True}}).sort([('uploadDate',int(-1))]): 
        fid = r['_id']
        label = r['hecke_orbit_label']
        q = C.find({'hecke_orbit_label':label}) #,'conrey_galois_orbit':{"$exists":False}})
        n = q.count()
        if n<=1:
            continue
        else:
            F = D.load_from_mongo('Newform_factors',fid)
            x = sage_character_to_conrey_character(F.character())
            if x.number() <> r['cchi'] or x.number() not in r['character_galois_orbit']:
                raise ValueError,"Check the record:{0}".format(label)
            wmf_logger.debug("Duplicates for {0} : {1}".format(r['hecke_orbit_label'],n))
            for x in q.sort([('uploadDate',pymongo.ASCENDING)]):
                if x['_id']==r['_id']:
                    continue
                F2 = D.load_from_mongo('Newform_factors',x['_id'])
                if F2 == F:
                    wmf_logger.debug("Removing {0} : {1}".format(x['hecke_orbit_label'],x['_id']))
                    C.delete_one({'_id':x['_id']})

def remove_faulty_records(D):
    for r in D._mongodb['webnewforms'].find({'version':float(1.3)}).sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)]):
        label = r['hecke_orbit_label']
        try:
            F = WebNewForm_computing(label,recompute=False,update_from_db=True)
        except Exception as e:
            wmf_logger.critical("ERROR:{0} {1}".format(label,str(e)))
            D._mongodb['webnewforms'].update({'_id':r['_id']},{"$set":{'fix':int(2)}})
            continue
        l = F.twist_info
        try:
            for x in F.twist_info[1]:
                if 'computing' in type(x):
                    wmf_logger.debug("Need to fix: {0}".format(label))
                    D._mongodb['webnewforms'].update({'_id':r['_id']},{"$set":{'fix':int(1)}})
                #F=WebNewForm_computing(label,recompute=True,update_from_db=True)
        except Exception as e:
            wmf_logger.debug("No twist info. Need to fix: {0}".format(label))
            D._mongodb['webnewforms'].update({'_id':r['_id']},{"$set":{'fix':int(1)}})


def remove_bad_factors(D,nmax=10,nmin=1):
    args = []
    for r in D._newform_factors.find({'checked':{"$ne":int(1)},'N':{"$gt":int(nmin-1),"$lt":int(nmax+1)}}).sort([('N',int(1)),('k',int(1))]):
        args.append(r['_id'])
    wmf_logger.debug("checking {0} records!".format(len(args)))
    res = list(remove_bad_factors_par(args))
    return True

@parallel(16)
def remove_bad_factors_par(fid):
  D = CheckingDB('/home/stromberg/data/modforms-db/',host='lmfdb')
  for r in D._newform_factors.find({'_id':fid}):
        label = r['hecke_orbit_label']
        remove = False
        if  r['cchi'] not in r['character_galois_orbit']:
            wmf_logger.critical("Problem with character and galois orbit! label={0}".format(label))
            continue
        F = D.load_from_mongo('Newform_factors',r['_id'])
        if F is None:
            wmf_logger.critical("Problem with factor is None! label={0} id={1}".format(label,r['_id']))
        if not F.is_cuspidal():
            wmf_logger.critical("Problem with factor is not cuspidal! label={0}".format(label))
            remove = True
        if not remove:
            M = D.get_ambient(r['N'],r['k'],r['cchi'])
            x = M.character()
            xc = sage_character_to_conrey_character(x)
            if xc.number()<>r['cchi']:
                wmf_logger.critical("Problem with ambient! character is not correct label={0}".format(label))
                remove=True
            if not F.is_submodule(M):
                wmf_logger.critical("Problem with factor is not submodule of ambient! label={0}".format(label))
                remove=True
        if remove:
            F1 = D._db.load_factor(r['N'],r['k'],r['cchi'],r['newform'])
            if F1 == F:
                wmf_logger.critical("Need to delete files at {0}".format(D._db.factor(r['N'],r['k'],r['cchi'],r['newform'])))
            wmf_logger.debug("Removing MongoDB record!")
            #D._newform_factors.remove({'_id':fid})
        else:
            D._newform_factors.update({'_id':fid},{"$set":{'checked':int(1)}})

#def check_directories(D):
    


def fix_cm(D,nmax=10,nmin=1,date=""):
    args = []
    #for r in D._newform_factors.find({'N':{"$gt":int(nmin-1),"$lt":int(nmax+1)}}).sort([('N',int(1)),('k',int(1))]):
    import dateutil.parser
    if date != "":
        d = dateutil.parser.parse(date)
        s = {"modification_date":{"$exists":False},'version':{"$gt":float(1.2)}}
        for r in D._mongodb['webnewforms'].find(s):
            args.append(r['hecke_orbit_label'])
        s = {"modification_date":{"$lt":d},'version':{"$gt":float(1.2)}}
        for r in D._mongodb['webnewforms'].find(s):
            args.append(r['hecke_orbit_label'])
    else:
        for r in D._mongodb['webnewforms'].find({'version':{"$gt":float(1.2)}}):
            args.append(r['hecke_orbit_label'])

    wmf_logger.debug("checking {0} records!".format(len(args)))
    res = list(fix_cm_par(args))
    return res

@parallel(32)
def fix_cm_par(label,recompute=True):
    F = WebNewForm_computing(label,recompute=recompute)
    F.set_is_cm()
    F.save_to_db()
    return F.is_cm

def add_oldspace_decompositions(D):
    args = []
    for x in D._mongodb['webmodformspace'].find({'_has_oldspace':int(0)}):
        label = x['space_label']
        #M.set_oldspace_decomposition()
        #M.save_to_db()
        #   return True
        args.append((label,D))
    return list(add_oldspace_par(args))
@parallel(32)
def add_oldspace_par(label,db):
    try:
        M = WebModFormSpace_computing(label)
    except RunTimeError as e:
        D._mongodb['webmodformspace_errors'].insert({'label':label,'error':str(e)})
    M.set_oldspace_decomposition()
    M._collection.update({'space_label':label},{"$set":{'_has_oldspace':int(1)}})
    M.save_to_db()
    return True

def fix_dimension_data(D):
    for d in D._mongodb['dimension_table'].find({"space_label":{"$exists":True}}):
        s = D._mongodb['webmodformspace'].find_one({"space_label":d['space_label']})
        if s is None:
            continue
        r = {}
        d1 = d['d_cusp']; d2 = s['dimension_cusp_forms']
        if  d1 != d2:
            print "{0} dimensions cusp are not matching {1} != {2}".format(d['space_label'],d1,d2)
            r['dimension_cusp_forms'] = d1
        d1 = d['d_mod']; d2 = s['dimension_modular_forms']
        if  d1 != d2:
            print "{0} dimensions mod are not matching {1} != {2}".format(d['space_label'],d1,d2)
            r['dimension_modular_forms'] = d1
        d1 = d['d_newf']; d2 = s['dimension_new_cusp_forms']
        if  d1 != d2:
            print "{0} dimensions newf are not matching {1} != {2}".format(d['space_label'],d1,d2)
            r['dimension_new_cusp_forms'] = d1
        if r != {}:
            D._mongodb['webmodformspace'].update({'_id':s['_id']},r)