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
from compmf import MongoMF,MongoMF,data_record_checked_and_complete,CompMF
from compmf.utils import multiply_mat_vec,convert_matrix_to_extension_fld
from sage.misc.cachefunc import cached_function

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
            if D._mongodb['webmodformspace'].find({'level':int(N),'weight':int(k),'character':int(cchi)}).count()>0:
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
                ('hecke_orbit_label',pymongo.ASCENDING)],
         'unique':True},
        {'name': 'webnewforms.files', 'index': [
                ('hecke_orbit_label',pymongo.ASCENDING)],
         'unique':True},
        {'name': 'webmodformspace' ,'index':[
                ('level',pymongo.ASCENDING),
                ('weight',pymongo.ASCENDING),
                ('chi',pymongo.ASCENDING)],
         'unique':False},
        {'name': 'webmodformspace', 'index':[
                ('space_label',pymongo.ASCENDING)],
         'unique':True},
        {'name':  'webmodformspace.files', 'index':[
                ('space_label',pymongo.ASCENDING)],
         'unique':True},
        {'name':  'webchar', 'index': [
                ('modulus',pymongo.ASCENDING),
                ('number',ASCENDING)],
         'unique':True},
        #       {'name':  'webchar', 'index': [
        #            ('label',ASCENDING)],
        #        'unique':True},
        {'name' : 'webeigenvalues', 'index':[
                ('hecke_orbit_label',ASCENDING)],
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
from compmf.character_conversions import  dirichlet_group_conrey_galois_orbits,conrey_character_from_number,dirichlet_group_conrey_galois_orbits_numbers,dirichlet_character_sage_from_conrey_character_number
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
def update_database_of_dimensions(D,nrange=[1,500],krange=[1,20]):
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
            num_in_db = len(D._mongodb['webmodformspace'].find({'level':int(n),'k':int(k)}).distinct('character'))

            r = {'gamma1_label':label,
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



def fix_pprec_to_nmax(D,nmax=10,ncpus=1,verbose=0):
    args = []
    for r in D._aps.find({'pmax':{"$exists":False},'N':{"$lt":int(nmax)}}).sort([('N',int(1)),('k',int(1))]):
        args.append(r['_id'])
    if ncpus>=32:
        return list(fix_pprec_parallel_32(args))
    elif ncpus>=8:
        return list(fix_pprec_parallel_8(args))
    else:
        l = []
        for fid in args:
            l.append(fix_pprec_parallel_one(fid),verbose=verbose)
        return list(l)

@parallel(ncpus=8)
def fix_pprec_parallel_8(fid):
    return fix_pprec_parallel_one(fid)


@parallel(ncpus=32)
def fix_pprec_parallel_32(fid):
    return fix_pprec_parallel_one(fid)

def fix_pprec_parallel_one(fid,verbose=0):
    from sage.all import prime_pi,nth_prime
    D = MongoMF(host='localhost',port=int(37010))
    s = {'_id':fid,'pmax':{"$exists":False}}
    r = D._aps.find_one(s)
    if r is None:
        return
    if verbose > 0:
        print "record = ",r
    try:
        E,v = D.load_from_mongo('ap',fid)
        if verbose > 0:
            wmf_logger.debug("Multiplying E and v")
        if E.base_ring() <> v.base_ring():
            EE = convert_matrix_to_extension_fld(E,v.base_ring())
            D.delete_from_mongo('ap',r['_id'])
            # Insert an updated version of E with changed base ring
            fs_ap = gridfs.GridFS(D._mongodb, 'ap')
            rr = deepcopy(r)
            rr.pop('_id')
            t = fs_ap.put(dumps( (EE,v)),**rr)
            # filename=r['filename'],
            #               N=r['N'],k=r['k'],chi=r['chi'],cchi=r['cchi'],
            #               character_galois_orbit=r['character_galois_orbit'],
            #               conrey_galois_orbit_number=r['conrey_galois_orbit_number'],
            #               newform=r['newform'],
            #               hecke_orbit_label=r['hecke_orbit_label'],
            #               nmin=r['nmin'],nmax=r['nmax'],
            #               cputime = r['cputime'],
            #               sage_version = r['sage_version'],
            #               ambient_id = r['ambient_id'])
            if t is not None:
                # delete old record
                D.delete_from_mongo('ap',r['_id'])
            c = EE*v
        else:
            c = E*v
    except Exception as e:
        wmf_logger.debug("Removing record {0} which has old class number field elements!".format(r['hecke_orbit_label']))
        if not 'out of memory' in str(e):
            D.delete_from_mongo('ap',r['_id'])
        raise ValueError,"Could not load E,v for {0}. Error:{1}".format(r['hecke_orbit_label'],e)

    # first check that it satisfies Ramanujan...
    pprec = r['prec']
    if c[0].parent() is QQ:
        a2 = abs(c[0])/2.0**(RR(r['k']-1)/RR(2))
    else:
        a2 = abs(c[0].complex_embedding())/2.0**(RR(r['k']-1)/RR(2))
    if abs(a2) > 2:
        wmf_logger.debug("Removing record {0} that does not satisfy Ramanujan! a2={1}".format(r['hecke_orbit_label'],a2))
        return D.delete_from_mongo('ap',r['_id'])
    n = len(c)
    nmax = int(nth_prime(n+1)-1)
    nmin = int(0)
    nn =nth_prime(n)
    wmf_logger.debug("Updating r={0} with id={4} from pprec:{1} to nmin:{2} and nmax={3}".format(r['hecke_orbit_label'],pprec,nmin,nmax,r['_id']))
    res = D._aps.update({'_id':r['_id']},{"$set":{'nmax':nmax,'nmin':nmin,'pmax':int(nn)}})
    #,"$unset":{'prec':''}})
    #wmf_logger.debug("updated: {0}".format(res))


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
    for q in D._aps.find({'N':{"$lt":int(nmax)+1,"$gt":int(nmin)}}).sort([('N',int(1)),('k',int(1))]):
        N=q['N']; k=q['k']; ci=q['cchi']; fid=q['_id']
        ambient_id = q['ambient_id']
        wmf_logger.debug("Checking:{0}".format(q['hecke_orbit_label']))
        M = D.load_from_mongo('Modular_symbols',ambient_id)
        if M is None:
            M = D.get_ambient(N,k,ci,sources=['mongo'])
        ok = False
        if not M is None:
            E,v=D.load_from_mongo('ap',fid)
            c = E*v
            K=v.base_ring()
            S = M.new_subspace().cuspidal_subspace().decomposition()[q['newform']]
            E1,v1=S.compact_system_of_eigenvalues([2])
            K1 = v1.base_ring()
            if len(v)<>len(v1):
                wmf_logger.critical("length are different! Need to remove!")
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
            wmf_logger.critical("Removing record for {0}".format(q['hecke_orbit_label']))
            i+=1
            if i > nlim and nlim > 0:
                return 
