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
import bson
from sage.all import parallel,dumps,Gamma1
from wmf import wmf_logger,WebNewForm_computing,WebModFormSpace_computing
from compmf import MongoMF,MongoMF,data_record_checked_and_complete
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
            for x in character_conversions.dirichlet_character_conrey_galois_orbits_reps(n):
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
    q = D._mongodb[webmodformspace].find().sort([('level',pymongo.ASCENDING),('weight',pymongo.ASCENDING)])
    for r in q:
        n = r['level']; k = r['weight']; i = r['character_orbit_rep']
        n = str(n); k=str(k); i=str(i)
        if i == 1:
            if not tbl0.has_key(n):
                tbl0[n] = {}
            if tbl0[n].has_key(k):
                d,t = tbl0[n][k]
            else:
                d = r['dimension_new_cusp_forms']
            tbl0[n][k] = (int(d),int(1))
        if not tbl1.has_key(n):
            tbl1[n] = {}
        if not tbl1[n].has_key(k):
            tbl1[n][k] = {}
        if tbl1[n][k].has_key(i):
            d,t = tbl1[n][k][i]
        else:
            d = r['dimension_new_cusp_forms']
        tbl1[n][k][i] = (int(d),int(1))        
        if not tbl1[n][k].has_key("-1"):
            tbl1[n][k]["-1"] = (int(dimension_new_cusp_forms(Gamma1(int(n)),int(k))),int(1))
    for n in tbl1.keys():
        for k in tbl1[n].keys():
            dtot = 0
            for i in tbl1[n][k].keys():
                if int(i) > 0:
                    dtot+=int(tbl1[n][k][i][0])
            d = tbl1[n][k].get("-1",(0,0))[0]
            if str(d)==str(dtot):
                tbl1[n][k]["-1"]=(int(dtot),int(1))
            else:
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
    for r in D._mongodb['webnewforms'].find({'version':{"$lt":float(1.3)}}).limit(llim):
        level = r['level']; weight=r['weight']; character = r['character']; label=r['label']
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
    coll = D._mongodb['webnewforms']
    file_collection = D._mongodb['webnewforms.files']
    args = []
    for r in coll.find({"zeta_orders":{"$exists":False}}):
        args.append((r['level'],r['weight'],r['character']))
        l =  generate_one_webmodform_space32(args,recompute=True)
        return len(list(l))
    #M=webmodformspace_computing(r['level'],r['weight'],r['character'],compute=False)
    #    M.get_zetas()
    #    M.save_to_db()
