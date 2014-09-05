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
from wmf import wmf_logger,WebNewForm_computing,WebModFormSpace_computing
from compmf import MongoMF

def generate_web_modform_spaces(level_range=[],weight_range=[],ncpus=1,host='localhost',port=int(37010)):
    r"""
    Compute and insert objects of type WebModFormSpace with levels in the given range.

    NOTE: We only compute forms which have an entry in the mongodb so you need to use the CompMF class first to generate these.

    """
    try: 
        D  = MongoMF(host,port,db)
    except pymongo.errors.ConnectionFailure as e:
        raise ConnectionFailure,"Can not connect to the database and fetch aps and spaces etc. Error: {0}".format(e.message)
    args = []; s={}
    if level_range <> []:
        if len(level_range)==1:
            s['N']=level_range[0]
        else:
            s['N']={"$gt":int(level_range[0]-1),"$lt":int(level_range[-1])}
    if weight_range <> []:
        if len(weight_range)==1:
            s['k']=weight_range[0]
        else:
            s['k']={"$gt":int(weight_range[0]-1),"$lt":int(weight_range[-1])}
    if chi_range <>[]:
        if len(chi_range)==1:
            s['chi']=chi_range[0]
        else:
            s['chi']={"$gt":int(chi_range[0]-1),"$lt":int(chi_range[-1])}
    s['complete']=int(2)
    q = D._modular_symbols.find(s).sort([('N',pymongo.ASCENDING),('k',pymongo.ASCENDING)])
    for r in q:
        N = r['N']; k=r['k']; chi=r['cchi']
        args.append(N,k,chi)
    if ncpus>=32:
        generate_one_webmodform_space32(args)
    elif ncpus>=16:
        generate_one_webmodform_space16(args)        
    elif ncpus>=8:
        generate_one_webmodform_space8(args)
    elif ncpus>=4:
        generate_one_webmodform_space4(args)
    else:
        generate_one_webmodform_space1(args)
        
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
    

def generate_one_webmodform_space1(level,weight,chi):    
    r"""
    Generates one modform space.

    """
    M = WebModFormSpace_computing(level,weight,chi)
