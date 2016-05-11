r"""

Move all tests/checks here... 
"""



import pymongo
import gridfs
import bson
from sage.all import parallel,dumps,Gamma1,QQ,prime_pi,RR,deepcopy,nth_prime,RR, prime_range, dimension_new_cusp_forms
from wmf import wmf_logger,WebNewForm_computing,WebModFormSpace_computing
from compmf import MongoMF,MongoMF,data_record_checked_and_complete,CompMF,CheckingDB
from compmf.utils import multiply_mat_vec,convert_matrix_to_extension_fld
from sage.misc.cachefunc import cached_function
from lmfdb.modular_forms.elliptic_modular_forms import emf_version,WebNewForm
from multiprocessing import Pool
from sage.all import ModularSymbols, ceil, RealField, previous_prime
from utils import orbit_label,orbit_index_from_label

def reset_expansions(l):
    for label in l:
        F = WebNewForm_computing(label)
        F._coefficients = {}
        F._embeddings = {}
        F.set_q_expansion()
        F.set_q_expansion_embeddings()
        try:
            F.save_to_db()
            with open("fixed.txt", "a") as fp:
                fp.write("'"+label+"',")
        except pymongo.errors.WTimeoutError:
            wmf_logger.critical("Timed out! label={0}".format(label))
            with open("failed.txt", "a") as fp:
                fp.write("'"+label+"',")
            pass
        
def test_forms_in_list(l):
    r"""
    Compare with sage coefficients (therefore possibly slow...)
    """
    for label in l:
        F = WebNewForm(label)
        M = ModularSymbols(F.character.sage_character,F.weight,sign=1)
        S = M.cuspidal_submodule().new_submodule()
        d =  orbit_index_from_label(F.label)
        nmax = max(F.level+1,20)
        f = S[d].q_eigenform(nmax+1,names='a')
        try:
            test = sum([ abs(F.coefficient(n) - f.padded_list()[n]) for n in range(nmax)])
        except IndexError as e:
            print "nmax=",nmax
            print "len(f.coefficients())=",len(f.coefficients())
            raise e
        except TypeError as e:
            print e
            test = 0
            for n in range(nmax):
                c1 = F.coefficient(n)
                c2 = f.padded_list()[n]
                if hasattr(c1,'complex_embeddings'):
                    emb1 = c1.complex_embeddings()
                    emb2 = c2.complex_embeddings()
                    for j in range(len(emb1)):
                        test += abs(emb1[j]-emb2[j])
                else:
                    test += abs(RR(c2)-RR(c1))
        if test != 0: #.count(0) != len(test):
            return False,label,F,f,test
    return True


bober_errors = ['11.11.10.b', '11.11.2.a', '11.12.1.a', '11.5.10.a', '11.7.10.b', '11.9.10.b', '12.11.5.b', '13.11.2.a', '13.12.3.a', '13.9.2.a', '14.10.9.a', '14.10.9.b', '14.11.13.a', '14.12.9.a', '14.12.9.b', '14.5.13.a', '15.10.4.a', '15.11.11.a', '15.11.14.a', '15.11.14.b', '15.12.2.a', '15.12.4.a', '15.5.11.a', '15.6.4.a', '15.7.11.a', '15.8.4.a', '15.9.11.a', '15.9.14.a', '15.9.14.b', '16.11.15.a', '16.11.15.b', '16.12.5.a', '16.7.15.a', '17.7.3.a', '17.9.3.a', '19.10.4.a', '19.11.18.a', '19.11.18.b', '19.11.2.a', '19.12.4.a', '19.12.7.a', '19.3.18.b', '19.4.7.a', '19.5.18.a', '19.5.18.b', '19.5.2.a', '19.6.4.a', '19.6.7.a', '19.7.18.a', '19.7.18.b', '19.7.2.a', '19.8.4.a', '19.9.18.a', '19.9.18.b', '19.9.2.a', '20.10.9.a', '20.11.11.a', '20.11.19.d', '20.4.9.a', '20.5.11.a', '20.6.9.a', '20.7.11.a', '20.7.19.d', '20.8.9.a', '20.9.11.a', '21.10.20.a', '21.10.20.b', '21.11.10.a', '21.11.10.b', '21.11.13.a', '21.11.2.a', '21.11.2.b', '21.11.8.a', '21.12.20.a', '21.5.10.a', '21.5.10.c', '21.5.13.a', '21.5.2.a', '21.5.2.b', '21.6.20.a', '21.7.10.a', '21.7.10.b', '21.7.13.a', '21.7.2.a', '21.7.8.a', '21.8.20.a', '21.9.10.a', '21.9.10.b', '21.9.13.a', '21.9.2.a', '21.9.2.b', '21.9.8.a', '23.11.22.b', '23.4.2.a', '23.5.5.a', '23.6.2.a', '23.7.22.a', '23.7.22.c', '23.9.22.b', '24.3.5.b', '25.10.24.c', '25.10.4.a', '25.12.24.b', '25.3.7.a', '25.4.4.a', '25.5.2.a', '25.6.4.a', '26.10.17.a', '26.10.3.a', '26.12.17.a', '26.12.3.a', '26.2.3.a', '26.4.3.a', '26.4.3.b', '26.6.17.a', '26.6.3.a', '26.8.17.a', '27.10.10.a', '27.11.26.a', '27.11.26.c', '27.11.26.d', '27.11.8.a', '27.12.10.a', '27.2.4.a', '27.3.26.a', '27.4.4.a', '27.5.2.a', '27.5.26.a', '27.6.4.a', '27.7.26.b', '27.7.26.c', '27.7.8.a', '27.8.10.a', '27.9.8.a', '28.10.27.b', '28.12.27.a', '28.12.9.a', '28.8.9.a', '29.11.12.a', '29.3.2.a', '29.4.1.b', '29.6.4.a', '29.9.12.a', '30.10.19.b', '30.12.19.b', '30.6.19.b', '30.8.19.b', '31.11.30.a', '31.11.30.b', '31.3.30.b', '31.4.1.b', '31.5.3.a', '31.5.30.b', '31.7.30.a', '31.7.30.b', '31.7.30.c', '31.9.30.a', '31.9.30.b', '32.10.17.a', '32.10.5.a', '32.11.31.a', '32.11.31.b', '32.12.17.a', '32.12.5.a', '32.5.31.b', '32.6.17.a', '32.6.5.a', '32.7.31.a', '32.7.31.b', '32.8.17.a', '32.8.5.a', '32.9.31.a', '32.9.31.b', '33.10.4.a', '33.10.4.b', '33.11.10.a', '33.11.23.a', '33.12.4.a', '33.12.4.b', '33.4.1.c', '33.5.10.a', '33.5.23.a', '33.6.2.a', '33.6.4.a', '33.6.4.b', '33.7.10.a', '33.7.23.a', '33.8.2.a', '33.8.4.a', '33.8.4.b', '33.9.10.a', '33.9.23.a', '34.10.13.a', '34.10.13.b', '34.10.33.a', '34.10.33.b', '34.10.9.a', '34.10.9.b', '34.12.13.a', '34.12.13.b', '34.12.33.b', '34.12.9.a', '34.12.9.b', '34.2.13.a', '34.2.13.b', '34.2.9.a', '34.4.1.a', '34.4.1.b', '34.4.13.a', '34.4.13.b', '34.4.9.a', '34.4.9.b', '34.5.3.a', '34.5.3.b', '34.6.13.a', '34.6.13.b', '34.6.33.a', '34.6.33.b', '34.6.9.a', '34.6.9.b', '34.7.3.a', '34.7.3.b', '34.8.13.a', '34.8.13.b', '34.8.33.a', '34.8.33.b', '34.8.9.a', '34.8.9.b', '35.4.1.a', '35.4.1.c', '37.10.36.a', '37.12.36.a', '37.2.36.a', '37.3.2.a', '37.5.8.a', '37.6.36.a', '37.7.8.a', '37.8.36.a', '38.10.7.a', '38.10.7.b', '38.11.27.a', '38.11.37.a', '38.12.7.a', '38.12.7.b', '38.2.7.a', '38.2.7.b', '38.3.27.a', '38.3.37.a', '38.4.1.c', '38.4.7.a', '38.4.7.b', '38.4.7.c', '38.5.37.a', '38.6.7.a', '38.6.7.b', '38.7.27.a', '38.7.37.a', '38.8.7.a', '38.8.7.b', '38.9.27.a', '38.9.37.a', '39.10.16.a', '39.10.16.b', '39.10.25.a', '39.10.25.b', '39.10.4.a', '39.10.4.b', '39.11.14.a', '39.11.17.a', '39.11.17.b', '39.11.29.a', '39.11.29.b', '39.12.25.a', '39.12.25.b', '39.12.4.a', '39.12.4.b', '39.2.16.a', '39.3.14.a', '39.3.17.a', '39.3.17.b', '39.3.29.a', '39.3.29.b', '39.4.1.a', '39.4.1.c', '39.4.16.a', '39.4.16.b', '39.4.16.c', '39.4.2.a', '39.4.4.a', '39.4.4.b', '39.4.4.c', '39.5.14.a', '39.5.17.a', '39.5.17.b', '39.5.29.a', '39.5.29.b', '39.6.16.a', '39.6.16.b', '39.6.2.a', '39.6.2.b', '39.6.25.a', '39.6.25.b', '39.6.4.a', '39.6.4.b', '39.7.14.a', '39.7.17.a', '39.7.29.a', '39.8.16.a', '39.8.16.b', '39.8.25.a', '39.8.25.b', '39.8.4.a', '39.8.4.b', '39.9.14.a', '39.9.17.a', '39.9.17.b', '40.10.21.a', '40.10.9.a', '40.11.17.a', '40.11.17.b', '40.11.19.a', '40.11.19.b', '40.11.19.c', '40.12.21.a', '40.5.17.a', '40.5.17.b', '40.5.17.c', '40.5.19.a', '40.5.19.b', '40.5.19.c', '40.6.21.a', '40.6.9.a', '40.7.17.a', '40.7.17.b', '40.7.19.a', '40.7.19.b', '40.7.19.c', '40.8.21.a', '40.9.17.a', '40.9.17.b', '40.9.19.a', '40.9.19.b', '40.9.19.c', '41.12.1.a', '43.4.9.a', '44.4.1.a', '45.4.19.b', '45.4.2.a', '47.11.46.c', '47.3.5.a', '47.8.1.a', '48.10.47.c', '48.11.17.e', '48.11.31.b', '48.12.47.b', '48.4.47.b', '48.5.31.a', '48.6.47.d', '48.7.31.b', '48.8.47.c', '48.9.31.c', '49.10.18.a', '49.10.18.b', '49.10.18.c', '49.10.18.d', '49.10.18.e', '49.10.18.g', '49.12.1.g', '49.12.18.d', '49.12.18.e', '49.12.18.i', '49.12.18.j', '49.3.3.a', '49.5.19.a', '49.5.19.b', '49.5.19.c', '49.7.19.c', '49.7.19.e', '49.8.1.e', '49.8.1.f', '49.9.19.b', '49.9.19.c', '49.9.19.d']


def check_data_for_Gamma1(max_level, max_weight, start_level=1, start_weight=1):
    C = CompMF('/mnt/data/stromberg/modforms-db')
    D = MongoMF()
    args = []
    for N in xrange(start_level, max_level+1):
        for k in [w for w in xrange(start_weight, max_weight+1)]:
            for r in D._mongodb['webmodformspace'].find({'level':N,'weight':k,'version':float(1.3)}):
                args.append(r['space_label'])
    return args

def spaces_from_form_labels(labels):
    res = []
    from lmfdb.modular_forms.elliptic_modular_forms.backend.emf_utils import parse_newform_label
  
    for label in labels:
        N,k,i,d=parse_newform_label(label)
        space_label = "{0}.{1}.{2}".format(N,k,i)
        if space_label not in res:
            res.append(space_label)
    res.sort()
    return res

def check_spaces_for_recomputation(args,do_recompute=False,**kwds):
    recompute_spaces = []
    for label in args:
        S = WebModFormSpace_computing(label, recompute=False)
        if check_one_space(S):
            continue
        recompute_spaces.append(S.space_label)
    print "Need to recompute {0} spaces!".format(len(recompute_spaces))
    ## Then do the recomputations. We can even try with parallell...
    if do_recompute:
        pool = Pool(processes=kwds.get('ncpus',1))
        chunksize=kwds.get('chunksize',20)
        results = pool.imap_unordered(recompute_space_completely,recompute_spaces,chunksize)
        return list(results)
    return recompute_spaces

# def check_forms_for_recomputation(args,do_recompute=False,**kwds):
#     recompute_spaces = []
#     for label in args:
#         F = WebNewForm_computing(label, recompute=False)
#         if check_one_form(S):
#             continue
#         recompute_forms.append(S.space_label)
#     print "Need to recompute {0} spaces!".format(len(recompute_spaces))
#     ## Then do the recomputations. We can even try with parallell...
#     if do_recompute:
#         pool = Pool(processes=kwds.get('ncpus',1))
#         chunksize=kwds.get('chunksize',20)
#         results = pool.imap_unordered(recompute_form_completely,recompute_forms,chunksize)
#         return list(results)
#     return recompute_forms


def check_one_space(S):
    success = True
    if not check_orbits(S):
        wmf_logger.info("space failed orbit check!")
        success= False
    if not check_if_updated(S):
        wmf_logger.info("space is not updated!")
        success= False        
    if not check_deligne(S):
        wmf_logger.info("space does not satisfy deligne's bound!")
        success= False
    if not check_coefficients(S):
        wmf_logger.info("a form in the space does not have correct coefficients!")
        success = False
    return success
    
def check_deligne(S):
    from sage.all import RealField
    for f in S.hecke_orbits.values():
        if not check_deligne_one_form(f):
            return False
    return True

def check_deligne_one_form(f):
    if f.dimension==0 and dimension_new_cusp_forms(f.character.sage_character,f.weight)==0:
        return True
    if f.max_cn()<2:
        return False
    for p in prime_range(f.max_cn()):
        try:
            cp = f.coefficient(p)
        except StopIteration:
            wmf_logger.info("Newform does not have coefficient {0}".format(p))
            return False
        if cp.norm()==0:
            continue
        t = 1000
        if cp.parent() <> QQ:
            prec_start = ceil(RR(abs(cp.norm())).log()/RR(p).log()/53.0)+1
            for mul_prec in range(prec_start,prec_start+20):
                RF = RealField(mul_prec*53)
                norm = RF(p)**((RF(f.weight)-RF(1))/RF(2))
                l = [x/norm for x in cp.complex_embeddings(53*mul_prec)]
                err = abs(sum(l) - cp.trace()/norm)
                if  err < 1e-8:
                    ### arbitrary test to ensure that the precision is sufficient
                    ###  Should be checked in theory...
                    t = max([abs(x) for x in l])
                    break
        else:  ## for a rational form 53 bits of precision should be ok...
            t = RR(abs(cp))/RR(p)**((f.weight-1.0)/2.0)
        if abs(t) > 2.0:
            wmf_logger.critical("The aps in the coefficients are incorrect for {0}. We got c({1})/{1}^(k-1)/2)={2} for c({1})={3} Please check!".format(f.hecke_orbit_label,p,t,cp))
            return False
    return True
    
        
def check_if_updated(S):
    r"""
    Return False if S or one of its orbits have not updated from db / fs
    """
    if not S.has_updated_from_db() and S.has_updated_from_fs():
        return False
    for f in S.hecke_orbits.values():
        if not f.has_updated_from_fs() and f.has_updated_from_db():
            return False
    return True
                
def check_orbits(S): ### This should essentially check more than Drew's check
    r"""
    Check consistency of orbits in the space S.
    Return True if the orbits are all there and have correct dimensions.
    """
    from sage.all import dimension_new_cusp_forms
    dim_true = dimension_new_cusp_forms(S.character.sage_character,S.weight)
    if S.dimension_new_cusp_forms != dim_true:
        S.set_dimensions()
    if S.dimension_new_cusp_forms != dim_true:
        return False
    dim_here = 0
    try:
        for label in S.hecke_orbits:
            F = S.hecke_orbits[label]
            FF = WebNewForm_computing(F.hecke_orbit_label)
            fact = FF.as_factor()
            if fact.character() != F.character.sage_character:
                wmf_logger.critical("The character is wrong!")
                return False
            if fact.dimension() != F.dimension:
                F.set_dimensions()
            dim_here += F.dimension
        if dim_here != dim_true:
            wmf_logger.info("The sum of dimensions doesn't add up!")
            return False
    except Exception as e:
        wmf_logger.critical("The orbit check failed because of {0}".format(e))
        return False
    return True
def check_coefficients(S):
    for label in S.hecke_orbits:
        f = S.hecke_orbits[label]
        if not check_coefficient_of_form(f):
            return False
    return True

def check_coefficient_of_form(F,nrange=[]):
    r"""
    Compare the coefficients of self with those of the form F.
    """

    if nrange == []:
        nrange = range(2,max(F.character.character.conductor()+1,F.max_cn()))
    if hasattr(F,'as_factor'):
        f = F.as_factor()
    else:
        f = WebNewForm_computing(F.hecke_orbit_label).as_factor()
    f = f.q_eigenform(max(nrange)+1,names='a')
    for i in nrange:
        c1 = f.padded_list()[i]
        c2 = F.coefficient(i)
        if c1 == c2:
            continue
        if hasattr(c1,'norm'):
            if c1.norm() == c2.norm() and c1.trace() == c2.trace():
                wmf_logger.debug("Norms and traces are equal so might be isomorphic parents!")
        wmf_logger.debug("Coefficients {0} differs!".format(i))
        return False
    return True

@parallel(ncpus=8)
def recompute_space_completely(label, path='/mnt/data/stromberg/modforms-db'):
    C = CheckingDB(path)
    D = MongoMF()
    from lmfdb.modular_forms.elliptic_modular_forms.backend.emf_utils import parse_space_label
    if hasattr(label,'space_label'):
        label = label.space_label
    N,k,ci = parse_space_label(label)
    C.compute_and_insert_one_space(N,k,ci)
    C.check_record(N,k,ci,check_content=True,recheck=True)
    cid = D.register_computation(level=N,weight=k,cchi=ci,typec='wmf')
    S = WebModFormSpace_computing(N,k,ci, recompute=True, update_from_db=False)
    S.save_to_db()
    D.register_computation_closed(cid)
