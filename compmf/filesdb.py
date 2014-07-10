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
#################################################################
#                                                               #
# Database of Classical GL2 Holomorphic Modular Forms over Q    #
# (c) William Stein, 2012                                       #
#                                                               #
#################################################################
r"""

Database classes for storing data for (holomorphic, cuspidal) modular forms on subgroups of the modular group and integer weight.

In this file we define three classes:
Filenames                            -- basic class for storing objects in files
FilenamesMFDB(Filenames)             -- specialized for storing modular forms related objects
FilenamesMFDBLoading(FilenamesMFDB)  -- class for accessing and storing modular forms data

These classes are based on the file compute.py in the project: https://github.com/williamstein/mfdb
by William Stein.




"""
import os
from sage.all import cached_function,prime_range,DirichletGroup,dimension_new_cusp_forms,save,load,trivial_character,Sequence,ModularSymbols,prime_pi,nth_prime

import sqlite3
from stat import S_ISDIR

try:
    import paramiko
    has_paramiko = True
except ImportError:
    have_paramiko = False
    pass # print "Note that remote files are not supported without paramiko installed!"

from compmf import clogger

class Filenames(object):
    def __init__(self, datadir,host='',db_file='',username=''):
        r"""
        If host is left empty we work with local files, otherwise via SFTP
        """
        self._host = host
        self._sftp = None
        self._ssh = None
        if self._host <> '':
            try:
                self._ssh = paramiko.SSHClient()
            except ImportError:
                raise ValueError,"If paramiko is not installed you must use local files! You set host:{0}".format(self._host)
            self._ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self._ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
            self._ssh.connect(self._host, username=username, password='')
            self._sftp = self._ssh.open_sftp()
        if not self.path_exists(datadir):
            raise RuntimeError, "please create the data directory '%s'"%datadir
        self._data = datadir
        if db_file <> '':
            self._known_db_file = db_file
        elif self._host=='':
            self._known_db_file = os.path.join(self._data, 'known.sqlite3')
        else:
            self._known_db_file = "{0}:{1}/known.sqlite3".format(self._host,self._data)

    ##
    ## file interface functions
    ##
    def path_exists(self,f):
        r"""
        Check if the path f exists as a local file or on the remote host.
        """
        if self._sftp == None:
            return os.path.exists(f)
        else:
            return remote_exists(self._sftp,f)

    def make_path_name(self,f,*args):
        r"""
        Make a path name from a list of arguments.
        """
        if self._sftp==None:
            s = os.path.join(self._data,f)
        else:
            s = "{0}/{1}".format(self._data,f)
        for fs in args:
            s+= "/{0}".format(fs)
        return s

    def makedirs(self,f):
        r"""
        Make directories.
        """
        if self._sftp==None:
            os.makedirs(f)
        else:
            self._sftp.mkdir(f)
    def listdir(self,f):
        r"""
        List files in directory.
        """
        try:
            if self._sftp == None:
                return os.listdir(f)
            else:
                return self._sftp.listdir(f)
        except OSError:
            return []
    def isdir(self,path):
        r"""
        Check if path is a directory.
        """
        if self._sftp==None:
            return os.path.isdir(path)
        else:
            try:
                return S_ISDIR(self._sftp.stat(path).st_mode)
            except IOError:
                # Path does not exist, so by definition not a directory
                return False

    def delete_file(self,path):
        r"""
        Delete file at path.
        """
        if self._sftp == None:
            os.unlink(path)
        else:
            self._sftp.remove(path)

class FilenamesMFDB(Filenames):
    r"""
    Class for working with files representing modular forms spaces.


    NOTE: In the descriptions below we use M(N,k,i) to denote the space of modular forms
    on Gamma0(N), weight k and with character chi = DirichletGroup(N).galois_orbits()[i][0]
    and similarly for ModSym(N,k,i)
    """

    def __init__(self,datadir,host='',db_file='',username=''):
        r"""
        Initialize self.
        """
        super(self,FilenamesMFDB).__init(datadir,host=host,db_file=db_file,username=username)

    ###
    ### Functions which decide the formats for filenames of stored objects
    ###


    def meta(self, filename):
        r"""
        The filename of an associated meta-data object.
        """
        base, ext = os.path.splitext(filename)
        return base + '-meta.sobj'

    def space_name(self, N, k, i):
        r"""
        Directory name files related to M(N,k,i)
        """
        return '%05d-%03d-%03d'%(N,k,i)

    def space(self, N, k, i, makedir=False):
        r"""
        Return the full directory name of M(N,k,i) and create it if it doesn't exist.
        """
        f = self.make_path_name(self.space_name(N,k,i))
        if makedir and not self.path_exists(f):
            self.makedirs(f)
        return f

    def M(self, N, k, i, makedir=False):
        r"""
        The name of the file containing ModSym(N,k,i)
        """
        return self.make_path_name(self.space(N,k,i,makedir=makedir), 'M.sobj')

    def ambient(self, N, k, i, makedir=False):
        r"""
        The name of the file containing the modular symbols space M(N,k,i) stored in dictionary format.
        """
        return self.make_path_name(self.space(N,k,i,makedir=makedir), 'ambient.sobj')

    def decomp_meta(self, N, k, i):
        r"""
        The name of the metadata for the decomposition.
        """
        return self.make_path_name(self.space(N,k,i), 'decomp-meta.sobj')

    def factor(self, N, k, i, d, makedir=False):
        r"""
        The name of the file containing the newform factor nr. d of the space M(N,k,i)
        """
        f = self.make_path_name(self.space(N,k,i,makedir), '%03d'%d)
        if makedir and not self.path_exists(f):
            self.makedirs(f)
        return f


    #
    ## Functions for returning the names of files for different objects.
    ##

    def factor_basis_matrix(self, N, k, i, d):
        r"""
        The basis matrix of the newform factor nr. d of M(N,k,i)
        """
        return self.make_path_name(self.factor(N,k,i,d), 'B.sobj')

    def factor_dual_basis_matrix(self, N, k, i, d):
        r"""
        The dual basis matrix of the newform factor nr. d of M(N,k,i)
        """
        return self.make_path_name(self.factor(N,k,i,d), 'Bd.sobj')

    def factor_dual_eigenvector(self, N, k, i, d, makedir=True):
        r"""
        The dual eigenvector of the newform factor nr. d of M(N,k,i)
        """
        return self.make_path_name(self.factor(N,k,i,d,makedir=makedir), 'v.sobj')

    def factor_eigen_nonzero(self, N, k, i, d):
        r"""
        The smallest i s.t the i-th entries of the entries of a basis for the dual vector space are not all 0.
        """
        return self.make_path_name(self.factor(N,k,i,d), 'nz.sobj')

    def factor_aplist(self, N, k, i, d, makedir, *args):
        r"""
        The list of a(p)'s for the factor with range determined by the *args.
        """
        if len(args)==1:
            args = (0,args[0])
        a = '-'.join('%05d'%x for x in args)
        return self.make_path_name(self.factor(N,k,i,d,makedir), 'aplist-%s.sobj'%a)

    def factor_atkin_lehner(self, N, k, i, d, makedir):
        r"""
        The Atkin-Lehner eigenvalues of the factor.
        """
        return self.make_path_name(self.factor(N,k,i,d,makedir), 'atkin_lehner.txt')

    ##
    ## Quering the files directly
    ##

    def number_of_known_factors(self, N, k, i):
        r"""
        Returning the number of factors for which we have computed data.
        TODO: add check for empty files.
        """
        fn = self.space(N,k,i)
        return len([d for d in self.listdir(fn) if d.isdigit() and self.isdir(self.make_path_name(fn, d)) ])

    def known_levels(self):
        r"""
        The number of levels in the database.
        """
        res = []
        for Nki in self.listdir(self._data):
            z = Nki.split('-')
            if len(z) == 3:
                N, k, i = parse_Nki(Nki)
                if N not in res:
                    res.append(N)
        return res

    def find_known(self):
        """
        Return iterator of 5-tuples of Python ints, defined as follows:

        (N, k, i, newforms, maxp)
        (37, 2, 0, 2, 10000)

        Here N = level, k = weight, i = character, newforms = number of newforms,
        maxp = integer such that a_p is known for p<=maxp.

        If no newforms are known but there are newforms (they just
        haven't been computed), then newforms is set to -1.
        """
        for Nki in self.listdir(self._data):
            z = Nki.split('-')
            if len(z) == 3:
                N, k, i = parse_Nki(Nki)
                if k==1: # weight 1 not implemented
                    continue
                newforms = [x for x in self.listdir(self.make_path_name(self._data, Nki)) if x.isdigit()]
                if len(newforms) == 0:
                    # maybe nothing computed?
                    if i == 0:
                        # program around a bug in dimension_new_cusp_forms: Trac 12640
                        d = dimension_new_cusp_forms(N,k)
                    else:
                        chi = character(N, i)
                        d = dimension_new_cusp_forms(chi, k)
                    if d == 0:
                        # definitely no newforms
                        yield (N,k,i,0,0)
                    else:
                        # we just don't know the newforms yet
                        yield (N,k,i,-1,0)
                else:
                    maxp = None
                    for n in newforms:
                        v = set([])
                        this_maxp = 0
                        for X in self.listdir(self.make_path_name(self._data, Nki, n)):
                            if X.startswith('aplist') and 'meta' not in X:
                                args = [int(a) for a in X.rstrip('.sobj').split('-')[1:]]
                                v.update(prime_range(*args))
                                this_maxp = max(this_maxp, max(args))
                        if len(v) != len(prime_range(this_maxp)):
                            # something missing!
                            clogger.debug("data ranges are missing in the aplist data for %s"%Nki)
                            maxp = 100
                        else:
                            maxp = this_maxp if maxp is None else min(this_maxp, maxp)

                    yield (N,k,i,len(newforms),maxp)

    ##
    ## Since the file system can be quite slow we also have a sqlite database
    ##
    def update_known_db(self):
        r"""
        Create the sqlite database and iterate through the known files.
        """
        # 1. create the sqlite3 database
        # Unlink is only necessary with loal files
        if os.path.exists(self._known_db_file):
            os.unlink(self._known_db_file)
        db = sqlite3.connect(self._known_db_file)
        cursor = db.cursor()
        schema = 'CREATE TABLE "known" (N, k, i, newforms, maxp)'
        cursor.execute(schema)
        db.commit()

        # 2. iterate over known 5-tuples inserting them in db
        for t in self.find_known():
            #print t
            cursor.execute("INSERT INTO known VALUES(?,?,?,?,?)", t)

        db.commit()


    def known(self, query):
        r"""
        Return the result of an SQL query against the datbase.
        """
        query = query.replace('prec','maxp')

        # 1. open database
        db = sqlite3.connect(self._known_db_file)
        cursor = db.cursor()

        # 2. return result of query
        # really stupid and dangerous?
        if ';' in query:
            query = query.split(';')[0]

        cmd = 'SELECT * FROM known '
        if query.strip():
            cmd += 'WHERE %s'%query
        cmd += ' ORDER BY N, k, i'
        return cursor.execute(cmd)





    def find_missing(self, Nrange, krange, irange, fields=None):
        """
        Return generator of

        {'N':N, 'k':k, 'i':i, 'missing':missing, ...},

        where missing is in the intersection of fields and the
        following strings (or all strings if fields is None):

        'M', 'decomp',
        'aplist-00100',  'aplist-00100-01000',  'aplist-01000-10000',
        'charpoly-00100','charpoly-00100-01000','charpoly-01000-10000',
        'zeros',
        'leading',
        'atkin_lehner'

        If the string is not 'M' or 'decomp' then the meaning is that
        the indicated data is not complete for *all* newforms in this
        space, i.e., at least one is missing.  If 'decomp' is given,
        it means the decomp isn't complete (it could be partial).

        Spaces with dimension 0 are totally ignored.

        Note that not every missing is listed, just the *next* one that
        needs to be computed, e.g., if (11,2,0,'M') is output, then
        nothing else for (11,2,0,*) will be output.
        """
        if fields is None:
            fields = set(['M', 'decomp',
                          'aplist-00100',  'aplist-00100-01000',  'aplist-01000-10000',
                          'charpoly-00100','charpoly-00100-01000','charpoly-01000-10000',
                          'zeros',
                          'leading',
                          'atkin_lehner'])
        else:
            assert isinstance(fields, (list, tuple, str))
            fields = set(fields)

        space_params = set(self.listdir(self._data))
        for k in rangify(krange):
            for N in rangify(Nrange):
                for ch in rangify(irange):

                    if isinstance(ch, str):
                        CHI = list(enumerate(characters(N)))
                        if ch == 'quadratic':
                            CHI = [(i,chi) for i,chi in CHI if chi.order()==2]
                        elif ch == 'all':
                            pass
                        else:
                            raise ValueError
                    else:
                        try:
                            CHI = [(ch, character(N, ch))]
                        except IndexError:
                            CHI = []

                    for i, chi in CHI:
                        N,k,i = int(N), int(k), int(i)
                        obj = {'space':(N,k,i)}
                        dim = dim_new(chi, k)
                        if dim > 0:
                            fields0 = list(fields)
                            if 'M' in fields0:
                                fields0.remove('M')
                            if 'decomp' in fields0:
                                fields0.remove('decomp')
                            Nki = self.space_name(N,k,i)
                            d3 = self.make_path_name(self._data, Nki)
                            if Nki not in space_params or not self.path_exists(self.make_path_name(d3, 'M-meta.sobj')):
                                if 'M' in fields:
                                    obj2 = dict(obj)
                                    obj2['missing'] = 'M'
                                    yield obj2
                                break
                            newforms = []
                            for fname in self.listdir(d3):
                                if fname.isdigit():
                                    # directory containing data about a newforms
                                    d2 = self.make_path_name(d3, fname)
                                    deg = self.make_path_name(d2, 'degree.txt')
                                    if self.path_exists(deg):
                                        degree = eval(open(deg).read())
                                    else:
                                        B_file = self.make_path_name(d2, 'B.sobj')
                                        if not self.path_exists(B_file):
                                            degree = None
                                        else:
                                            degree = load(B_file).nrows()
                                            open(self.make_path_name(d2, 'degree.txt'),'w').write(str(degree))
                                    f = {'fname':fname, 'degree':degree,
                                         'other':set([x.split('.')[0] for x in self.listdir(d2)])}
                                    newforms.append(f)
                            degs = [f['degree'] for f in newforms]
                            if None in degs:
                                sum_deg = None
                            else:
                                sum_deg = sum(degs)
                            if 'decomp' in fields and (sum_deg is None or sum_deg != dim):
                                obj2 = dict(obj)
                                obj2['missing'] = 'decomp'
                                obj2['newforms'] = newforms
                                if sum_deg > dim:
                                    obj2['bug'] = 'sum of degrees (=%s) is too big (should be %s) -- internal consistency error!'%(sum_deg, dim)
                                yield obj2
                                break
                            missing = []
                            for other in fields0:
                                for j in range(len(newforms)):
                                    if other not in newforms[j]['other']:
                                        missing.append(other)
                                        break
                            if missing:
                                missing.sort()
                                obj2 = dict(obj)
                                obj2['missing'] = 'other'
                                obj2['other'] = missing
                                yield obj2




    def show_existing_files_db(self):
        r"""
        Display an overview of the existing files.
        """
        self.update_known_db()
        res = {}
        for t in self.find_known():
            N,k,i,newf,d = t
            if N not in res:
                res[N] = {}
            if k not in res[N]:
                res[N][k] = {}
            res[N][k][i]=newf,d
        Ns = res.keys(); Ns.sort()
        s=""
        for N in Ns:
            #s="{0:0>4} \t : \n".format(N)
            ks = res[N].keys(); ks.sort()
            #print "ks=",ks
            for k in ks:
                s+="{0:0>3} {1:0>3} \t ".format(N,k)
                i_s = res[N][k].keys(); i_s.sort()
                for i in i_s:
                    no,nap = res[N][k][i]
                    s+="{0}\t{1} \t {2:0>3}\t".format(i,no,nap)
                s+="\n"
        return s



################################################

class FilenamesMFDBLoading(FilenamesMFDB):
    r"""
    Class for loading / saving and converting modular forms objects.
    """
    def __init__(self,data):
        super(FilenamesMFDB,self).__init__(data)

    def save_ambient_space(self,M, i):
        r"""
        Save the ambient space.
        """
        N = M.level()
        k = M.weight()
        fname = self.ambient(N, k, i, makedir=True)
        if self.path_exists(fname):
            clogger.debug("%s already exists; not recreating"%fname)
            return
        clogger.debug("Creating  {0}".format(fname))
        clogger.debug("space {0}".format(M))
        ambient = ambient_to_dict(M, i)
        clogger.debug("as dict: {0}".format(ambient))
        save(ambient,fname)

        
    def load_M(self,N, k, i):
        r"""
        Load the ambient space as a modular symbols space.
        """
        return load(self.M(N, k, i, makedir=False))

    def convert_M_to_ambient(self,N, k, i):
        r"""
        Take an M.sobj file and create an ambient.sobj file containing the correspoding dictionary.
        """
        self.save_ambient_space(self.load_M(N,k,i), i)

    def convert_all_M_to_ambient(self):
        r"""
        Go through the database and convert all M.sobj files.
        """
        d = self._data

        for X in self.listdir(d):
            p = self.make_path_name(d, X)
            if self.isdir(p):
                f = set(self.listdir(p))
                if 'M.sobj' in f and 'ambient.sobj' not in f:
                    print X
                    #try:
                    self.convert_M_to_ambient(*parse_Nki(X))
                    #except:
                    #    clogger.debug("ERROR!")

    def delete_all_M_after_conversion(self):
        r"""
        Gto through the database and delete all converted M.sobj files.
        """
        d = self._data
        for X in self.listdir(d):
            p = self.make_path_name(d, X)
            if self._db.isdir(p):
                f = set(self.listdir(p))
                if 'M.sobj' in f and 'ambient.sobj' in f:
                    clogger.debug("Remove: \n {0}".format(X))
                    self.delete_file(self.make_path_name(p, 'M.sobj'))


    def load_ambient_space(self,N, k, i):
        r"""
        Load the ambient space as a modular forms space. If the M.sobj exist we load from that, otherwise we initialize from ambient.sobj.

        """
        fnameM = self.M(N, k, i, makedir=False)
        if self.path_exists(fnameM):
            return load(fnameM)
        fname = self.ambient(N, k, i, makedir=False)
        if self.path_exists(fname):
            ambient = load(fname)
            if not ambient is None:
                M = dict_to_ambient(ambient)
                #clogger.debug("Constructed space: {0}".format(M))
                return M
        s = "ambient space (%s,%s,%s) not yet computed"%(N,k,i)
        clogger.warning(s)
        raise ValueError,s

    def load_factor(self,N, k, i, d, M=None):
        r"""
        Load the factor nr. d in M(N,k,i)
        """
        import sage.modular.modsym.subspace
        if M is None:
            M = self.load_ambient_space(N, k, i)
        f = self.factor(N, k, i, d, makedir=False)
        if not self.path_exists(f):
            raise RuntimeError, "no such factor (%s,%s,%s,%s)"%(N,k,i,d)
        try:
            B = load(self.factor_basis_matrix(N, k, i, d))
            Bd = load(self.factor_dual_basis_matrix(N, k, i, d))
            v = load(self.factor_dual_eigenvector(N, k, i, d))
            nz = load(self.factor_eigen_nonzero(N, k, i, d))
        except IOError:
            raise RuntimeError,"Data is incomplete for factor ({0}) at {1}".format((N,k,i,d),f)
        if B._cache is None:
            B._cache = {}
        if Bd._cache is None:
            Bd._cache = {}
        B._cache['in_echelon_form'] = True
        Bd._cache['in_echelon_form'] = True
        # These two lines are scary, obviously, since they depend on
        # private attributes of modular symbols.
        A = sage.modular.modsym.subspace.ModularSymbolsSubspace(M, B.row_module(), Bd.row_module(), check=False)
        A._HeckeModule_free_module__dual_eigenvector = {('a',nz):(v,False)}
        A._is_simple = True
        A._HeckeModule_free_module__decomposition = {(None,True):Sequence([A], check=False)}
        return A

    def load_aps(self,N, k, i, d, ambient=None,numc=-1):
        r"""
        Load aps for a given factor. If numc > 0 we return all sets.
        """
        import sage.modular.modsym.subspace
        clogger.debug("Load_aps of {0}".format((N,k,i,d)))
        F = self.load_factor(N, k, i, d, ambient)
        factor_dir = self.factor(N, k, i, d, makedir=False)
        try:
            tmp = self.listdir(factor_dir)
        except OSError:
            return
        aplist_files = []; aplist_meta_files = []
        for fname in tmp:
            if "aplist" in fname:
                if "meta" in fname:
                    numap = int(fname.split("-")[2])
                    aplist_meta_files.append((numap,fname))
                else:
                    numap = int(fname.split("-")[2].split(".")[0])
                    aplist_files.append((numap,fname))
        if aplist_files == []:
            return
        if numc == 'max' or numc == -1: # Find max no. of coeffs.
            numap,fname = max(aplist_files)
        elif numc == 'min' or numc == 0: # Find min. no. of coeffs.
            numap,fname = min(aplist_files)
        else:
            numap,fname = aplist_files[0]
            for n,c in aplist_files[1:]:
                if c >= numc:
                    numap,fname = n,c
                    break
        clogger.debug("want numc={0} and have: {1}".format(numc,numap))
        if numc <> 'max' and numc<> 'min' and numc>0 and numap < numc:
            maxc = max(aplist_files)[0]
            print "aplist_files=",aplist_files
            raise ValueError,"We do not have {0} coefficients! At most: {1}".format(numc,maxc)
        metaname = fname.split(".")[0]+"-meta.sobj"
        try:
            E = load("{0}/{1}".format(factor_dir,fname))
            if isinstance(E,tuple):
                E = E[0]
            if prime_pi(numap) > E.nrows():
                clogger.debug("Have only {0} aps. Claim that we have {1}".format(E.nrows(),numap))
                ## Rename the file:
                new_prec = nth_prime(E.nrows()+1)-1
                save(E,self.factor_aplist(N,k,i,d,False,0,new_prec))
        except Exception as e:
            raise ValueError,"Could not load factor: {0}/{1}. Error:{2}".format(factor_dir,fname,e.message)
        try: 
            v = load("{0}/v.sobj".format(factor_dir))
        except Exception as e:
            raise ValueError,"Could not load factor: {0}/{1}. Error:{2}".format(factor_dir,'v.sobj',e.message)      
        try: 
            meta = load("{0}/{1}".format(factor_dir,metaname))
        except Exception as e:
            metaname1 = "aplist-00100-meta.sobj".format(numap)
            try:
                meta = load("{0}/{1}".format(factor_dir,metaname1))                
            except Exception as e:
                clogger.critical("Could not load factor: {0}/{1} or /{2}. Error:{3}".format(factor_dir,metaname,metaname1,e.message))
                meta = {}
            
        return E,v,meta


###### Helper routines ############



@cached_function
def characters(N):
    """
    Return representatives for the Galois orbits of Dirichlet characters of level N.
    """
    return [X[0] for X in DirichletGroup(N).galois_orbits()]

def character_to_int(eps):
    """
    Return integer corresponding to given character.
    """
    if eps.is_trivial():
        return 0
    N = eps.modulus()
    X = characters(N)
    try:
        return X.index(eps)
    except IndexError:
        # very unlikely -- would have to be some weird character
        # not got from character(N,i)
        for i, Y in enumerate(DirichletGroup(N).galois_orbits()):
            if X in Y:
                return i
        raise RuntimeError

def character(N, i):
    if i==0:
        return trivial_character(N)
    #clogger.debug("character(%s, %s)"%(N,i))
    return characters(N)[i]



def ambient_to_dict(M, i=None):
    """
    Data structure of dictionary that is created:

    - 'space': (N, k, i), the level, weight, and and integer
    that specifies character in terms of table; all Python ints
    - 'eps': (character) list of images of gens in QQ or cyclo
    field; this redundantly specifies the character
    - 'manin': (manin_gens) list of 2-tuples or 3-tuples of
    integers (all the Manin symbols)
    - 'basis': list of integers -- index into above list of Manin symbols
    - 'rels': relation matrix (manin_gens_to_basis) -- a sparse
    matrix over a field (QQ or cyclotomic)
    - 'mod2term': list of pairs (n, c), such that the ith element
    of the list is equivalent via 2-term relations only
    to c times the n-th basis Manin symbol; these
    """
    clogger.debug("Make ambient from M={0}, i={1}".format(M,i))        
    if i is None:
        i = character_to_int(M.character())
    return {'space':(int(M.level()), int(M.weight()), int(i)),
            'eps':list(M.character().values_on_gens()),
            'manin':[(t.i,t.u,t.v) for t in M._manin_generators],
            'basis':M._manin_basis,
            'rels':M._manin_gens_to_basis,
            'mod2term':M._mod2term}

def dict_to_ambient(modsym):
    r"""
    Convert a dictionary (as stored in ambient.sobj) to a modular symbols space.
    """
    #clogger.debug("Converting dict to ambient!")
    N,k,i = modsym['space']
    eps   = modsym['eps']
    manin = modsym['manin']
    basis = modsym['basis']
    rels  = modsym['rels']
    mod2term  = modsym['mod2term']
    F = rels.base_ring()
    if i == 0:
        eps = trivial_character(N)
    else:
        eps = DirichletGroup(N, F)(eps)

    from sage.modular.modsym.manin_symbols import ManinSymbolList, ManinSymbol
    manin_symbol_list = ManinSymbolList(k, manin)

    def custom_init(M):
        # reinitialize the list of Manin symbols with ours, which may be
        # ordered differently someday:
        syms = M.manin_symbols()
        ManinSymbolList.__init__(syms, k, manin_symbol_list)
        
        M._manin_generators = [ManinSymbol(syms, x) for x in manin]
        M._manin_basis = basis
        M._manin_gens_to_basis = rels
        M._mod2term = mod2term
        return M

    return ModularSymbols(eps, k, sign=1, custom_init=custom_init, use_cache=False)



def remote_exists(sftp, path):
    """os.path.exists for paramiko's SCP object
    """
    try:
        sftp.stat(path)
    except IOError, e:
        if e[0] == 2:
            return False
        raise
    else:
        return True


def parse_Nki(s):
    return tuple(int(x) for x in s.split('-'))


def rangify(v):
    return [v] if isinstance(v, (int, long, Integer, str)) else v
