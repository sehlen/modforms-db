


from sage.modular.modsym.space import is_ModularSymbolsSpace
from sage.all import (ModularSymbols, DirichletGroup, trivial_character,
                      dimension_new_cusp_forms,
                      save, load,
                      cputime,
                      fork, parallel,
                      Integer,
                      prime_range, prime_divisors,
                      version,
                      Sequence,
                      cached_function,
                      Integer)

from filesdb import rangify, FilenamesMFDBLoading




    
 

class ComputeMFData(object):
    r"""
    Use Williams methods to compute tables of modular forms data
    """
    def __init__(self,db,compute_ambient=False):
        r"""
        compute_ambient = True if you want to compute ambient spaces when missing
        """
        if isinstance(db,str):
            db = FilenamesMFDBLoading(db)
        self._compute_ambient = compute_ambient
        self._files = db  ## Should be instance of, e.g. FilenamesMFDB
        #self._collection = self._db

    def files(self):
        return self._files
    # Compute objects
    
    @fork    
    def compute_ambient_space(self,N, k, i,**kwds):#
        r"""
        Compute the ambient space M(N,k,i) and save as file. 
        """
        if i == 'all':
            G = DirichletGroup(N).galois_orbits()
            sgn = (-1)**k
            for j, g in enumerate(G):
                if g[0](-1) == sgn:
                    self.compute_ambient_space(N,k,j,**kwds)
            return
        if i == 'quadratic':
            G = DirichletGroup(N).galois_orbits()
            sgn = (-1)**k
            for j, g in enumerate(G):
                if g[0](-1) == sgn and g[0].order()==2:
                    self.compute_ambient_space(N,k,j,**kwds)
            return
        filename = self.files().ambient(N, k, i)
        if self.files().path_exists(filename):
            return
        eps = character(N, i)
        t = cputime()
        M = kwds.get('M',None)
        if M ==  None or M.sign()<>N or M.weight()<>k or M.level()<>N or M.character()<> eps or not is_ModularSymbolsSpace(M):
            t = cputime()
            M = ModularSymbols(eps, weight=k, sign=1)
            tm = cputime(t)
        else:
            tm = kwds.get('tm')
        self.files().save_ambient_space(M,i)
        #save(M, filename)
        meta = {'cputime':tm, 'dim':M.dimension(), 'M':str(M), 'version':sage.version.version}
        save(meta, self.files().meta(filename))
        if verbose>0:
            print "save {0} to {1}".format(meta,filename)


    def compute_ambient_spaces(self,Nrange, krange, irange, ncpu=1,**kwds):
        r"""
        Compute a list of ambient spaces in parallel. 
        """ 
        @parallel(ncpu)
        def f(N,k,i):
            self.compute_ambient_space(N,k,i,**kwds)

        v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
        for X in f(v):
            print X
   
       
    @fork    
    def compute_decompositions(self,N, k, i,verbose=0,**kwds):
        if i == 'all':
            G = DirichletGroup(N).galois_orbits()
            sgn = (-1)**k
            numo = 0
            for j, g in enumerate(G):
                if g[0](-1) == sgn:
                    numo+=self.compute_decompositions(N,k,j,**kwds)
            return numo

        if i == 'quadratic':
            G = DirichletGroup(N).galois_orbits()
            sgn = (-1)**k
            numo = 0
            for j, g in enumerate(G):
                if g[0](-1) == sgn and g[0].order()==2:
                    num0+=self.compute_decompositions(N,k,j,**kwds)
            return numo

        filename = self.files().ambient(N, k, i)        
        if verbose>0:
            print "filename=",filename
        if not self.files().path_exists(filename):
            print "Ambient space ({0},{1},{2}) not computed. filename={3}".format(N,k,i,filename)
            #return
            self.compute_ambient_space(N, k, i,**kwds)
        if not self.files().path_exists(filename):
            return 0
        eps = DirichletGroup(N).galois_orbits()[i][0]
        # check if the factor already exists by checking for 
        if verbose>0:
            print "check if path exists {0}".format(self.files().factor_eigen_nonzero(N, k, i, 0))
        if self.files().path_exists(self.files().factor_eigen_nonzero(N, k, i, 0)):
            return self.files().number_of_known_factors(N,k,i) 
        t = cputime()
        M = self.files().load_ambient_space(N, k, i)
        if verbose>0:
            print "M=",M
        if kwds.get('D') is None:
            D = M.cuspidal_subspace().new_subspace().decomposition()
        else:
            D = kwds['D']
        if verbose>0:
            print "D=",D
        for d in range(len(D)):
            self.files().factor(N,k,i,d,makedir=True)
            f = self.files().factor_basis_matrix(N, k, i, d)
            if self.files().path_exists(f):
                continue
            A = D[d]
            B  = A.free_module().basis_matrix()
            Bd = A.dual_free_module().basis_matrix()
            v  = A.dual_eigenvector(names='a', lift=False)    # vector over number field
            nz = A._eigen_nonzero()
            name = self.files().factor_basis_matrix(N, k, i, d)
            save(B, name)
            save(Bd, self.files().factor_dual_basis_matrix(N, k, i, d))
            save(v, self.files().factor_dual_eigenvector(N, k, i, d))
            save(nz, self.files().factor_eigen_nonzero(N, k, i, d))

        tm = cputime(t)
        meta = {'cputime':tm, 'number':len(D), 'version':sage.version.version}
        save(meta, self.files().decomp_meta(N, k, i))
        return len(D)
    
    def compute_decomposition_ranges(self,Nrange, krange, irange, ncpu,verbose=0):
        @parallel(ncpu)
        def f(N,k,i):
            self.compute_decompositions(N,k,i,verbose=verbose)

        v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
        for X in f(v):
            print X

    def compute_ambient_space_ranges(self,Nrange, krange, irange, ncpu):
        @parallel(ncpu)
        def f(N,k,i):
            self.files().compute_ambient_space(N,k,i)
        v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
        for X in f(v):
            print X
    # atkin_lehner
    @fork    
    def compute_atkin_lehner(self,N, k, i,M=None,m=None,**kwds):
        verbose = kwds.get('verbose',0)
        filename = self.files().ambient(N, k, i)
        if not self.files().path_exists(filename):
            print "Ambient (%s,%s,%s) space not computed. Filename=%s "%(N,k,i,filename)
            return -1
            #compute_ambient_space(N, k, i)
        if verbose>0:
            print "computing atkin-lehner for (%s,%s,%s)"%(N,k,i)
        if m is None:
            m = self.files().number_of_known_factors(N, k, i)
        if M is None:
            M = self.files().load_ambient_space(N, k, i)
        for d in range(m):
            atkin_lehner_file = self.files().factor_atkin_lehner(N, k, i, d, False)
            if self.files().path_exists(atkin_lehner_file):
                if verbose>0:
                    print "skipping computing atkin_lehner for (%s,%s,%s,%s) since it already exists"%(N,k,i,d)
                # already done
                continue
            # compute atkin_lehner
            if verbose>0:
                print "computing atkin_lehner for (%s,%s,%s,%s)"%(N,k,i,d)
            t = cputime()
            A = self.files().load_factor(N, k, i, d, M)
            #print "A=",A
            #print "character=",A.character()
            #print "character order=",A.character().order()
            al = ' '.join(['+' if a > 0 else '-' for a in atkin_lehner_signs(A)])
            #print al
            open(atkin_lehner_file, 'w').write(al)
            tm = cputime(t)
            meta = {'cputime':tm, 'version':sage.version.version}
            save(meta, self.files().meta(atkin_lehner_file))

    # aplists
    @fork    
    def compute_aplists(self,N, k, i, *args):
        if i == 'all':
            G = DirichletGroup(N).galois_orbits()
            sgn = (-1)**k
            for j, g in enumerate(G):
                if g[0](-1) == sgn:
                    self.compute_aplists(N,k,j,*args)
            return

        if i == 'quadratic':
            G = DirichletGroup(N).galois_orbits()
            sgn = (-1)**k
            for j, g in enumerate(G):
                if g[0](-1) == sgn and g[0].order()==2:
                    self.compute_aplists(N,k,j,*args)
            return

        if len(args) == 0:
            args = (100, )

        filename = self.files().ambient(N, k, i)
        if not self.files().path_exists(filename):
            print "Ambient ({0},{1},{2}) space not computed. Filename:{3}".format(N,k,i,filename)
            return -1
            #compute_ambient_space(N, k, i)
        if verbose>0:
            print "computing aplists for (%s,%s,%s)"%(N,k,i)

        m = self.files().number_of_known_factors(N, k, i)
        if verbose>0:
            print "m=",m
        if m == 0:
            # nothing to do
            return

        M = self.files().load_ambient_space(N, k, i)
        if verbose>0:
            print "M,m=",M,m
        for d in range(m):
            aplist_file = self.files().factor_aplist(N, k, i, d, False, *args)
            if self.files().path_exists(aplist_file):
                if verbose>0:
                    print "skipping computing aplist(%s) for (%s,%s,%s,%s) since it already exists"%(args, N,k,i,d)
                # already done
                continue

            # compute aplist
            if verbose>0:
                print "computing aplist(%s) for (%s,%s,%s,%s)"%(args, N,k,i,d)
            t = cputime()
            A = self.files().load_factor(N, k, i, d, M)
            aplist, _ = A.compact_system_of_eigenvalues(prime_range(*args), 'a')
            #print aplist, aplist_file
            save(aplist, aplist_file)
            tm = cputime(t)
            meta = {'cputime':tm, 'version':sage.version.version}
            if verbose>0:
                print "meta=",meta
            save(meta, self.files().meta(aplist_file))

    def compute_aplists_ranges(self,Nrange, krange, irange, ncpu, *args):
        @parallel(ncpu)
        def f(N,k,i):
            self.compute_aplists(N,k,i,*args)

        v = [(N,k,i) for N in rangify(Nrange) for k in rangify(krange) for i in rangify(irange)]
        for X in f(v):
            print X

    def phase1_goals(self,stages, fields=None):
        """
          Trivial character S_k(G0(N)):
             0. k=2   and N<=4096= 2^12
             1. k<=16 and N<=512 = 2^9
             2. k<=32 and N<=32  = 2^5
          Nontrivial character S_k(N, chi):
             3. k<=16, N<=128 = 2^7, all chi quadratic
             4. k=2,   N<=128 = 2^7, all chi!=1, up to Galois orbit
        """
        if 0 in stages:
            for X in self.files().find_missing(range(1,4096+1), [2], 0, fields=fields):
                yield X

        if 1 in stages:
            for X in self.files().find_missing(range(1,513),
                                            range(2,17,2), 0, fields=fields):
                yield X

        if 2 in stages:
            for X in self.files().find_missing(range(1,33),
                                            range(2,33,2), 0, fields=fields):
                yield X

        if 3 in stages:
            for X in self.files().find_missing(range(1,129),
                                     range(2,17), 'quadratic', fields=fields):
                yield X

        if 4 in stages:
            for X in self.files().find_missing(range(1,129), 2, 'all', fields=fields):
                yield X


    def phase1_goals0(self,stages, fields=None):
        v = []
        for X in phase1_goals(stages, fields):
            print X
            v.append(X)
        return v

    # suggestion for collection is: mfself.compute.nsql.missing.phase0
    def phase1_goals_self(self, stages, fields=None):
        for X in self.phase1_goals(stages, fields):
            print X
            self._collection.insert(X)

    def generate_computations_missing_M(self):
        v = []
        for X in self._collection.find_missing(fields="M"):
            N,k,i = X['space']
            v.append('compute_ambient_space(%s,%s,%s)'%(N,k,i))
        return list(sorted(set(v)))

    def generate_computations_missing_decomp(self):
        v = []
        for X in self._collection.find_missing(field="decomp"):
            N,k,i = X['space']
            v.append('compute_decompositions(%s,%s,%s)'%(N,k,i))
        return list(sorted(set(v)))        

    def generate_computations_missing_aplists(self):
        v = []
        for X in self._collection.find_missing(fields="other"):
            N,k,i = X['space']
            if 'aplist-00100' in X['other']:
                v.append('compute_aplists(%s, %s, %s, 100)'%(N,k,i))
            if 'aplist-00100-01000' in X['other']:
                v.append('compute_aplists(%s, %s, %s, 100, 1000)'%(N,k,i))
            if 'aplist-01000-10000' in X['other']:
                v.append('compute_aplists(%s, %s, %s, 1000, 10000)'%(N,k,i))
    ##         if 'charpoly-00100' in X['other']:
    ##             v.append('compute_charpolys(%s, %s, %s, 100)'%(N,k,i))
    ##         if 'charpoly-00100-01000' in X['other']:
    ##             v.append('compute_charpolys(%s, %s, %s, 100, 1000)'%(N,k,i))
    ##         if 'charpoly-01000-10000' in X['other']:
    ##             v.append('compute_charpolys(%s, %s, %s, 1000, 10000)'%(N,k,i))
        return list(sorted(set(v)))        



import sage.all
def parallel_eval(v, ncpu, do_fork=True):
    """
    INPUT:
    - v -- list of strings that can be eval'd in this scope.
    """
    @parallel(ncpu)
    def f(X):
        if do_fork:
            @fork
            def g():
                return eval(X)
            return g()
        else:
            return eval(X)
    for Z in f(v):
        print Z
        
    

def atkin_lehner_signs(A):
    N = A.level()
    return [A.atkin_lehner_operator(p).matrix()[0,0] for p in prime_divisors(N)]




    
