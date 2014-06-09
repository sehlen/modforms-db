Descriptions of databases with modular forms related information in our project 

=== File system based database ===

In the directory /data we store spaces of modular forms computed by the mfdb project.

The space of modular forms of level N, weight k and character chi (in the sage ordering, i.e. the first character in the chi'th Galois orbit in DirichletGroup(N) is stored in the directory

/N-k-chi/

In this directory there are files:  

ambient.sobj            # A modular symbols space descripbed in terms of a dictionary: 
                          {'basis':  [], 
                           'manin': [(.,.,.)...]
                           'eps':   [...]
                           'mod2term': [(.,.),...]
                           'rels': 
                           'space': (N,k,chi)}
ambient-meta.sobj       # meta data 
decomp-meta.sobj        # meta data about the decomposition 

[and optionally (only included to make some computations and re-computations of ap's easier):]
M.sobj                  # ModularSymbols space, i.e. sage: ModularSymbols(N,k,sign=1) 
M-meta.sobj

and subdirectories of the form:

000/
001/

each such directory is containing information about a Hecke orbit (ordered in the way that ModularSymbols.new_subspace().cuspidal_subspace().decomposition()  does it) given by an isrreducible modular symbols subspace. The data about it is given in the following files: 

B.sobj                          # basis matrix      -- d x n matrix over K [A.free_module().basis_matrix()]
Bd.sobj                         # dual_basis_matrix -- d x n matrix over K [A.dual_free_module().basis_matrix()]
v.sobj                          # dual_eigenvector -- degree n vector over relative extension of K [ A.dual_eigenvector(names='a', lift=False)]
nz.sobj                         # smallest i s.t the i-th entries of the entries of a basis for the dual vector space are not all 0 [A._eigen_nonzero()]


and optionally:

degree.txt                      # Information about the degree? 
aplist-A--B.sobj                # Contains a(p) for p in prime_range(A,B)
aplist-A--B-meta.sobj           # Meta information about the computed a(p)'s
atkin_lehner.txt                # Contains a list of "+" and "-" corresponding to eigenvalue 1 or -1 of W_p for p in prime_divisors(N)
                                # This is only given when either the character is trivial or when character is quadratic and p=N.

For more info see design.txt in the mfdb repo

=== Mongo database === 

On lmfdb.warwick.ac.uk:37010 we have the mongo database 'modularforms2' containing the following collections : 

Modular_symbols.[chunks/files]   # Contains  (ambient, i.e. not necessarily irreducible) modular symbols spaces
Newform_factors.[chunks/files]   # Contains newform factors (i.e. Hecke irreducible modular symbols spaces)
Atkin_Lehner.[chunks/files]      # Atkin-Lehner eigenvalues in the same format as above
ap.[chunks/files]                # Lists of ap's in the format (E,v) where 
vector_on_basis.[chunks/files]

 u'dimension_table',            # Dimension table. Format? 
 u'dimensions',                 # Dimension table. Format? 


+ other temporary / unrelated collections 



