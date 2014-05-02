from scipy.sparse import lil_matrix, csr_matrix
#from scipy.sparse.linalg import expm
from numpy import zeros , sqrt 
from bitops import bget, bflip

def buildH(Nsites, Hx, Jzmat, alpha=0.0):
    '''
    Jz mat is a csr sparse matrix 
    alpha 0.0 only transverse field, 1.0 only Jz couplings 
    '''

    Nstates = 1 << Nsites
    #H = lil_matrix((Nstates,Nstates),dtype=float)
    H = zeros((Nstates, Nstates), float)

    # loop over all basis states
    for i in range(Nstates):
        #print bin(i)
        for site in range(Nsites):
            for j in range(Jzmat.indptr[site], Jzmat.indptr[site+1]):
                sitej = Jzmat.indices[j]
                if sitej > site: 
                    # 0-> -1, 1-> 1 
                    H[i,i] +=  -(alpha) * Jzmat.data[j]* (2*bget(i, sitej)-1) \
                                                       * (2*bget(i, site) -1)
        
        # set the transverse field 
        for site in range(Nsites):
            j = bflip(i, site)
            H[i, j] = -Hx * (1.-alpha) 
    
    return H #.tocsr()

def energy(Nsites, Jzmat, i):
    res = 0.0
    for site in range(Nsites):
        for j in range(Jzmat.indptr[site], Jzmat.indptr[site+1]):
            # 0-> -1, 1-> 1 
            #print site, Jzmat.indices[j], -Jzmat.data[j], (2*bget(i, site) -1 ),  (2*bget(i, Jzmat.indices[j])-1)
            res +=  - Jzmat.data[j]* (2*bget(i, Jzmat.indices[j])-1) \
                                   * (2*bget(i, site) -1 )

    return res/2.0

if __name__=='__main__':
    from numpy import array , zeros , linspace 
    from numpy.linalg import eigh 
    import sys 
    from load_couplings import load_couplings
    
    inputname = '5site.input'
    Jzmat = load_couplings(inputname, comments='#')  
    Nsites = Jzmat.shape[0]
    
    Hx = -1.0 

    for alpha in linspace(0.0, 1.0, 11):
        Hmat = buildH(Nsites, Hx, Jzmat , alpha)
    #print Hmat #.todense()
        w, v = eigh(Hmat)
        print alpha, w[0], w[1], w[2], w[3], w[4]
    #print v[:,0], sqrt(1./len(w))

    sys.exit(1)
    
    i = int(str(11),2)
    print energy(Nsites, Jzmat, i)
    print "######"

    print i , bin(i)[2:].zfill(Nsites)
    for site in range(Nsites):
        print site, i, bget(i, site)

