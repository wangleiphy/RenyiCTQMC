from numpy import zeros  , log 
from bitops import btest, bget, bset, bclear, fermionParity

def buildH(Kmat, V):
    '''
    -t*hopping +V*(n-0.5) (n-0.5)
    '''
    
    Nsites = Kmat.shape[0]
    Nstates = 1 << Nsites
    H = zeros((Nstates, Nstates), float)

    # loop over all basis states
    for state in range(Nstates):
        #print bin(state)

        for s in range(Nsites):
            for j in range(Kmat.indptr[s], Kmat.indptr[s+1]):
                sj = Kmat.indices[j]  # a neighboring site sj 

                #interaction , if to avoid double counting 
                if sj > s: H[state, state] += V* (bget(state, s) -0.5) *  (bget(state, sj) -0.5)

                #hopping from site s to site sj 
                if btest(state, s):
                    parity1 = fermionParity(state, s)
                    state2 = bclear(state,s)
                    if not btest(state2,sj):
                        parity2 = parity1 ^ fermionParity(state2,sj)
                        state3=bset(state2,sj)

                        H[state, state3] += Kmat.data[j] * (1-2*parity2)
    
    return H 

def S2(Nsite, LA, beta, w, v):
    '''
    compute renyi EE S2 
    '''

    weights = exp(-beta * w)
    Z = weights.sum()  
    weights = weights/Z 

    rho = dot(dot(v , diag(weights)  ) , v.transpose().conjugate())

    #print rho 

    Nstates = len(w)

    rhoA = zeros((1<<LA, 1<<LA), float)
    # loop over all basis states
    for statei in range(Nstates):
        for statej in range(Nstates): # this is less efficient 

            iB = int("{0:b}".format(statei).zfill(Nsite)[:-LA], 2)
            jB = int("{0:b}".format(statej).zfill(Nsite)[:-LA], 2)
            if (iB!=jB): continue

            #state index first translate to binarty (with fixed site Nsite), we then just take the lowest LA bits then turn it to an interger
            iA = int("{0:b}".format(statei).zfill(Nsite)[-LA:], 2)
            jA = int("{0:b}".format(statej).zfill(Nsite)[-LA:], 2)
    
            rhoA[iA, jA] += rho[statei, statej]
    
    #print rhoA 
    return  -log( dot(rhoA, rhoA).trace() )


if __name__=='__main__':
    from numpy import array , zeros , linspace , dot , exp , diag 
    import scipy.sparse as sps 
    from numpy.linalg import eigh 
    import sys 

    L = 8
    V = 1.0

    Thop = 1.0
    Kmat = zeros((L, L),float)
    for s in range(L-1):   
        Kmat[s, s+1] = -Thop 
        Kmat[s+1, s] = -Thop 

    #print Kmat 

    Kmat = sps.csr_matrix(Kmat)

    Hmat = buildH(Kmat, V)
    #print Hmat 

    w, v = eigh(Hmat)

    beta = 1.

    for LA in range(1,L):
        print LA, S2(L, LA, beta, w, v)



