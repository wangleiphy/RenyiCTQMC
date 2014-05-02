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

def S2(Kmat, V, LA, beta):
    '''
    compute renyi EE S2 
    '''
    
    Nsite = Kmat.shape[0]
    Hmat = buildH(Kmat, V)

    w, v = eigh(Hmat)

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
    from numpy import array , zeros , linspace , dot , exp , diag , arange 
    import scipy.sparse as sps 
    from numpy.linalg import eigh 
    import sys 

    L = 8
    Thop = 1.0

    #Kinetic energy matrix 
    Kmat = zeros((L, L),float)
    for s in range(L):   
        Kmat[s, (s+1)%L] = -Thop 
        Kmat[(s+1)%L, s] = -Thop 
    Kmat = sps.csr_matrix(Kmat)

     
    beta = 2.
    LA = L/2
    Snonint = S2(Kmat, 0.0, LA, beta)

    for V in arange(0.1, 1.1, 0.1):
        print V, S2(Kmat, V, LA, beta) - Snonint  , Snonint  


