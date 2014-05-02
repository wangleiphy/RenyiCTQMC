from numpy import array, exp , diag , dot , zeros , log 
from scipy.linalg import eigh 

t = 1.
V = 1.
beta = 10.

#-t*hop +V*(n-0.5) (n-0.5)
Ham = array([[0.25*V, 0., 0., 0. ], 
             [0., -0.25*V, -t, 0.],
             [0., -t, -0.25*V, 0.],
             [0., 0.,  0., 0.25*V]]
           )

#print Ham

e, v = eigh(Ham) 

#print e 
weights = exp(-beta * e)
Z = weights.sum()  
weights = weights/Z 

rho = dot(dot(v , diag(weights)  ) , v.transpose().conjugate())

rhoA = zeros((2,2), float)
rhoA[0,0]= rho[0,0] + rho[1,1]
rhoA[0,1]= rho[0,2] + rho[1,3]
rhoA[1,0]= rho[2,0] + rho[3,1]
rhoA[1,1]= rho[2,2] + rho[3,3]

print rho 
print rhoA 

print -log( dot(rhoA, rhoA).trace() )


