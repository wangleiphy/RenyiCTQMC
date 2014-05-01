import pyalps 
import pyalps.plot
import matplotlib.pyplot as plt 

prefix = '../data/fulldiag/honeycomblatticeL2_W2_V0.5'
resultFiles = pyalps.getResultFiles(prefix=prefix)

print resultFiles 

data = pyalps.evaluateFulldiagVersusT(resultFiles,DELTA_T=0.1, T_MIN=0.1, T_MAX=1.)

print data 

#for s in pyalps.flatten(data):
#  plt.figure()
#  plt.title("Antiferromagnetic Heisenberg chain")
#  pyalps.plot.plot(s)
#plt.show()
