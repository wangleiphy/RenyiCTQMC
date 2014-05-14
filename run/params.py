import os.path 
import subprocess
import re , sys 
from math import cos , sin , pi 

def writeParameterFile(fname,parms):
    """ This function writes a text input file for simple ALPS applications like DMFT
    
        The arguments are:
        
          filename: the name of the parameter file to be written
          parms: the parameter dict
    """
    f = file(fname,'w')
    for key in parms:
      value = parms[key]
      if type(value) == str:
        f.write(str(key)+' = "' + value + '"\n')
      else:
        f.write(str(key)+' = ' + str(value) + '\n')
    f.close()
    return fname

'''
write parameters for main 
'''

def params(lattice, L, W, NA0, NA1, V=1.0, T = 0.1, Maxorder = 2048, Ntau=500, Nscratch=5000, SWEEPS=1000000, THERMALIZATION=100000, NSKIP=200, Add=0.1, Remove = 0.2, ZtoW=0.1, WtoZ =0.1, eta=1., folder='../data/', textoutput=0):
    
    key = lattice.replace(' ','') 
    key += 'L' + str(L)\
           +'W' + str(W)\
           +'NA0' + str(NA0)\
           +'NA1' + str(NA1)\
           +'V'+str(V)\
           +'T' + str(T)\
           +'Ntau'+str(Ntau)\
           +'MAXORDER'+ str(Maxorder)\
           +'Nscratch'+str(Nscratch)\
           +'Therm'+str(THERMALIZATION)\
           +'Sweeps'+str(SWEEPS) \
           +'Skip'+str(NSKIP)\
           +'Add'+str(Add)\
           +'Remove'+str(Remove)\
           +'ZtoW'+str(ZtoW)\
           +'WtoZ'+str(WtoZ)\
           +'eta'+str(eta)\

    totalprob = Add + Remove+ ZtoW + WtoZ

    if (abs(totalprob-1.)> 1E-8):
        print 'summation of probs ', totalprob , ' are you sure ?'
        sys.exit(1)

    inputname = '../jobs/'+ key +'.in'
    outputname = folder + key +'.dat'
    
    parms ={'LATTICE_LIBRARY' : "../input/mylattices.xml"
            # above we should not change 

            ,'LATTICE'  : lattice
            ,'filename' : outputname
            ,'textoutput' :textoutput 
            ,'L'  : L 
            ,'W'  : W
            ,'NA0'  : NA0
            ,'NA1'  : NA1
            ,'FLAVORS' : 1
            ,'MAX_ORDER' : Maxorder
            ,'N_TAU' : Ntau
            ,'V' : V
            ,'TEMPERATURE' : T
            ,'RECALC_PERIOD' : Nscratch 
            ,'THERMALIZATION' : THERMALIZATION
            ,'SWEEPS' : SWEEPS 
            ,'MEASUREMENT_PERIOD'  : NSKIP
            ,'Add': Add
            ,'Remove': Remove
            ,'ZtoW': ZtoW
            ,'WtoZ': WtoZ 
            ,'eta' : eta 
            }

    writeParameterFile(inputname, parms)

    return inputname 
