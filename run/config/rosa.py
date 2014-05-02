import subprocess 
import time 
from numpy import arange 

#nickname = 'M4worm'
#nickname = 'fixT'
#nickname = 'scaleT'
#nickname = 'smartshift'
#nickname = 'smartM4worm'
#nickname = 'farthest'
#nickname = 'supercell'
#nickname = 'largestep'
#nickname = 'unequaltime'
#nickname = 'goldfine'
nickname = 'largeV'
measure_unequaltime = 0

Zupdate = 0.1
Wupdate = 0.1

ZtoW2 = 0.1
W2toZ = 0.1

W2toW4 = 0.06
W4toW2 = 0.12

ZtoW4 = 0.02
W4toZ = 0.12

#Nlist = [72]
#Tlist = [3./(4*6.)]
#expectedM2 = 0.13; expectedM4 = 0.035

#Nlist = [162]
#Tlist = [3./(4*9.)]
#expectedM2 = 0.08; expectedM4 = 0.012

#Nlist = [288]
#Tlist = [3/(4.*12.)]
#expectedM2 = 0.05; expectedM4 = 0.005

Nlist = [450]
Tlist = [3/(4.*15.)]
expectedM2 = 0.04; expectedM4 = 0.004

#Tlist = arange(0.44, 0.49, 0.01)
#Tlist = [0.51, 0.53, 0.55, 0.57, 0.59] 
#Vlist = arange(1.2, 1.5, 0.1)
#Vlist = [1.2, 1.3, 1.4]
#Vlist = arange(1.3, 1.52, 0.02)
#Vlist = [1.24, 1.26, 1.28]

Vlist = arange(1.4, 1.9, 0.1)
Nneighbors = 3
nmaxlist = [3] 

THERMALIZATION = 10**5
SWEEPS = 10**6
NSKIP = 100
Ntau = 500

wtime = '24:00:00'
tmin = 300
tmax = 600
ncores = 640
prog = '../bin/worm'

resfolder = '/scratch/rosa/lewang/spinlessdata/' + nickname  + '/'
h, m, s = [int(i) for i in wtime.split(':')]

Tlimit = max(3600*h + 60*m + s - int(tmax*2.) , 0)
prog += ' -i '+ str(tmin) + ' -a ' + str(tmax) + ' -T ' + str(Tlimit) 

def submitJob(bin,args,jobname,wtime,run=False,ncores=None, wait = None):

#SBATCH --ntasks-per-node=16
#SBATCH --mem=2048
            #prepare the job file 
            job='''#!/bin/bash
#SBATCH --ntasks=%g
#SBATCH --time=%s
#SBATCH --account=s395
#SBATCH --job-name=%s
#SBATCH --output=%s
#SBATCH --error=%s\n'''%(ncores,wtime,jobname,jobname+'.log',jobname+'.log')

            job +='aprun -n '+  str(ncores)+' '+ str(bin) + ' '
            for key, val in args.items():
                job += str(key) +' '+ str(val) + ' '
            
            #print job
            jobfile = open("jobfile", "w")
            jobfile.write("%s"%job)
            jobfile.close()
            
            #submit the job 
            if run:
                cmd = ['sbatch','jobfile']
            else:
                cmd = ['cat','jobfile']

            subprocess.check_call(cmd)
            time.sleep(0.05)

