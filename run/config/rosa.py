import subprocess 
import time 
from numpy import arange 

Add  = 0.15
Remove = 0.15

ZtoW = 0.35
WtoZ = 0.35
eta = 0.5

#latticename = 'open chain lattice'
#latticename = 'open honeycomb lattice'
#latticename = 'cylindrical honeycomb lattice'
#latticename = 'honeycomb lattice'
latticename = 'square lattice'
#######################################
nickname = 'ratiotrick'

Llist = [12]
Wlist = [12]
NAstep = 12

#Tlist = [0.5]
Tlist =[0.9, 1.0]
#Tlist = arange(0.6, 2.1, 0.1)
#Vlist = arange(0.5, 3.0, 0.5)
#Vlist = arange(2., 11., 1.)
Vlist = [2.]
#Vlist = [1.1, 1.2, 1.3, 1.4, 1.5]

NSKIP = 100
THERMALIZATION = 10**5
SWEEPS = 10**6
Nscratch = 1000
Ntau = 1000

wtime = '24:00:00'
tmin = 300
tmax = 600
ncores = 256 
prog = '../bin/ratiotrick'
#######################################

resfolder = '/scratch/rosa/lewang/renyidata/' + nickname  + '/'
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

