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

Llist = [4]
Wlist = [4]
NAstep = 4 

Tlist = [0.1]
#Tlist = arange(0.2, 1.1, 0.1)
Vlist = arange(0.2, 1.6, 0.2)
#Vlist = arange(0.5, 10., 0.5)
#Vlist = [1.0]

NSKIP = 200
THERMALIZATION = 10**5
SWEEPS = 10**6
Nscratch = 200
Ntau = 1000

wtime = '1:00:00'
tmin = 300
tmax = 600
ncores = 160
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

