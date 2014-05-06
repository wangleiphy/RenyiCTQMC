import subprocess 
from numpy import arange 

Add  = 0.15
Remove = 0.15

ZtoW = 0.35
WtoZ = 0.35
eta = 0.5

latticename = 'open chain lattice'
#latticename = 'open honeycomb lattice'
#latticename = 'cylindrical honeycomb lattice'
###############################
nickname = 'logweight'

Llist = [8]
NAlist = [8]

Tlist = [0.1, 0.2, 0.3, 0.4, 0.5]
Vlist = arange(0.1, 1.6, 0.1)
#Vlist = arange(2., 11., 1.)

NSKIP = 100 
THERMALIZATION = 10**5
SWEEPS = 10**7
Nscratch = 500
##############################

tmin = 60
tmax = 300
ncores = 16
wtime = '8:00'

resfolder = '/cluster/work/scr6/lewang/renyidata/' + nickname  + '/'
#h, m = [int(i) for i in wtime.split(':')]
#Tlimit = max(3600*h + 60*m - int(tmax*2) , 0)

prog = 'mpirun ../bin/main -i '+ str(tmin) + ' -a ' + str(tmax) 

def submitJob(bin,args,jobname,wtime,run=False,ncores=None, wait=[]):

        if run:
            cmd = ['bsub']
        else:
            cmd = ['echo', 'bsub']

        if ncores != None:
            cmd += ['-n', str(ncores)]

        cmd += ['-W',str(wtime)] 
        cmd += ['-J', jobname]  # LSF jobname
        cmd += ['-oo', jobname+'.log'] #log file 
        
        if wait is not None:
            if len(wait) > 0:
                conds='ended('+wait[0]+')'
                for w in wait[1:]: conds += '&&ended('+w+')'
                cmd += ['-w',conds]

        #cmd += ['-R', '\"rusage[mem=4096]\"'] #memory usage 

        cmd += [bin] 

        for key, val in args.items():
            cmd += [key, val]

        subprocess.check_call(cmd)
        #time.sleep(2)
        return None 
