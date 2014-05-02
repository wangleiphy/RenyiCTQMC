import subprocess 
from numpy import arange 

measure_unequaltime = 0

Zupdate = 0.1
Wupdate = 0.1

ZtoW2 = 0.1
W2toZ = 0.1
W2toW4 = 0.1
W4toW2 = 0.1

ZtoW4 = 0.02
W4toZ = 0.1

#nickname = 'firsttry'
#nickname = 'secondtry'
#nickname = 'scanTnew'
#nickname = 'directmeasure'
#nickname = 'worm'
#nickname = 'M4worm_new'
#nickname = 'scalelowerT'
#nickname = 'largestep'
#nickname = 'unequaltime'
nickname = 'Kappa'

nmaxlist = [2]

#Nlist = [72]
#Tlist = [3/(6.)]
#expectedM2 = 0.1; expectedM4 = 0.02

#Nlist = [162]
#Tlist = [3/(9.)]
#expectedM2 = 0.04; expectedM4 = 0.005

Nlist = [288]
Tlist = [3/(12.)]
expectedM2 = 0.02; expectedM4 = 0.002

Vlist = arange(0.1, 1.5, 0.1)

THERMALIZATION = 10**5
SWEEPS = 10**6
NSKIP = 100
Nneighbors = 4

tmin = 60
tmax = 300
ncores = 16
wtime = '4:00'

resfolder = '/cluster/work/scr6/lewang/spinlessdata/' + nickname  + '/'
#h, m = [int(i) for i in wtime.split(':')]
#Tlimit = max(3600*h + 60*m - int(tmax*2) , 0)

prog = 'mpirun ../bin/worm -i '+ str(tmin) + ' -a ' + str(tmax) 

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
