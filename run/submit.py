#!/usr/bin/env python
from numpy import linspace 
from params import params 
import os.path 
import sys 

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-run", action='store_true', help="Run or not")
    parser.add_argument("-waitfor", type=int, help="wait for this job for finish")
    input = parser.parse_args()
    #default parameters 

    jobdir='../jobs/'
    
    textoutput = 1 # 

    Maxorder = 1024

    Nscratch = 1000 # period to rebuilt M 
    THERMALIZATION = 10**5
    SWEEPS = 10**6
    NSKIP = 100 # skip when doing measurement 
    
    #this import might overwrite the above default parameters 
    #########################################################
    import socket
    machinename = socket.gethostname()
    if 'brutus' in machinename:
        from config.brutus import * 
    elif 'rosa' in machinename:
        from config.rosa import * 
    elif 'monch' in machinename:
        from config.monch import * 
    else:
        print 'where am I ?', machinename 
        sys.exit(1)
    #########################################################

    cmd = ['mkdir', '-p', resfolder]
    subprocess.check_call(cmd)
    

    for L, W in zip(Llist, Wlist):
        for V in Vlist:
            for T in Tlist:

                jobid = input.waitfor 
                for NA0, NA1 in zip(NA0list, NA1list):
                #for NA0 in range(0, L*W, NAstep):
                #        NA1 = NA0 + NAstep 

                        inputfile = params(latticename, L , W, NA0, NA1,   
                                           V=V, T= T, 
                                           Maxorder = Maxorder, Ntau = Ntau, Nscratch = Nscratch, 
                                           SWEEPS=SWEEPS, THERMALIZATION=THERMALIZATION, NSKIP=NSKIP, 
                                           Add = Add, Remove = Remove, 
                                           ZtoW = ZtoW, WtoZ = WtoZ, 
                                           folder=resfolder, textoutput=textoutput)


                        bin = prog + ' ' + inputfile 

                        args = {}
                        jobname = jobdir + os.path.basename(inputfile).replace('.in','')

                        jobid = submitJob(bin,args,jobname,wtime,ncores=ncores,run=input.run, wait=jobid)
