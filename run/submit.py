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

    Maxorder = 2048

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
    
    jobid = input.waitfor 

    for L, W, NA in zip(Llist, Wlist, NAlist):

                for V in Vlist:
                    for T in Tlist:
           
                        inputfile = params(latticename, L , W, NA,   
                                           V=V, T= T, 
                                           Maxorder = Maxorder, Ntau = Ntau, Nscratch = Nscratch, 
                                           SWEEPS=SWEEPS, THERMALIZATION=THERMALIZATION, NSKIP=NSKIP, 
                                           Add = Add, Remove = Remove, 
                                           ZtoW = ZtoW, WtoZ = WtoZ, 
                                           eta = eta, WLSteps = WLSteps, 
                                           folder=resfolder, textoutput=textoutput)


                        bin = prog + ' ' + inputfile 

                        args = {}
                        jobname = jobdir + os.path.basename(inputfile).replace('.in','')

                        jobid = submitJob(bin,args,jobname,wtime,ncores=ncores,run=input.run, wait=jobid)
