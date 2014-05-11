'''
python load_corr.py -f /cluster/work/scr6/lewang/PQMCDATA/nncorrhofstadter3latticep1periodicL6W6N36U -s
'''
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import sys ,os  
import subprocess 
import socket
import subprocess

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-logscale", action='store_true',  help="logscale")

parser.add_argument("-copydata", action='store_true',  help="copy data")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
resultFiles = list(set(resultFiles))

#filter resultFilies
#for f in list(resultFiles):
#    if ('L3' in f) or ('Nneighbors3' not in f) :
#    if ('Therm1000000S' not in f) :
#    if ('Sweeps10000000Skip1000' not in f):
#    if ('Skip20' not in f):
#    if ('Sweeps1000000' in f):
#     if ('V1.42' in f):
#        resultFiles.remove(f)

print resultFiles 

nncorr = []
for resultFile in resultFiles:

    data = []
    r = 0
    while True:
        try:
            d = pyalps.loadMeasurements([resultFile], 'Weight_%g'%r)[0][0]
            d.props['k'] = r
            d.props['observable'] = 'weight'
            data.append(d)
            r +=1
        except IndexError:
            break

    nncorr.append(pyalps.collectXY(data,'k','weight'))


nncorr = pyalps.flatten(nncorr)


for d in nncorr:
    d.props['xlabel'] = r'$k$'

    V = d.props['V']
    T = d.props['TEMPERATURE']
 
    d.props['label'] =  r'$V=%g,T=%g$'%(V,T)

print nncorr 

if args.copydata:
    for resultFile in resultFiles:
        cmd = ['cp', resultFile, '../data/']
        subprocess.check_call(cmd)

plt.figure()
nncorr = list(nncorr) # we need a list for propsort 
#pyalps.propsort(nncorr,'V') 


if args.logscale:
    plt.gca().set_yscale('log')



pyalps.plot.plot(nncorr)
plt.legend(loc='lower left')


#plt.figure()
#pyalps.plot.plot(farthest)
#plt.legend(loc='upper left')


if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    recipient = "lewang@phys.ethz.ch"
    message = 'Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
    message += '\n' + pyalps.plot.convertToText(nncorr)
    subject = 'Figure: ' + args.outname

    machinename = socket.gethostname()
    if 'brutus' in machinename or 'monch' in machinename:
        pyalps.sendmail(recipient    # email address of recipients 
                       , subject = subject 
                       , message = message 
                       , attachment= args.outname 
                       )
    else:
        cmd = ['sendmail.py', '-t', recipient+',', '-s', 'Automatic email message from ALPS. '+ subject , '-m', message, '-a', args.outname]
        subprocess.check_call(cmd)

    os.system('rm '+args.outname)


