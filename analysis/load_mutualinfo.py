'''
python load_results.py -f ../data/pscanLhofstadter3latticep1periodicL*Ntau2000 -s 
'''
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import os, sys 
import subprocess 
import socket
import argparse
from numpy import array , linspace , sqrt , arange , log 
import re 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-copydata", action='store_true',  help="copy data")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)

resultFiles = list(set(resultFiles))


Afiles = []
ABfiles = []
#filter resultFilies
for f in list(resultFiles):
    L = int(re.search('L([0-9]*)W',f).group(1)) 
    W = int(re.search('W([0-9]*)NA',f).group(1)) 
    NA = int(re.search('NA([0-9]*)V',f).group(1)) 

    if (NA == L*W):  
        ABfiles.append(f)
    else:
        Afiles.append(f)

Adata = pyalps.loadMeasurements(Afiles, 'S2')
ABdata = pyalps.loadMeasurements(ABfiles, 'S2')

S2 = pyalps.collectXY(ABdata, x='TEMPERATURE', y='S2',  foreach = ['L','V'])
S2A = pyalps.collectXY(Adata, x='TEMPERATURE', y='S2',  foreach = ['L','V'])

print S2 
print S2A 

MI = []
for d0, d1 in zip(S2, S2A):
    d = pyalps.DataSet()

    d.props = d0.props
    d.props['label'] = '$I_2$'

    d.x = d0.x 
    d.y = 2.*d1.y - d0.y 

    MI.append(d)

print pyalps.plot.convertToText(MI)
pyalps.plot.plot(MI)

plt.legend(loc='upper left')


if args.copydata:
    for resultFile in resultFiles:
        cmd = ['cp', resultFile, '../data/']
        subprocess.check_call(cmd)

if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    recipient = "lewang@phys.ethz.ch"
    message = 'Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
    message += '\n' + pyalps.plot.convertToText(res)
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
