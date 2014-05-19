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
from numpy import array , linspace , sqrt , arange , log , cumsum 
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


#Afiles = []
#ABfiles = []
#filter resultFilies
#for f in list(resultFiles):
#    L = int(re.search('L([0-9]*)W',f).group(1)) 
#    W = int(re.search('W([0-9]*)NA',f).group(1)) 
#    NA = int(re.search('NA([0-9]*)V',f).group(1)) 

#    if (NA == L*W):  
#        ABfiles.append(f)
#    else:
#        Afiles.append(f)

data = pyalps.loadMeasurements(resultFiles, 'S2')
data = pyalps.flatten(data)

#first collect using NA1 and perform the cumulative sum 
data = pyalps.collectXY(data, x='NA1', y='S2',  foreach = ['TEMPERATURE','L','V'])
for d in data:
    d.y = cumsum(d.y)

#print data 

MI = []
for d in data:
    r = pyalps.DataSet()

    r.props = d.props
    r.props['label'] = '$I_2$'
    r.props['observable'] = 'I2'
    L = int(d.props['L'])

    r.y = array([2.*d.y[len(d.y)/2-1] - d.y[-1] ])/(2.*L)

    MI.append(r)

MI = pyalps.collectXY(MI, x='TEMPERATURE', y='I2',  foreach = ['L','V'])

print MI 
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
