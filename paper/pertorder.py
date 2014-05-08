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
from numpy import array , linspace , sqrt , arange , log , loadtxt , zeros_like 
from config import * 
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

#filter resultFilies
#for f in list(resultFiles):
#    Temp= float(re.search('T([0-9]*\.?[0-9]*)Ntau',f).group(1)) 

data = []
print resultFiles 

data = pyalps.loadMeasurements(resultFiles, ['PertOrder_0', 'PertOrder_1'])
data = pyalps.flatten(data)
print data 

res0 = pyalps.collectXY(data, x='V', y='PertOrder_0',  foreach = ['TEMPERATURE'])
res1 = pyalps.collectXY(data, x='V', y='PertOrder_1',  foreach = ['TEMPERATURE'])

pyalps.propsort(res0,'TEMPERATURE')
res0 = res0[::-1]

pyalps.propsort(res1,'TEMPERATURE')
res1 = res1[::-1]

icolor = 0
for d0, d1 in zip(res0, res1):
    Temp = float(d0.props['TEMPERATURE'])
    d0.props['label'] = r'$\langle k\rangle_{\mathcal{Z}^2},T=%g$'%(Temp)
    d1.props['label'] = r'$\langle k\rangle_{\mathcal{Z}^A},T=%g$'%(Temp)


    d1.props['ylabel'] = 'Perturbation Orders'
    d1.props['xlabel'] = '$V/t$'

    d0.props['color'] = colors[icolor]
    d0.props['line'] = 'o'

    d1.props['color'] = colors[icolor]
    d1.props['line'] = 's'

    icolor += 1 

pyalps.plot.plot(res0)
pyalps.plot.plot(res1)

print pyalps.plot.convertToText(res0)
print pyalps.plot.convertToText(res1)

plt.legend(loc='upper right')
#plt.subplots_adjust(left=0.15)


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
