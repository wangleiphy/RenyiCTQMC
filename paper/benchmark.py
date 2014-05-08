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

#parser.add_argument("-x", default="TEMPERATURE", help="variable")
parser.add_argument("-y", default="dS2", help="observable")
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
for f in list(resultFiles):
    Temp= float(re.search('T([0-9]*\.?[0-9]*)Ntau',f).group(1)) 
    if Temp not in [0.2, 0.4, 0.6, 0.8, 1.0]:
        resultFiles.remove(f)


data = []
print resultFiles 

data = pyalps.loadMeasurements(resultFiles, args.y)
data = pyalps.flatten(data)
print data 

res = pyalps.collectXY(data, x='V', y=args.y,  foreach = ['TEMPERATURE'])
#res = pyalps.collectXY(data, x='L', y=args.y,  foreach = ['TEMPERATURE','V'])

print res 
pyalps.propsort(res,'TEMPERATURE')
res = res[::-1]

icolor = 0
for d in res:
    Temp = float(d.props['TEMPERATURE'])
    d.props['label'] = '$T=%.1f$'%(Temp)
    if args.y=='dS2':
        d.props['ylabel'] = '$dS_2$'
    else:
        d.props['ylabel'] = '$S_2$'

    d.props['xlabel'] = '$V/t$'

    d.props['color'] = colors[icolor]
    d.props['line'] = 'o'

    V, S2, dS2 =  loadtxt('../data/L8LA4T'+str(Temp)+'_exact.dat', unpack=True, usecols = (0,1,2))

    if args.y=='dS2':
        plt.plot(V, dS2, '-', c=colors[icolor])
    else:
        plt.plot(V, S2, '-', c='k')

    icolor += 1 

print pyalps.plot.convertToText(res)
pyalps.plot.plot(res)

plt.legend(loc='upper right')
#plt.subplots_adjust(left=0.15)

plt.ylim([0.2, 2.0])


plt.plot(V, zeros_like(V)+ log(2.), '--', c ='k')
#ugly hack to add one yticks 
ax2 = plt.gca().twinx()
plt.ylim([0.2, 2.0])
plt.yticks([log(2.)])
lab = plt.gca().get_yticks().tolist()
lab[-1] = '$ln2$'
plt.gca().set_yticklabels(lab)
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')


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
