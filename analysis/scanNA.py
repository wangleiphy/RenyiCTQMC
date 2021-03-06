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
from config import * 
import re 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

#parser.add_argument("-x", default="TEMPERATURE", help="variable")
parser.add_argument("-y", default="S2", help="observable")
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
#    L = int(re.search('L([0-9]*)W',f).group(1)) 
#    W = int(re.search('W([0-9]*)N',f).group(1)) 
#    N = int(re.search('N([0-9]*)U',f).group(1)) 
#    V= float(re.search('latticeV([0-9]*\.?[0-9]*)T',f).group(1)) 
#    if (N!= L*W):  
#        resultFiles.remove(f)

#    elif ('L3' in f) or ('L9' in f) or ('L15' in f):
#    elif ('Theta155' in f):
#        resultFiles.remove(f)
#    if V not in [ 1.24,  1.26,  1.28,  1.3 ,  1.32,  1.34,  1.36,  1.38,  1.4 , 1.42]:
#        resultFiles.remove(f)

#filter resultFilies
#for f in list(resultFiles):
    #Temp= float(re.search('T([0-9]*\.?[0-9]*)Ntau',f).group(1)) 
    #if Temp not in [0.4, 0.6, 0.8, 1.0, 1.2]:
    #    resultFiles.remove(f)

#    V = float(re.search('V([0-9]*\.?[0-9]*)T',f).group(1)) 
#    if V not in [0.5, 1.0, 1.5, 2.0, 2.5]:
#        resultFiles.remove(f)

data = []
print resultFiles 

data = pyalps.loadMeasurements(resultFiles, args.y)
data = pyalps.flatten(data)
print data 

#first collect using NA1 and perform the cumulative sum 
res = pyalps.collectXY(data, x='NA1', y=args.y,  foreach = ['TEMPERATURE','L','V'])

print res 
pyalps.propsort(res,'TEMPERATURE')
res = res[::-1]

icolor = 0
for r in res:
    r.y = cumsum(r.y)

    Temp = float(r.props['TEMPERATURE'])
    V = float(r.props['V'])
    #r.props['label'] = '$T=%.2f$'%(Temp)
    r.props['label'] = '$V=%.1f,T=%.1f$'%(V,Temp)


    r.props['ylabel'] = '$S_2$'
    r.props['xlabel'] = '$N_A$'

    r.props['color'] = colors[icolor]
    r.props['line'] = '-o'
    icolor = (icolor+1)%len(colors)

print res 

#res = pyalps.collectXY(res, x='V', y='S2',  foreach = ['TEMPERATURE','L'])
#res = pyalps.collectXY(data, x='L', y=args.y,  foreach = ['TEMPERATURE','V'])
#res = pyalps.collectXY(data, x='TEMPERATURE', y=args.y,  foreach = ['L','V'])

print pyalps.plot.convertToText(res)
pyalps.plot.plot(res)
#plt.xlim([0,0.18])
#plt.ylim([0,0.4])

#plt.xticks([8, 16, 24, 32, 40, 48, 56, 64], ['8', '16', '24', '32', '40', '48', '56', '64'])
plt.legend(loc='upper left')


###################################################
#ugly hack to add one yticks 
#ax2 = plt.gca().twinx()
#plt.yticks([log(2.)])
#lab = plt.gca().get_yticks().tolist()
#lab[-1] = '$ln2$'
#plt.gca().set_yticklabels(lab)
###################################################


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
