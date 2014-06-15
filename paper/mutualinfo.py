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
from config import * 
from numpy import array , linspace , sqrt , arange , log , cumsum 
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import re 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-f1", nargs='+', default="params", help="fileheaders")
parser.add_argument("-f2", nargs='+', default="params", help="fileheaders")

parser.add_argument("-copydata", action='store_true',  help="copy data")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles1 = []
for fileheader in args.f1:
    resultFiles1 += pyalps.getResultFiles(prefix=fileheader)
resultFiles1 = list(set(resultFiles1))

resultFiles2 = []
for fileheader in args.f2:
    resultFiles2 += pyalps.getResultFiles(prefix=fileheader)
resultFiles2 = list(set(resultFiles2))

#filter resultFilies
#for f in list(resultFiles2):
    #L = int(re.search('L([0-9]*)W',f).group(1)) 
    #W = int(re.search('W([0-9]*)NA',f).group(1)) 
    #NA = int(re.search('NA([0-9]*)V',f).group(1)) 
    #Temp= float(re.search('T([0-9]*\.?[0-9]*)Ntau',f).group(1)) 
#    if L==12 and ("eta0.5" in f):
#        resultFiles2.remove(f)

for f in list(resultFiles2):
    #L = int(re.search('L([0-9]*)W',f).group(1)) 
    #W = int(re.search('W([0-9]*)NA',f).group(1)) 
    #NA = int(re.search('NA([0-9]*)V',f).group(1)) 
    Temp= float(re.search('T([0-9]*\.?[0-9]*)Ntau',f).group(1)) 
    if Temp not in [0.8, 0.9, 1.0, 1.1, 1.2]:
        resultFiles2.remove(f)

#################################################################

plt.figure(figsize=(8, 5))
ax1 = plt.subplot(121)

#########################
at = AnchoredText("a",prop=dict(size=20), frameon=True,loc=1,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

data = pyalps.loadMeasurements(resultFiles1, 'S2')
data = pyalps.flatten(data)

#first collect using NA1 and perform the cumulative sum 
data = pyalps.collectXY(data, x='NA1', y='S2',  foreach = ['TEMPERATURE','L','V'])
for d in data:
    d.y = cumsum(d.y)


MI = []
for d in data:
    r = pyalps.DataSet()

    r.props = d.props
    r.props['observable'] = 'I2'

    L = int(d.props['L'])
    r.y = array([2.*d.y[len(d.y)/2-1] - d.y[-1]])/(2.*L)

    MI.append(r)

MI = pyalps.collectXY(MI, x='V', y='I2',  foreach = ['L','TEMPERATURE'])

pyalps.propsort(MI,'L')
icolor = 0
for d in MI:

    L = int(d.props['L'])
    d.props['label'] = '$L=%g$'%(L)
    d.props['color'] = colors[icolor]
    d.props['line'] = '-o'
    d.props['xlabel'] = '$V/t$'
    d.props['ylabel'] = '$I_2/\ell_A$'
    icolor += 1 

print pyalps.plot.convertToText(MI)
pyalps.plot.plot(MI)

plt.legend(loc='lower left')
plt.xlim([1.0, 1.5])


#################################################################
ax2 = plt.subplot(122, sharey=ax1)

#########################
at = AnchoredText("b",prop=dict(size=20), frameon=True,loc=1,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

data = pyalps.loadMeasurements(resultFiles2, 'S2')
data = pyalps.flatten(data)

#first collect using NA1 and perform the cumulative sum 
data = pyalps.collectXY(data, x='NA1', y='S2',  foreach = ['TEMPERATURE','L','V'])
for d in data:
    d.y = cumsum(d.y)

MI = []
for d in data:
    r = pyalps.DataSet()

    r.props = d.props
    r.props['observable'] = 'I2'

    L = int(d.props['L'])
    r.y = array([2.*d.y[len(d.y)/2-1] - d.y[-1]])/(2.*L)

    MI.append(r)

MI = pyalps.collectXY(MI, x='TEMPERATURE', y='I2',  foreach = ['L','V'])
pyalps.propsort(MI,'L')

icolor = 0
for d in MI:

    L = int(d.props['L'])
    d.props['label'] = '$L=%g$'%(L)
    d.props['color'] = colors[icolor]
    d.props['line'] = '-o'
    d.props['xlabel'] = '$T/t$'
    d.props['ylabel'] = '$I_2/\ell_A$'
    icolor += 1 

print pyalps.plot.convertToText(MI)
pyalps.plot.plot(MI)

plt.legend(loc='upper left')
plt.xlim([0.8, 1.2])
plt.xticks([0.8, 0.9, 1.0, 1.1, 1.2], ['0.8', '0.9', '1.0', '1.1', '1.2'])

#################################################################
plt.setp(ax2.get_yticklabels(), visible=False)
plt.ylabel('')
plt.subplots_adjust(wspace=0.15, bottom = 0.12)


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
    message += '\n' + pyalps.plot.convertToText(MI)
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
