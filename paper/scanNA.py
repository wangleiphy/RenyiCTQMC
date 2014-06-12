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
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from config import * 
import re 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-f1", nargs='+', default="params", help="fileheaders")
parser.add_argument("-f2", nargs='+', default="params", help="fileheaders")


parser.add_argument("-y", default="S2", help="observable")
parser.add_argument("-copydata", action='store_true',  help="copy data")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

################################################################
resultFiles1 = []
for fileheader in args.f1:
    resultFiles1 += pyalps.getResultFiles(prefix=fileheader)
resultFiles1 = list(set(resultFiles1))

#filter resultFilies
for f in list(resultFiles1):
    V = float(re.search('V([0-9]*\.?[0-9]*)T',f).group(1)) 
    if V not in [0.5, 1.0, 1.5, 2.0, 2.5]:
        resultFiles1.remove(f)
################################################################
resultFiles2 = []
for fileheader in args.f2:
    resultFiles2 += pyalps.getResultFiles(prefix=fileheader)
resultFiles2 = list(set(resultFiles2))


#filter resultFilies
for f in list(resultFiles2):
    Temp= float(re.search('T([0-9]*\.?[0-9]*)Ntau',f).group(1)) 
    if Temp not in [0.4, 0.6, 0.8, 1.0, 1.2]:
        resultFiles2.remove(f)
################################################################

plt.figure(figsize=(8,8))
ax1 = plt.subplot(211)
#########################
at = AnchoredText("a",prop=dict(size=20), frameon=True,loc=3,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

data = pyalps.loadMeasurements(resultFiles1, args.y)
data = pyalps.flatten(data)
print data 

#first collect using NA1 and perform the cumulative sum 
res = pyalps.collectXY(data, x='NA1', y=args.y,  foreach = ['TEMPERATURE','L','V'])

print res 
pyalps.propsort(res,'V')

icolor = 0
for r in res:
    r.y = cumsum(r.y)

    Temp = float(r.props['TEMPERATURE'])
    V = float(r.props['V'])
    r.props['label'] = '$V=%.1f$'%(V)

    r.props['ylabel'] = '$S_2$'
    r.props['xlabel'] = '$N_A$'

    r.props['color'] = colors[icolor]
    r.props['line'] = '-o'
    icolor += 1 

#print pyalps.plot.convertToText(res)
pyalps.plot.plot(res)
plt.legend(loc='upper left')
plt.setp( ax1.get_xticklabels(), visible=False)
plt.xlabel('')
plt.title('$T=0.5t$')

################################################################
ax2 = plt.subplot(212, sharex=ax1)
#########################
at = AnchoredText("b",prop=dict(size=20), frameon=True,loc=3,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

data = pyalps.loadMeasurements(resultFiles2, args.y)
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
    r.props['label'] = '$T=%.1f$'%(Temp)

    r.props['ylabel'] = '$S_2$'
    r.props['xlabel'] = '$N_A$'

    r.props['color'] = colors[icolor]
    r.props['line'] = '-o'
    icolor += 1 

#print pyalps.plot.convertToText(res)
pyalps.plot.plot(res)
plt.title('$V=2t$')
################################################################

plt.xticks([8, 16, 24, 32, 40, 48, 56, 64], ['8', '16', '24', '32', '40', '48', '56', '64'])
plt.legend(loc='upper left')
#plt.subplots_adjust(hspace=0.1)


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
