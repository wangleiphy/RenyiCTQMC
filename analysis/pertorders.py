import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import os , sys

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheader",default="params", help="fileheader")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = pyalps.getResultFiles(prefix=args.fileheader)

pertorder = pyalps.loadTimeSeries(resultFiles[0], 'PertOrder')
pertorder0 = pyalps.loadTimeSeries(resultFiles[0], 'PertOrder_0')
pertorder1 = pyalps.loadTimeSeries(resultFiles[0], 'PertOrder_1')

pertmin = min(pertorder0.min(), pertorder0.min())
pertmax = max(pertorder1.max(), pertorder1.max())

plt.figure()
n, bins, patches = plt.hist(pertorder, range=(pertmin, pertmax), bins =20)
n, bins, patches = plt.hist(pertorder0, range=(pertmin, pertmax), bins =20)
n, bins, patches = plt.hist(pertorder1, range=(pertmin, pertmax), bins =20)


if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    pyalps.sendmail('lewang@phys.ethz.ch'    # email address of recipients 
            , message='Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
                   , attachment= args.outname 
                   , subject='Figure: ' + args.outname 
                   )
    os.system('rm '+args.outname)
