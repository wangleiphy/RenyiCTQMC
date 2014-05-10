import h5py 
import pyalps 
import argparse
from numpy import array 

import matplotlib.pyplot as plt
import pyalps.plot
import os , sys

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-filename", default="params", help="fileheaders")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-show", action='store_true',  help="show figure right now")
    group.add_argument("-outname", default="result.pdf",  help="output pdf file")

    args = parser.parse_args()


    #load kilst and rhok from file 
    h5file = h5py.File(args.filename, 'r')
    for s in range(2):
        hist = array(h5file['lng'][str(s)]) 
        #hist = array(h5file['pertorder_histogram'][str(s)]) 
        print hist 

        plt.plot(hist, '-o')

    h5file.close()

    plt.xlim([0,20])

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
      
