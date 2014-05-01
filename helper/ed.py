import pyalps
from subprocess import call 

L=2
W=2
t0 = 1.0
binname = 'fulldiag'

parms = []
for V0 in [0.5, 1., 2., 4.]:
    mu = V0 * 1.5 
    parms.append({
          'LATTICE_LIBRARY'           : '../input/mylattices.xml',
          'MODEL_LIBRARY'             : '../input/mymodels.xml',
   	      'LATTICE'                   : "honeycomb lattice", 
          'MODEL'                     : "spinless fermions",
          'L'                         : L,
          'W'                         : W,
          't0'                        : t0,  
          'V0'                        : V0,
          'mu'                        : mu, 
	      'CONSERVED_QUANTUMNUMBERS'  : 'N',
          'TRANSLATION_SYMMETRY'      : 'true', 
#         'MEASURE_LOCAL[nloc]'       : 'n', # can not measure local when have translational invarance 
#         'TOTAL_MOMENTUM'            : "0 0",
#          'MEASURE_CORRELATIONS[nncorr]': 'n', 
#	      'N_total'                   : N,
# 	      'INITIAL_SITE'              : 0, 
          'DENSITIES'                 : 0 # not devide by Nsite 
#         'MEASURE_CORRELATIONS[corr]': 'cdag:c',
#         'INITIAL_SITE'              : 0
          #"PRINT_EIGENVECTORS"        : 1
        })


#input_file = pyalps.writeInputFiles('honeycomb',parms)
#res = pyalps.runApplication('fulldiag',input_file,writexml=False)


folder = '../data/'+ binname + '/'
for p in parms:
    #parmname = folder+p['LATTICE'].replace(" ", "")+'L'+str(p['L'])+'_W'+str(p['W'])+'_N'+str(p['N_total'])+'_V'+str(p['V0']) 
    parmname = folder+p['LATTICE'].replace(" ", "")+'L'+str(p['L'])+'_W'+str(p['W'])+'_V'+str(p['V0']) 
    input_file = pyalps.writeInputFiles(parmname, [p])
    pyalps.runApplication(binname,input_file ,writexml=False)#,MPI=2)
    #check_call(['bsub', '-n', '32','-oo',input_file.replace('.in.xml','.log'),'-W','60','mpirun', 'sparsediag', '--mpi', input_file])
