#include <Eigen/Dense>
#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/ngs.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <iostream>

#include "buildK.h"
#include "nonintS2.h"


int main(int argc, char** argv){
   alps::parseargs options(argc, argv); 
   alps::params params(options.input_file);

   double beta = 1./boost::lexical_cast<double>(params["TEMPERATURE"]); 

   for (unsigned L =3; L<=45; L+=3){// scan system size 

       //update system size and rebuild lattice 
       params["L"] = L ; 

       alps::Parameters Params(make_deprecated_parameters(params));
       alps::graph_helper<> lattice(Params); 

       unsigned NA =lattice.num_sites()/2; 

       //hamiltonian 
       Eigen::MatrixXd K(buildK(lattice)); 

       double S2 = nonintS2(K, NA, beta); 
 
       std::cout << L << " " << NA << " " << S2 << std::endl ; 
   }
   return 0; 
}



