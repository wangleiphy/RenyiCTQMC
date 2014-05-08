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

   for (unsigned W =3; W<=45; ++W){// scan system width  

       //update system size and rebuild lattice 
       params["W"] = W ; 
       params["L"] = W ; // also change length 

       alps::Parameters Params(make_deprecated_parameters(params));
       alps::graph_helper<> lattice(Params); 

       //hamiltonian 
       Eigen::MatrixXd K(buildK(lattice)); 

       unsigned NA =lattice.num_sites()/2; 
       double S2A = nonintS2(K, NA, beta); 
    
       NA = lattice.num_sites(); 
       double S2 = nonintS2(K, NA, beta); 
 
       std::cout << W << " " << S2A << " "<< 2.*S2A - S2 << std::endl ; 
   }
   return 0; 
}



