#include "interaction_expansion.hpp"
#include <iostream>

void InteractionExpansion::print(std::ostream &os) const{
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"***                           CTQMC calculation of Renyi-EE for spinless fermions                       ***"<<std::endl;
  os<<"***                                   Lei Wang, ETH Zurich, 2013-2014                                   ***"<<std::endl;
  os<<"***                                         lewang@phys.ethz.ch                                         ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"n_bond: "<< n_bond << ",\tNA: " << NA << std::endl; 
  os<<"noninteracting S2: "<< S2 << std::endl; 
  os<<"n_tau: "<<n_tau << ",\tmax order: "<<max_order << std::endl; 
  os<<"mc steps: "<<mc_steps << ",\ttherm steps: "<<therm_steps << std::endl;
  os<<"recalc period: "<<recalc_period<<",\tmeasurement period: "<< measurement_period << std::endl; 
  os<<"T: "<<temperature<<",\tV: "<<V<<std::endl;
  
 {
  os<<"probs: ";
  std::ostream_iterator<double> out_it (os,", ");
  copy(probs.begin(), probs.end(), out_it );
 }
  os<<"\teta: "<< eta <<std::endl;
}


void InteractionExpansion::update_params(alps::params &parms) const{
    //std::cout << "before: " << parms.size() << std::endl; 

//    parms["Nsite"] = n_site; 
    parms["Nbond"] = n_bond; 
//    parms["Ncell"] = n_cell; 

    //std::cout << "after: " << parms.size() << std::endl; 
}
