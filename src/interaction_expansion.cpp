#include "interaction_expansion.hpp"
#include <ctime>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/model/hamiltonian_matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include "buildK.h"
#include "nonintS2.h"

InteractionExpansion::InteractionExpansion(alps::params &parms, int node)
:alps::mcbase(parms,node),
Params(make_deprecated_parameters(parms)), 
lattice(Params),
max_order(boost::lexical_cast<unsigned>(parms["MAX_ORDER"])), // pert_order = [0,  max_order) 
//n_flavors(boost::lexical_cast<unsigned int>(parms["FLAVORS"])),
NA(boost::lexical_cast<unsigned>(parms["NA"])),  
NB(lattice.num_sites()-NA),  
K_(buildK(lattice)), 
KAB_(buildKAB(lattice, NA)), 
KABprime_(buildKABprime(lattice, NA)), 
//n_site(KAB.rows()),
n_bond(lattice.num_bonds()),
//n_cell(n_site/2),
//distmap(get_distmap(lattice)), 
//disttable(get_disttable(distmap, n_site)), 
//neighbors(get_neighbors(distmap, n_site, boost::lexical_cast<unsigned int>(parms["Nneighbors"]))),
//shellsize(get_shellsize(distmap)), 
//neighborshellsize(get_neighborshellsize(shellsize)), 
n_tau(boost::lexical_cast<unsigned>(parms["N_TAU"])),
//n_taumeasure(boost::lexical_cast<unsigned int>(parms["N_TAUMEASURE"])),
//n_self(parms["NSELF"] | boost::lexical_cast<unsigned int>(10*n_tau)),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
wanglandau_steps((unsigned long)parms["WL_STEPS"]),
wanglandau_convg((double)parms["WL_CONVG"]),
temperature(boost::lexical_cast<double>(parms["TEMPERATURE"])),                        
beta(1./temperature),  
timestepinv(n_tau*temperature),  
timestep(1./timestepinv),// beta/Ntau 
V(boost::lexical_cast<double>(parms["V"])),                        
recalc_period(parms["RECALC_PERIOD"] | 500),
measurement_period(parms["MEASUREMENT_PERIOD"] | 200),
M(2),  // there are two copies of M 
Msuper(), 
bare_green_itime(n_tau+1, K_, beta, timestep),
super_bare_green_itime(16, KAB_, KABprime_, beta),
sweeps(0),
eta(boost::lexical_cast<double>(parms["eta"])),
logweight(0.),
sign(1.),
Add(boost::lexical_cast<double>(parms["Add"])),
Remove(boost::lexical_cast<double>(parms["Remove"])),
ZtoW(boost::lexical_cast<double>(parms["ZtoW"])),
WtoZ(boost::lexical_cast<double>(parms["WtoZ"])),
probs(),// empty vector 
sector(0), // initialy we are in Z space 
table(),
S2(nonintS2(K_, NA, beta)), 
pertorder_hist(2),
lng(2),
lnf(2, 1.),
wanglandau_scalingfactor(2)
{
   probs.push_back(Add); 
   probs.push_back(Add+Remove); 
   probs.push_back(Add+Remove+ZtoW); 
   probs.push_back(Add+Remove+ZtoW+WtoZ); 

   for (unsigned i=0; i< 2; ++i){
       pertorder_hist[i].resize(max_order, 0.); 
       lng[i].resize(max_order, 0.); 
       wanglandau_scalingfactor[i].resize(max_order, 1.); 
   }


   //initialize ALPS observables
   initialize_observables();
    
   if(node==0) {
       print(std::cout); // print parameters to screen 
       update_params(parms); //write back a few generated params 
   }
    
   //perform wang-landau to get lng to flat the pertorder histogram  
   if (wanglandau_steps > 0)
       wanglandau(node); 

   reset(); // reset matrix, weight , sweeps ...  
    
   //set the wang-landau scaling factor exp(lng)
   for (unsigned s=0; s< 2; ++s){
     for (unsigned i=0; i<max_order; ++i){
           wanglandau_scalingfactor[s][i]  = exp(lng[s][i]); 
      }
   }

}


void InteractionExpansion::update()
{
  for(unsigned int i=0;i<measurement_period;++i){
    sweeps++;
    interaction_expansion_step();                
    if(sweeps % recalc_period ==0)
      reset_perturbation_series();
  }
}

void InteractionExpansion::measure(){
  if (sweeps  <  therm_steps) 
   {
    //do nothing 
   }else{
    measure_observables();
   } 
}


double InteractionExpansion::fraction_completed() const {
    return (sweeps < therm_steps ? 0. : ( sweeps - therm_steps )/double(measurement_period)/ double(mc_steps));
}

void InteractionExpansion::reset(){

     Msuper.clear(); 
     M[0].clear(); 
     M[1].clear();  
     table.clear(); 
     
     logweight = 0.;  

     sweeps = 0; 
     sector = 0; 
} 


