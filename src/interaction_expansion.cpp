#include "interaction_expansion.hpp"
#include <ctime>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/model/hamiltonian_matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include "buildK.h"
//#include "bgl.hpp"

InteractionExpansion::InteractionExpansion(alps::params &parms, int node)
:alps::mcbase(parms,node),
Params(make_deprecated_parameters(parms)), 
lattice(Params),
max_order(boost::lexical_cast<unsigned>(parms["MAX_ORDER"])),
//n_flavors(boost::lexical_cast<unsigned int>(parms["FLAVORS"])),
NA(boost::lexical_cast<unsigned>(parms["NA"])),  
NB(boost::lexical_cast<unsigned>(parms["NB"])),  
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
super_bare_green_itime(2*n_tau+1, KAB_, KABprime_, beta, timestep),
sweeps(0),
eta(boost::lexical_cast<double>(parms["eta"])),
weight(1.),
sign(1.),
Add(boost::lexical_cast<double>(parms["Add"])),
Remove(boost::lexical_cast<double>(parms["Remove"])),
ZtoW(boost::lexical_cast<double>(parms["ZtoW"])),
WtoZ(boost::lexical_cast<double>(parms["WtoZ"])),
probs(),// empty vector 
sector(0), // initialy we are in Z space 
table()
//n_max(boost::lexical_cast<unsigned int>(parms["n_max"]))
//measure_unequaltime(boost::lexical_cast<bool>(parms["MEASURE_UNEQUALTIME"]) | false)
{
   probs.push_back(Add); 
   probs.push_back(Add+Remove); 
   probs.push_back(Add+Remove+ZtoW); 
   probs.push_back(Add+Remove+ZtoW+WtoZ); 

   //initialize ALPS observables
   initialize_observables();

   if(node==0) {
       print(std::cout); // print parameters to screen 
       update_params(parms); //write back a few generated params 
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
