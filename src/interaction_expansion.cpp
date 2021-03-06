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
NA(2),  
NB(2),  
K_(buildK(lattice)), 
//KAB_(buildKAB(lattice, NA)), 
//KABprime_(buildKABprime(lattice, NA)), 
//n_site(KAB.rows()),
n_bond(lattice.num_bonds()),
//n_cell(n_site/2),
//distmap(get_distmap(lattice)), 
//disttable(get_disttable(distmap, n_site)), 
//neighbors(get_neighbors(distmap, n_site, boost::lexical_cast<unsigned int>(parms["Nneighbors"]))),
//shellsize(get_shellsize(distmap)), 
//neighborshellsize(get_neighborshellsize(shellsize)), 
//n_tau(boost::lexical_cast<unsigned>(parms["N_TAU"])),
//n_taumeasure(boost::lexical_cast<unsigned int>(parms["N_TAUMEASURE"])),
//n_self(parms["NSELF"] | boost::lexical_cast<unsigned int>(10*n_tau)),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
temperature(boost::lexical_cast<double>(parms["TEMPERATURE"])),                        
beta(1./temperature),  
//timestepinv(n_tau*temperature),  
//timestep(1./timestepinv),// beta/Ntau 
V(boost::lexical_cast<double>(parms["V"])),                        
recalc_period(parms["RECALC_PERIOD"] | 500),
measurement_period(parms["MEASUREMENT_PERIOD"] | 200),
Msuper(2), 
super_bare_green_itime(), 
sweeps(0),
eta(1.), 
logweight(0.),
sign(1.),
Add(boost::lexical_cast<double>(parms["Add"])),
Remove(boost::lexical_cast<double>(parms["Remove"])),
ZtoW(boost::lexical_cast<double>(parms["ZtoW"])),
WtoZ(boost::lexical_cast<double>(parms["WtoZ"])),
probs(),// empty vector 
sector(0), // initialy we are in Z space 
S2(0.),
estimating(true), // initially we do estimation of eta 
parity_(2)
{
   probs.push_back(Add); 
   probs.push_back(Add+Remove); 
   probs.push_back(Add+Remove+ZtoW); 
   probs.push_back(Add+Remove+ZtoW+WtoZ); 

   //build super_green_function 
   NA[0] = boost::lexical_cast<unsigned>(parms["NA0"]); 
   NA[1] = boost::lexical_cast<unsigned>(parms["NA1"]); 

   for (unsigned i=0; i< 2; ++i){

       NB[i] = lattice.num_sites() - NA[i];

       Eigen::MatrixXd KAB = buildKAB(lattice, NA[i]); 
       Eigen::MatrixXd KABprime = buildKABprime(lattice, NA[i]); 
       super_bare_green_itime.push_back(super_green_function(16, KAB, KABprime, beta)); 
    
   //parity for this sector: we copy NB to the extension part  
   for (site_t s = 0; s< lattice.num_sites(); ++s)
       parity_[i].push_back(lattice.parity(s)); 
   for (site_t s=lattice.num_sites(); s<lattice.num_sites()+NB[i]; ++s) 
       parity_[i].push_back(lattice.parity(s-NB[i])); 
   }

   //initialize ALPS observables
   initialize_observables();

   S2 = nonintS2(K_, NA[1], beta) - (NA[0]==0? 0. : nonintS2(K_, NA[0], beta)); 

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
  if (estimating || sweeps  >= therm_steps) 
    measure_observables();
}


double InteractionExpansion::fraction_completed() const {
    return (sweeps < therm_steps ? 0. : ( sweeps - therm_steps )/double(measurement_period)/ double(mc_steps));
}


void InteractionExpansion::estimate_done(double neweta) {
    sweeps = 0; Msuper[0].clear(); Msuper[1].clear(); sector = 0; logweight= 0., sign = 1., eta= neweta;  measurements.reset(false); estimating = false; 

    //std::cout << "eta: " <<  eta << std::endl; 
    //std::cout << "sweeps: " <<  sweeps << std::endl; 
}

void InteractionExpansion::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);

    std::string context = ar.get_context();
    ar.set_context("/simulation/realizations/0/clones/0/checkpoint");
    
    //copy tlist and vlist to vectors 
    std::vector<itime_t> vt; 
    std::vector<site_t> vs; 

    for (unsigned i=0; i< 2*Msuper[0].num_vertices(); ++i){

        vt.push_back(Msuper[0].creators()[i].t()); 
        vs.push_back(Msuper[0].creators()[i].s());  
    }
    
    ar["eta"] << eta;
    ar["sign"] << sign;
    ar["estimating"] << estimating;

    ar["sweeps"] << sweeps;
    ar["sector"] << sector; 
    ar["logweight"] << logweight; 

    ar["vt"] << vt;
    ar["vs"] << vs;

    ar.set_context(context);
}

void InteractionExpansion::load(alps::hdf5::archive & ar) {

    mcbase::load(ar);

    std::vector<itime_t> vt; 
    std::vector<site_t> vs; 

    std::string context = ar.get_context();
    ar.set_context("/simulation/realizations/0/clones/0/checkpoint");

    ar["eta"] >> eta;
    ar["sign"] >> sign;
    ar["estimating"] >> estimating;

    ar["sweeps"] >> sweeps;
    ar["sector"] >> sector; 
    ar["logweight"] >> logweight; 

    ar["vt"] >> vt;
    ar["vs"] >> vs;

    ar.set_context(context);

    //copy vectors to tlist and vlist
    for (unsigned sec = 0; sec< 2; ++sec) {
        for (unsigned i=0; i< vt.size(); ++i){
            Msuper[sec].creators().push_back(creator(vs[i], vt[i])); 
        }
        Msuper[sec].num_vertices() = vt.size()/2; 
    }

    build_matrix(); 
}
