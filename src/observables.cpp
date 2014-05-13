#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("ZtoW")
               << alps::ngs::RealObservable("WtoZ")
               << alps::ngs::RealObservable("PertOrder")
               ; 

 measurements  << alps::ngs::RealObservable("Z")
               << alps::ngs::RealObservable("W")
               //<< alps::ngs::RealObservable("IntE")
               ; 

 for (unsigned i=0; i<2; ++i){
   {
    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   }

   {
    std::stringstream obs_name;
    obs_name<<"VertexRemoval_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   }

   {
    std::stringstream obs_name;
    obs_name<<"PertOrder_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   }
 }


 /*
 for (unsigned i=0; i< 10; ++i){
    std::stringstream obs_name;
    obs_name<<"Weight_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
 }
 */

}

//this function is called whenever measurements should be performed.
void InteractionExpansion::measure_observables() 
{
    measurements["Sign"]<<sign;
    unsigned pert_order = Msuper[sector].num_vertices(); 

    measurements["PertOrder"] << double(pert_order); 
    
    {
        std::stringstream obs_name;
        obs_name<<"PertOrder_"<<sector;
        measurements[obs_name.str().c_str()] << double(pert_order);  
    }

    /*
    if (pert_order < 10){
        std::stringstream obs_name;
        obs_name<<"Weight_"<<pert_order;
        measurements[obs_name.str().c_str()] << exp(logweight+ lng[0][pert_order]-lng[1][pert_order]);  
    }
    */


  if (sector==0){
    measurements["Z"] << 1.;
    measurements["W"] << 0.;

  }else{
    measurements["Z"] << 0.;
    measurements["W"] << 1.;
  }
}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){
     results.insert("dS2", log(eta*results["Z"]/results["W"])); // dS2 = -log(W/Z)
     results.insert("S2", S2 + results["dS2"]); 
}
