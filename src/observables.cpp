#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("PertOrder")
               << alps::ngs::RealObservable("ZtoW")
               << alps::ngs::RealObservable("WtoZ")
               ; 

 measurements  << alps::ngs::RealObservable("Z")
               << alps::ngs::RealObservable("W")
               //<< alps::ngs::RealObservable("IntE")
               ; 

 for (unsigned int i=0; i<2; ++i){
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
 }

}

//this function is called whenever measurements should be performed.
void InteractionExpansion::measure_observables() 
{
  measurements["Sign"]<<sign;
  measurements["PertOrder"] << double(Msuper.num_vertices()); // the pert order in total 

  if (sector==0){
    measurements["Z"] << 1.;
    measurements["W"] << 0.;
    //measurements["IntE"] << double(M.num_vertices());// the pert order measured in Z space 

  }else{
    measurements["Z"] << 0.;
    measurements["W"] << 1.;
    //measurements["IntE"] << 0.; 
  }
}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){
     //results["IntE"] = (-1./beta)*results["IntE"]/results["Z"];

     results.insert("Zratio", eta*results["W"]/results["Z"]);
}
