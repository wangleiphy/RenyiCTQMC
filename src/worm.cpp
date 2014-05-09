#include "interaction_expansion.hpp"

//transition between W and Z sector 

void InteractionExpansion::Z_to_W()
{

  if (sector==1) return; 
  int pert_order = Msuper.num_vertices(); 
  double log_metropolis_weight =  log ((WtoZ/ZtoW)*eta) + logweight  + (lng[0][pert_order]-lng[1][pert_order]);

  if(exp(log_metropolis_weight) > random()){
    measurements["ZtoW"]<< 1.;
    sector = 1; 

  }else{

    measurements["ZtoW"]<< 0.;

  }
}

void InteractionExpansion::W_to_Z()
{

  if (sector==0) return; 
  int pert_order = Msuper.num_vertices(); 

  double log_metropolis_weight = log( (ZtoW/WtoZ)/eta) - logweight - (lng[0][pert_order]-lng[1][pert_order]);

  if(exp(log_metropolis_weight) > random()){
    measurements["WtoZ"]<< 1.;
    sector = 0; 

  }else{
    measurements["WtoZ"]<<0.;
  }
}
