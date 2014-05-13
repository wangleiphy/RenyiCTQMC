#include "interaction_expansion.hpp"

//transition between W and Z sector 

void InteractionExpansion::Z_to_W()
{

  if (sector==1) return; 
  double log_metropolis_weight =  log ((WtoZ/ZtoW)*eta) + logweight;

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
  double log_metropolis_weight = log( (ZtoW/WtoZ)/eta) - logweight;

  if(exp(log_metropolis_weight) > random()){
    measurements["WtoZ"]<< 1.;
    sector = 0; 

  }else{
    measurements["WtoZ"]<<0.;
  }
}
