#include "interaction_expansion.hpp"

//transition between W and Z sector 

void InteractionExpansion::Z_to_W()
{

  if (sector==1) return; 

  double metropolis_weight = (WtoZ/ZtoW)*eta*weight;

  if(fabs(metropolis_weight) > random()){
    measurements["ZtoW"]<< 1.;

    sector = 1; 
    sign*=metropolis_weight<0.?-1.:1.;

  }else{

    measurements["ZtoW"]<< 0.;

  }
}

void InteractionExpansion::W_to_Z()
{

    if (sector==0) return; 

    double metropolis_weight = (ZtoW/WtoZ)/eta/weight;

    if(fabs(metropolis_weight) > random()){ //do the actual update
      measurements["WtoZ"]<< 1.;

      sector = 0; 
      sign*=metropolis_weight<0.?-1.:1.;

    }else{
      measurements["WtoZ"]<<0.;
    }
}
