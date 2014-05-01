#include "interaction_expansion.hpp"

void InteractionExpansion::Z_to_W()
{

  if (not Zflag) return; 


  double metropolis_weight = (WtoZ/ZtoW)/weight;

  if(fabs(metropolis_weight) > random()){
    measurements["ZtoW"]<< 1.;

    Zflag = false; 
    sign*=metropolis_weight<0.?-1.:1.;

  }else{

    measurements["ZtoW"]<< 0.;

  }
}

void InteractionExpansion::W_to_Z()
{

    if (Zflag) return; 

    double metropolis_weight = (ZtoW/WtoZ)*weight;

    if(fabs(metropolis_weight) > random()){ //do the actual update
      measurements["WtoZ"]<< 1.;

      Zflag = true; 

      sign*=metropolis_weight<0.?-1.:1.;

    }else{
      measurements["WtoZ"]<<0.;

    }
}


void InteractionExpansion::W_to_W()
{

  if (Zflag) return; 
  measurements["WtoW"]<< 1.;

}
