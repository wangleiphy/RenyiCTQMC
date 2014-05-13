#include "interaction_expansion.hpp"

///Compute the Green's function G0 (the BARE) green's function between two points
double InteractionExpansion::super_green0_spline(const unsigned sec, const creator &cdagger, const creator &c) const
{

  return super_green0_spline(sec, cdagger.t(), c.t(), cdagger.s(), c.s());  
}


///Compute the bare green's function for a given site, and imaginary time.
//this is actually exact propagation but not bilinear interpolation 
double InteractionExpansion::super_green0_spline(const unsigned sec, const itime_t tau1, const itime_t tau2, const site_t site1, const site_t site2) const
{  

  site_t s1 = site1; 
  site_t s2 = site2; 
  
  if (tau1>=beta && site1 >= NA[sec])
       s1 = site1 + NB[sec]; 
  
  if (tau2>=beta && site2 >= NA[sec])
       s2 = site2 + NB[sec];  
    
  double res; // we always propagate from small to larger 
  if (tau1>=tau2)
    res = super_bare_green_itime[sec](tau1, tau2, s1, s2); 
  else
    res = -super_bare_green_itime[sec](tau2, tau1, s2, s1)*lattice.parity(site1)*lattice.parity(site2); 

  //compare with direct calculation 
  //double diff =  super_bare_green_itime.fromscratch(tau1, tau2)(s1, s2) -res;
  //if (fabs(diff)>1E-6){
  //    std::cout << "gf:" << tau1 << " " << tau2 << " " << site1 << " " << site2 << " " << res << " " << res + diff << std::endl;
  //    abort(); 
  //}

  return res ;
}
