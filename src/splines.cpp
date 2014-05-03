#include "interaction_expansion.hpp"

//as the name says: linearly interpolate between two points
//pass dxinv avoid division 
template<class X, class Y> inline Y linear_interpolate(const X x0, const X dxinv, const Y y0, const Y y1, const X x)
{
  //dxinv = 1/(x1-x0)
  return y0 + (x-x0)*dxinv*(y1-y0);
}

/*
//http://en.wikipedia.org/wiki/Bilinear_interpolation
template<class X, class F> inline F bilinear_interpolate(const X dxinv, const X dyinv, 
                                                         const X x1, const X y1, 
                                                         const X x2, const X y2, 
                                                         const F f11, const F f12, const F f21, const F f22, 
                                                         const X x, const X y
                                                         )
{
  return dxinv * dyinv * (f11*(x2-x)*(y2-y) + f21*(x-x1)*(y2-y) + f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1));
}
*/

///Compute the Green's function G0 (the BARE) green's function between two points
double InteractionExpansion::green0_spline(const creator &cdagger, const creator &c) const
{
  //for the M matrix we need the bare temporal Green's function. 
  //For the dressed frequency GF we need the bare frequency GF.
  //here we get the bare temporal GF, interpolated, as needed for the M matrix.
  //we will receive this as an input into our solver later.
  itime_t delta_t = cdagger.t()-c.t();
  site_t site1 = cdagger.s();
  site_t site2 = c.s();

  //std::cout << "t1, t2, s1, s2 " << cdagger.t() << " " << c.t() << " " << site1 << " " << site2 << std::endl; 
  return green0_spline(delta_t, site1, site2);  
}


///Compute the bare green's function for a given site, and imaginary time.
double InteractionExpansion::green0_spline(const itime_t delta_t, const site_t site1, const site_t site2) const
{
  if(delta_t>=0.){
    int time_index = static_cast<int>(std::floor(delta_t*timestepinv));
    //std::cout << "detla_t, it1, it2 " << delta_t << " " << time_index_1 << " " << time_index_2 << std::endl; 
    return linear_interpolate(bare_green_itime.tau(time_index), timestepinv,
                              bare_green_itime(time_index,site1, site2),
                              bare_green_itime(time_index+1,site1, site2),delta_t);
  }else{
    //delta_t could be a small negative number, careful about precision issue 
    int time_index = static_cast<int>(std::floor(delta_t*timestepinv))+n_tau;
    //std::cout << "detla_t, it1, it2 " << delta_t << " " << time_index_1 << " " << time_index_2 << std::endl; 
    return -linear_interpolate(bare_green_itime.tau(time_index), timestepinv,
                               bare_green_itime(time_index,site1,site2),
                               bare_green_itime(time_index+1,site1,site2),delta_t+beta);
  }
}


///Compute the Green's function G0 (the BARE) green's function between two points
double InteractionExpansion::super_green0_spline(const creator &cdagger, const creator &c) const
{

  return super_green0_spline(cdagger.t(), c.t(), cdagger.s(), c.s());  
}


///Compute the bare green's function for a given site, and imaginary time.
//this is actually exact propagation but not bilinear interpolation 
double InteractionExpansion::super_green0_spline(const itime_t tau1, const itime_t tau2, const site_t site1, const site_t site2) const
{  

  site_t s1 = site1; 
  site_t s2 = site2; 
  
  if (tau1>=beta && site1 >= NA)
       s1 = site1 + NB; 
  
  if (tau2>=beta && site2 >= NA)
       s2 = site2 + NB;  
    
  double res; // we always propagate from small to larger 
  if (tau1>=tau2)
    res = super_bare_green_itime(tau1, tau2, s1, s2); 
  else
    res = -super_bare_green_itime(tau2, tau1, s2, s1)*lattice.parity(site1)*lattice.parity(site2); 

  //compare with direct calculation 
  //double diff =  super_bare_green_itime.fromscratch(tau1, tau2)(s1, s2) -res;
  //if (fabs(diff)>1E-6){
  //    std::cout << "gf:" << tau1 << " " << tau2 << " " << site1 << " " << site2 << " " << res << " " << res + diff << std::endl;
  //    abort(); 
  //}

  return res ;
}
