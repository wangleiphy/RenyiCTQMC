#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include "types.h"
#include <Eigen/Dense>
#include <vector>

class green_function{
public:
  typedef Eigen::MatrixXd  Mat;

  ///constructor: how many time slices, how many sites
  green_function(unsigned int ntime, const Mat& K, const itime_t beta, const itime_t timestep)
  :nt_(ntime)
  ,ns_(K.rows())
  ,tau_(ntime)
  ,gf_(ntime)
  {

   Eigen::SelfAdjointEigenSolver<Mat> ces;
   ces.compute(K);
   
   /* 
   std::cout << "K:\n" << K << std::endl; 
   std::cout << "single-particle eigenvalues:\n" << ces.eigenvalues() << std::endl; 
   double E0 = 0.;
   for (int l=0; l< nsite; ++l)
       if (ces.eigenvalues()(l)<0.)
         E0 += ces.eigenvalues()(l); 
   std::cout << "E0:" << E0 << std::endl; 
   */
 
   //calculate bare_green function in imaginary time 
   Eigen::VectorXd v(ns_); // always double   // this two 
   for(itime_index_t it=0; it<ntime; ++it){
      tau_[it] = (double)it*timestep; 
      for(site_t l=0; l<ns_; ++l) 
          v(l) = gf(ces.eigenvalues()(l), beta, tau_[it]); 
      gf_[it] = ces.eigenvectors() * v.asDiagonal() * ces.eigenvectors().adjoint(); 
   }

   /*output gf for test 
   for(itime_index_t it=0; it<ntime; ++it){
       std::cout << it << " " << tau_[it] << " " << gf_[it](0,0) 
                                          << " " << gf_[it](0,1) 
                                          << " " << gf_[it](1,0) 
                                          << " " << gf_[it](0,2) 
                                          << " " << gf_[it](2,0) << std::endl; 
   }
   abort();
   */

   //std::cout << "green_function done" << std::endl; 

  }

  ///destructor
  ~green_function(){
  }

  inline double tau(const itime_index_t it) const {return tau_[it];}

  inline double operator()(const unsigned int it, const unsigned int s1, const unsigned int s2)const{return gf_[it](s1, s2);}
  inline const Mat& equaltime()const {return gf_[0];}

  //size information
  unsigned int nsite()const{return ns_;}
  ///return # of imaginary time values
  unsigned int ntime()const{return nt_;}
  
private:

  //const values
  const unsigned int nt_; ///imag time points
  const unsigned int ns_; ///number of sites
  std::vector<double> tau_; 
  std::vector<Mat> gf_; 

  double gf(const double E, const double beta, const double tau){ //<c(t) c^{+}(0)> with t>0 
         return (E>0.) ? exp(-E*tau)/(1.+exp(-beta*E)) :  exp((beta-tau)*E)/(1.+exp(beta*E)) ; // it is samething, to avoid overflow 
  }

};

#endif
