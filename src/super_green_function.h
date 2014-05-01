#ifndef SUPER_GREEN_FUNCTION_H
#define SUPER_GREEN_FUNCTION_H

#include "types.h"
#include <Eigen/Dense>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
#include "boost/multi_array.hpp"

class super_green_function{
public:
  typedef Eigen::MatrixXd  Mat;

  ///constructor: how many time slices, how many sites
  green_function(unsigned int ntime, const Mat& KAB, const Mat& KABprime, const itime_t beta, const itime_t timestep)
  :nt_(ntime)
  ,ns_(KAB.rows())
  ,tau_(ntime)
  ,gf_(boost::extents[ntime][ntime])
  {
   
   //std::cout << "K:\n" << K << std::endl; 

   //calculate bare_green function in imaginary time 
   for(itime_index_t it1=0; it1<ntime; ++it1){
      double tau1 =  (double)it1*timestep; 
      tau_[it1] = tau1; 
      for(itime_index_t it2=0; it2<ntime; ++it2){
        double tau2 =  (double)it1*timestep; 
        if (tau1 > tau2)   
            gf_[it1][it2] = B(tau1, tau2)*gf(tau2); 
        else
            gf_[it1][it2] = (gf(tau1)- Eigen::MatrixXd::Indentity(ns_, ns_))*  Binverse(tau2, tau1); 
      }
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

  }

  ///destructor
  ~green_function(){
  }

  inline double tau(const itime_index_t it) const {return tau_[it];}

  inline double operator()(const unsigned it1, const unsigned it2, const unsigned s1, const unsigned s2)const{
      return gf_[it1][it2](s1, s2);}

  //size information
  unsigned int nsite()const{return ns_;}
  ///return # of imaginary time values
  unsigned int ntime()const{return nt_;}

  
private:

  //const values
  const unsigned int nt_; //imag time points
  const unsigned int ns_; //totoal number of sites ABBprime 
  std::vector<double> tau_; // from 0 to 2*beta 
  boost::multi_array<Mat, 2> gf_;  // gf_[tau1][tau2] is the gf matrix 


  Mat B(const double tau1, const double tau2){
      assert(tau1 > tau2); 

      if (tau1>beta){
        if (tau2>beta) 
            return expm(-(tau1-tau2)*KABprime) ;  
        else
            return expm(-(tau1-beta)*KABprime) * expm(-(beta-tau2)*KAB); 
      }else{
        return expm(-(tau1-tau2)*KAB);
    }
  }

  Mat Binv(const double tau1, const double tau2){
      assert(tau1 > tau2); 

      if (tau1>beta){
        if (tau2>beta) 
            return expm((tau1-tau2)*KABprime) ;  
        else
            return expm((beta-tau2)*KAB)*exp((tau1-beta)*KABprime); 
      }else{
        return expm((tau1-tau2)*KAB);
    }
   }

  Mat gf(const double tau){
    Mat res = Eigen::MatrixXd::Identity(ns_, ns_) + B(tau, 0.) * B(2.*beta, tau); 
    return res.inverse(); 
  }

  Mat expm(const Mat& A){return A.exp();}

};

#endif
