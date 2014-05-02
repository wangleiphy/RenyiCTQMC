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
  super_green_function(unsigned ntime, const Mat& KAB, const Mat& KABprime, const itime_t beta, const unsigned NA, const unsigned NB,  const itime_t timestep)
  :nt_(ntime)
  ,ns_(KAB.rows())
  ,tau_(ntime)
  ,gf_(boost::extents[ntime][ntime])
  ,beta_(beta)
  ,NA_(NA)
  ,NB_(NB)
  ,KAB_(KAB)
  ,KABprime_(KABprime)
  {
   
   //std::cout << "K:\n" << K << std::endl; 
    
   /*
   //calculate bare_green function in imaginary time 
   for(itime_index_t it1=0; it1<ntime; ++it1){
      double tau1 =  (double)it1*timestep; 
      tau_[it1] = tau1; 
      for(itime_index_t it2=0; it2<ntime; ++it2){
        double tau2 =  (double)it1*timestep; 
        if (tau1 > tau2)   
            gf_[it1][it2] = B(tau1, tau2)*G(tau2); 
        else
            gf_[it1][it2] = (G(tau1)- Eigen::MatrixXd::Identity(ns_, ns_))*  Binv(tau2, tau1); 
      }
   }
   */

   std::cout << "super_green_function done" << std::endl; 

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
  ~super_green_function(){
  }

  inline double tau(const itime_index_t it) const {return tau_[it];}

  inline double operator()(const unsigned it1, const unsigned it2, const unsigned s1, const unsigned s2)const{
      return gf_[it1][it2](s1, s2);}


  //direct compute 
  inline double gf(const double tau1, const double tau2, const unsigned site1, const unsigned site2)const{

        site_t s1 = site1; 
        site_t s2 = site2; 
        Mat tmp; 
       
        if (tau1>=beta_ && site1 >= NA_)
             s1 = site1 + NB_; 
       
        if (tau2>=beta_ && site2 >= NA_)
             s2 = site2 + NB_;  

        if (tau1 > tau2)   
            tmp = B(tau1, tau2)*G(tau2); 
        else
            tmp = (G(tau1)- Eigen::MatrixXd::Identity(ns_, ns_))*  Binv(tau2, tau1); 

        return tmp(s1, s2); 
  }


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

  const double beta_; 
  const unsigned NA_, NB_; 
  const Mat& KAB_; 
  const Mat& KABprime_; 

  //helper functions 
  Mat B(const double tau1, const double tau2) const { // B(tau1) ... Btau(tau2)
      assert(tau1>=tau2); 

      if (tau1>=beta_){
        if (tau2>=beta_) 
            return expm(-(tau1-tau2)*KABprime_) ;  
        else
            return expm(-(tau1-beta_)*KABprime_) * expm(-(beta_-tau2)*KAB_); 
      }else{
        return expm(-(tau1-tau2)*KAB_);
    }
  }

  Mat Binv(const double tau1, const double tau2) const {
      assert(tau1>=tau2); 

      if (tau1>=beta_){
        if (tau2>=beta_) 
            return expm((tau1-tau2)*KABprime_) ;  
        else
            return expm((beta_-tau2)*KAB_)*expm((tau1-beta_)*KABprime_); 
      }else{
        return expm((tau1-tau2)*KAB_);
    }
   }

  Mat G(const double tau) const {
    Mat res = Eigen::MatrixXd::Identity(ns_, ns_) + B(tau, 0.) * B(2.*beta_, tau); 
    return res.inverse(); 
  }

  Mat expm(const Mat& A) const {return A.exp();}

};

#endif
