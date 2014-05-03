#ifndef SUPER_GREEN_FUNCTION_H
#define SUPER_GREEN_FUNCTION_H

#include "types.h"
#include <Eigen/Dense>
#include <vector>
//#include <unsupported/Eigen/MatrixFunctions>
#include "boost/multi_array.hpp"

class super_green_function{
public:
  typedef Eigen::MatrixXd  Mat;

  ///constructor: how many time slices, how many sites
  super_green_function(unsigned ntime, const Mat& KAB, const Mat& KABprime, const itime_t beta)
  :nt_(ntime)
  ,ns_(KAB.rows())
  ,tau_(ntime)
  ,gf_(boost::extents[ntime][ntime])
  ,beta_(beta)
  ,dtau_(2.*beta/ntime)
//  ,NA_(NA)
//  ,NB_(NB)
  ,KAB_(KAB)
  ,KABprime_(KABprime)
  ,wKAB()
  ,wKABprime()
  ,uKAB()
  ,uKABprime()
  {
   

   Eigen::SelfAdjointEigenSolver<Mat> ces;

   ces.compute(KAB_);
   wKAB = ces.eigenvalues();  
   uKAB = ces.eigenvectors(); 

   ces.compute(KABprime_);
   wKABprime = ces.eigenvalues();  
   uKABprime = ces.eigenvectors(); 

   //calculate bare_green function in imaginary time use hirsch's direct inversion method 
   for(itime_index_t it=0; it<nt_; ++it)
      tau_[it] = (double)it*dtau_;
 
   Mat hirsch(Mat::Identity(nt_*ns_, nt_*ns_)); 
   hirsch.block(0, (nt_-1)*ns_, ns_, ns_) = B(tau_[nt_-1]+dtau_, tau_[nt_-1]); 
   for(itime_index_t it=1; it<nt_; ++it)
      hirsch.block(it*ns_, (it-1)*ns_, ns_, ns_) = -B(tau_[it], tau_[it-1]); 

   hirsch = hirsch.inverse().eval(); 

   for(itime_index_t it1=0; it1<nt_; ++it1){
    for(itime_index_t it2=0; it2<nt_; ++it2){
        gf_[it1][it2] = hirsch.block(it1*ns_, it2*ns_, ns_, ns_); 
      }
   }

   //output gf for test 
   for(itime_index_t it1=0; it1<ntime; ++it1){
    for(itime_index_t it2=0; it2<ntime; ++it2){
        std::cout << tau_[it1] << " "<< tau_[it2] << " " << gf_[it1][it2](0, 5) << " " << fromscratch(tau_[it1], tau_[it2], 0, 5)  << std::endl; 
   }
    std::cout << std::endl; 
   }
   abort();

   //std::cout << "super_green_function done" << std::endl; 
}

  ///destructor
  ~super_green_function(){
  }

  inline double tau(const itime_index_t it) const {return tau_[it];}


  //direct compute 
  inline double fromscratch(const double tau1, const double tau2, const unsigned site1, const unsigned site2)const{

        Mat res; 
        if (tau1 >= tau2)   
            res = B(tau1, tau2)*G(tau2); 
        else
            res = (G(tau1)- Eigen::MatrixXd::Identity(ns_, ns_))*  Binv(tau2, tau1); 

        return res(site1, site2); 
  }


  inline double operator()(const double tau1, const double tau2, const unsigned site1, const unsigned site2)const{
    //find the nearest point then propagate  
    int it1 = static_cast<int>(std::floor(tau1/dtau_));
    int it2 = static_cast<int>(std::floor(tau2/dtau_));

    return B(tau1, tau_[it1]).row(site1) * gf_[it1][it2] * Binv(tau2, tau_[it2]).col(site2); 
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
  const double dtau_; 

//  const unsigned NA_, NB_; 
  const Mat& KAB_; 
  const Mat& KABprime_; 

  //eigen value and vectors of KAB and KABprime 
  Eigen::VectorXd wKAB, wKABprime; 
  Mat uKAB, uKABprime; 

  //helper functions 
  Mat B(const double tau1, const double tau2) const { // B(tau1) ... Btau(tau2)
      assert(tau1>=tau2); 
      
      if (tau1>=beta_){
        if (tau2>=beta_) 
            return expmKABprime(tau1-tau2) ;  
        else
            return expmKABprime(tau1-beta_) * expmKAB(beta_-tau2); 
      }else{
        return expmKAB(tau1-tau2);
      }
  }

  Mat Binv(const double tau1, const double tau2) const {
      assert(tau1>=tau2); 

      if (tau1>=beta_){
        if (tau2>=beta_) 
            return expmKABprime(-(tau1-tau2)) ;  
        else
            return expmKAB(-(beta_-tau2))*expmKABprime(-(tau1-beta_)); 
      }else{
        return expmKAB(-(tau1-tau2));
      }
   }

  Mat G(const double tau) const {
    Mat res = Eigen::MatrixXd::Identity(ns_, ns_) + B(tau, 0.) * B(2.*beta_, tau); 
    return res.inverse(); 
  }

  Mat expmKAB(const double tau) const {// exp(-tau * KAB)
      Eigen::VectorXd v(ns_); 
      for(site_t l=0; l<ns_; ++l) 
          v(l) = exp(-tau * wKAB(l)); 
      return uKAB * v.asDiagonal() * uKAB.adjoint(); 
  }

  Mat expmKABprime(const double tau) const {//exp(-tau* KABprime)
      Eigen::VectorXd v(ns_); 
      for(site_t l=0; l<ns_; ++l) 
          v(l) = exp(-tau * wKABprime(l)); 
      return uKABprime * v.asDiagonal() * uKABprime.adjoint(); 
  }

};

#endif
