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
  ,wKAB()
  ,wKABprime()
  ,uKAB()
  ,uKABprime()
  {
   
   //std::cout << "super_green_function done" << std::endl; 

   Eigen::SelfAdjointEigenSolver<Mat> ces;

   ces.compute(KAB_);
   wKAB = ces.eigenvalues();  
   uKAB = ces.eigenvectors(); 

   ces.compute(KABprime_);
   wKABprime = ces.eigenvalues();  
   uKABprime = ces.eigenvectors(); 

   //calculate bare_green function in imaginary time 
   for(itime_index_t it_avg=0; it_avg<ntime; ++it_avg){
      double tau_avg = (double)it_avg*timestep; 
      tau_[it_avg] = tau_avg; 

      for(itime_index_t it_diff=0; it_diff<ntime; ++it_diff){
        double tau_diff =  (double)it_diff*timestep;  // >=0 

        double tau1 = tau_avg + 0.5*tau_diff; 
        double tau2 = tau_avg - 0.5*tau_diff; 

        gf_[it_avg][it_diff] = B(tau1, tau2)*G(tau2); 
      }
   }

   /*
   //output gf for test 
   for(itime_index_t it1=0; it1<ntime; ++it1){
    for(itime_index_t it2=0; it2<ntime; ++it2){
        std::cout << tau_[it1] << " "<< tau_[it2] << " " << gf_[it1][it2](0,5) << std::endl; 
    }
    std::cout << std::endl; 
   }
   abort();
   */
}

  ///destructor
  ~super_green_function(){
  }

  inline double tau(const itime_index_t it) const {return tau_[it];}

  inline double operator()(const unsigned it_avg, const unsigned it_diff, const unsigned s1, const unsigned s2)const{
      return gf_[it_avg][it_diff](s1, s2);}


  //direct compute 
  inline double gf(const double tau1, const double tau2, const unsigned site1, const unsigned site2)const{

        site_t s1 = site1; 
        site_t s2 = site2; 
        Mat res; 
       
        if (tau1>=beta_ && site1 >= NA_)
             s1 = site1 + NB_; 
       
        if (tau2>=beta_ && site2 >= NA_)
             s2 = site2 + NB_;  

        if (tau1 >= tau2)   
            res = B(tau1, tau2)*G(tau2); 
        else
            res = (G(tau1)- Eigen::MatrixXd::Identity(ns_, ns_))*  Binv(tau2, tau1); 

        return res(s1, s2); 
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
