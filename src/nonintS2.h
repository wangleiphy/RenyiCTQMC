#ifndef NONINTS2_HPP
#define NONINTS2_HPP

#include <Eigen/Dense>

double RenyiEntropy(Eigen::VectorXd eta, const int n){
     //eta are denmat eigen values 
     //1211.0006
     //1204.4731
 
     double xi, m, w; 

     double res = 0.0; 
     for (unsigned i=0; i< eta.size(); ++i){
         if (eta(i)>0. && eta(i)<1.) {
             xi = log((1.-eta(i))/eta(i)); 
 
             m = (1.+exp(-n*xi));
             w = (1.+exp(-xi));
 
             res += (log(m) - n*log(w)); 
         }
     }
 
     return res/(1.-n);
}

double fermif(const double x){
    return x>0 ? 1./(exp(x)+1.)  : exp(-x)/(exp(-x)+1.);  
}

double nonintS2(const Eigen::MatrixXd& K, const unsigned NA, const double beta){
 
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces;
    ces.compute(K);
    
    //finite temperature noninteracting density matrix 
    Eigen::VectorXd v(K.rows()); 
    for(site_t l=0; l<K.rows(); ++l) 
        v(l) = fermif(beta * ces.eigenvalues()(l)); 
    Eigen::MatrixXd rho = ces.eigenvectors() * v.asDiagonal() * ces.eigenvectors().adjoint(); 
    
    //reduced density matrix 
    Eigen::MatrixXd rhoA = rho.block(0, 0, NA, NA); 
    ces.compute(rhoA);

    return RenyiEntropy(ces.eigenvalues(),2); 
}

#endif 
