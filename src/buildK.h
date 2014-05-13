#ifndef BUILDK_H
#define BUILDK_H 

#include <boost/foreach.hpp>
#include <alps/lattice.h>
#include "types.h"


Eigen::MatrixXd buildK(const alps::graph_helper<>& lattice){

    //construct the hamiltonian K 
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(lattice.num_sites(), lattice.num_sites()); 
    BOOST_FOREACH(const alps::graph_helper<>::site_descriptor& s1, lattice.sites()) {
        BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& s2, lattice.neighbors(s1)) {
            K(s1, s2) = -1.0; 
        }
    }

   return K; 
}


Eigen::MatrixXd buildKAB(const alps::graph_helper<>& lattice, const unsigned NA){
    
    unsigned NB = lattice.num_sites() - NA; 

    //construct the hamiltonian 
    Eigen::MatrixXd KAB = Eigen::MatrixXd::Zero(lattice.num_sites()+NB, lattice.num_sites()+NB); 
    BOOST_FOREACH(const alps::graph_helper<>::site_descriptor& s1, lattice.sites()) {
        BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& s2, lattice.neighbors(s1)) {
            KAB(s1, s2) = -1.0; 
        }
    }

   return KAB; 
}


Eigen::MatrixXd buildKABprime(const alps::graph_helper<>& lattice, const unsigned NA){
    
    unsigned NB = lattice.num_sites() - NA; 

    //construct the hamiltonian 
    Eigen::MatrixXd KABprime = Eigen::MatrixXd::Zero(lattice.num_sites()+NB, lattice.num_sites()+NB); 
    BOOST_FOREACH(const alps::graph_helper<>::site_descriptor& s1, lattice.sites()) {
        BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& s2, lattice.neighbors(s1)) {
            KABprime(s1, s2) = -1.0; 
        }
    }

    //swap B and Bprime  
    for (unsigned i=NA; i<lattice.num_sites(); ++i) {
        KABprime.row(i).swap(KABprime.row(i+NB));
        KABprime.col(i).swap(KABprime.col(i+NB));
    }

   return KABprime; 

}


#endif 
