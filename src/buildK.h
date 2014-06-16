#ifndef BUILDK_H
#define BUILDK_H 

#include <boost/foreach.hpp>
#include <alps/lattice.h>
#include "types.h"


Eigen::MatrixXd buildK(const alps::graph_helper<>& lattice){

    //construct the hamiltonian K 
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(lattice.num_sites(), lattice.num_sites()); 
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()){
         K(lattice.source(b), lattice.target(b)) = (lattice.bond_type(b)==0) ? -1. : (lattice.bond_type(b)==1) ? 1. : 0. ;
         K(lattice.target(b), lattice.source(b)) = K(lattice.source(b), lattice.target(b));
    }

   return K; 
}


Eigen::MatrixXd buildKAB(const alps::graph_helper<>& lattice, const unsigned NA){
    
    unsigned NB = lattice.num_sites() - NA; 

    //construct the hamiltonian 
    Eigen::MatrixXd KAB = Eigen::MatrixXd::Zero(lattice.num_sites()+NB, lattice.num_sites()+NB); 
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()){
         KAB(lattice.source(b), lattice.target(b)) = (lattice.bond_type(b)==0) ? -1. : (lattice.bond_type(b)==1) ? 1. : 0. ;
         KAB(lattice.target(b), lattice.source(b)) = KAB(lattice.source(b), lattice.target(b));
    }

   return KAB; 
}


Eigen::MatrixXd buildKABprime(const alps::graph_helper<>& lattice, const unsigned NA){
    
    unsigned NB = lattice.num_sites() - NA; 

    //construct the hamiltonian 
    Eigen::MatrixXd KABprime = Eigen::MatrixXd::Zero(lattice.num_sites()+NB, lattice.num_sites()+NB); 
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()){
         KABprime(lattice.source(b), lattice.target(b)) = (lattice.bond_type(b)==0) ? -1. : (lattice.bond_type(b)==1) ? 1. : 0. ;
         KABprime(lattice.target(b), lattice.source(b)) = KABprime(lattice.source(b), lattice.target(b));
    }

    //swap B and Bprime  
    for (unsigned i=NA; i<lattice.num_sites(); ++i) {
        KABprime.row(i).swap(KABprime.row(i+NB));
        KABprime.col(i).swap(KABprime.col(i+NB));
    }

   return KABprime; 

}


#endif 
