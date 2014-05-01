//#define uint unsigned int 
#include<sys/types.h>
/// @file types.h
/// @brief definition for basic types like vectors that are used throughout the solvers.

#ifndef TYPES_H_
#define TYPES_H_
#include <vector>
#include <complex>
#include <map>

/// the array (vector) type to store frequence dependent Green's functions and Weiss fields
//typedef std::vector<double> vector_type;
/// two arrays are needed for a single site problem: one for up and one for down spins
//typedef std::pair<vector_type,vector_type> multiple_vector_type;
/// variant if you want complex vectors
//typedef std::vector<std::complex<double> > complex_vector_type;
/// variant for two complex vectors.
//typedef std::pair<complex_vector_type,complex_vector_type> multiple_complex_vector_type;

///dense matrix
//typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
///complex matrix
//typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> complex_matrix;

///same as complex_vector_type
//typedef std::vector<std::complex<double> > complex_vector;
///same as vector_type
//typedef std::vector<double> double_vector;
///vector of ints
//typedef std::vector<int> int_vector;

///enum for spin up and spin down
//enum  {up=0, down=1} ;
///addressing type for site indices (cluster)
typedef unsigned int site_t;
///addressing type for spin indices
typedef unsigned int spin_t;
///type of imaginary time values
typedef double itime_t;
///addressing type of imaginary time indices (discretized)
typedef unsigned int itime_index_t;
///addressing type of matsubara frequency
//typedef unsigned int frequency_t;

typedef alps::graph_helper<>::site_descriptor site_descriptor; 
typedef std::map<unsigned int, std::vector<site_descriptor> > DistanceMap;

#endif /*TYPES_H_*/
