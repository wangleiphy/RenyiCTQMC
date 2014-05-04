#ifndef WEAK_COUPLING_H
#define WEAK_COUPLING_H

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/lattice.h>
#include <alps/alea.h>

#include <cmath>
#include "green_function.h"
#include "super_green_function.h"
#include "types.h"
#include "mmatrix.hpp"
#include "operator.hpp"
#include <boost/chrono.hpp>
#include <boost/tuple/tuple_io.hpp>

/*types*/
class c_or_cdagger;

class InteractionExpansion: public alps::mcbase
{
public:
  typedef boost::chrono::high_resolution_clock clock;

  InteractionExpansion(alps::params& p, int rank);
  ~InteractionExpansion() {}

  void update();
  void measure();
  double fraction_completed() const;
    
  //perterbation orders in super and 0, 1 sectors 
  boost::tuple<unsigned, unsigned, unsigned> pertorders() const 
      {return boost::make_tuple(Msuper.num_vertices(), M[0].num_vertices(), M[1].num_vertices());}
  double get_weight() const {return exp(logweight);}
  unsigned get_sector() const {return sector;}

  //print progress 
  unsigned long progress() const {return sweeps;}

  void build_matrix(); 
  void test(); 

  void evaluate(results_type& results);


private:
  
  /*functions*/
  // in file io.cpp
  void print(std::ostream &os) const; //print parameters 
  void update_params(alps::params &parms) const; 
  
  /*green's function*/
  // in file spines.cpp
  double green0_spline(const creator &cdagger, const creator &c) const;
  double green0_spline(const itime_t delta_t, const site_t site1, const site_t site2) const;

  double super_green0_spline(const creator &cdagger, const creator &c) const;
  double super_green0_spline(const itime_t tau1, const itime_t tau2, const site_t site1, const site_t site2) const;
  
  /*the actual solver functions*/
  // in file solver.cpp
  void interaction_expansion_step();
  void reset_perturbation_series();

  // in file model.cpp 
  // add and remove vertices 
  void add();
  void remove();

  // in file worm.cpp 
  // create and destroy worm 
  void Z_to_W(); 
  void W_to_Z(); 
  void W_to_W(); 

  //add and remove vertex in worm space 
  //void wormadd(); 
  //void wormremove(); 
  //shift worm 
  //void wormshift(); 
  //create/destroy worm by open a vertex 
  //void wormopen(); 
  //void wormclose(); 
 

  // in file update.cpp:
  // add or remove vertex  
  std::vector<double> add_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight); 
  std::vector<double> remove_impl(const unsigned vertex, const bool compute_only_weight);

  //partition funciton sector
  double Zadd_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight); 
  double Zremove_impl(const unsigned vertex, const bool compute_only_weight);
  //worm sector 
  double Wadd_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight); 
  double Wremove_impl(const unsigned vertex, const bool compute_only_weight);

  //in file wormupdate.cpp:
  //create or destroy the worm
  //double wormcreate_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight);
  //double wormdestroy_impl(const bool compute_only_weight, const unsigned int Morder);
  
  //add or remove vertex in the presence of worm
  //double wormadd_impl(const std::vector<double>& taus, const std::vector<site_t>& sites, const bool compute_only_weight);
  //double wormremove_impl(const std::vector<unsigned int>& vertices, const bool compute_only_weight);
  //shift a worm 
  //double wormshift_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight);
  //open/close a vertex to get worm  
  //double wormopen_impl(const unsigned int vertex_nr, const bool compute_only_weight); 
  //double wormclose_impl(const bool compute_only_weight);
  
  /*measurement functions*/
  // in file observables.cpp
  void initialize_observables();
  void measure_observables();

  //in file measure.cpp
  //void measure_local(); 
  //void measure_nncorrelation(); 
  //in file unequaltime.cpp 
  //void measure_gf();     // <c(t) c^{+}(0)>
  //void measure_ntaun();  // <(n(t)-1/2)(n(0)-1/2)>
 
  /*private member variables, constant throughout the simulation*/
  const alps::Parameters Params;
  const alps::graph_helper<> lattice; 
  const unsigned int max_order;                        
  //const spin_t n_flavors;                          //number of flavors 
 
  const unsigned NA;                                 // number of sites in partition A
  const unsigned NB;                                 // 
  Eigen::MatrixXd K_;                                //the kinetic energy matrix of ABB' system 
  Eigen::MatrixXd KAB_;
  Eigen::MatrixXd KABprime_;

  //const site_t n_site;                               //number of sites
  const site_t n_bond;                               //number of *interaction* bond (fine when n.n. hopping and V )
  //const site_t n_cell;                             //number of unit cells = n_site/2 
  
  //graph stuff 
  //std::vector<DistanceMap> distmap;             //  vector<map(dist:vector<sites>)>
  //Eigen::MatrixXi          disttable;           //  table(si, sj) = dist  
  //std::vector<std::vector<site_t> > neighbors;  //neighbor list for each site within Nneighbors hoppings 
  //std::vector<unsigned int> shellsize;          //number of sites in dist steps 
  //std::vector<unsigned int> neighborshellsize;  //number of sites in dist+1 and dist-1 steps

//  const frequency_t n_matsubara;        //number of matsubara freq
//  const frequency_t n_matsubara_measurements;        //number of measured matsubara freq
  const itime_index_t n_tau;                        //number of imag time slices
//  const itime_index_t n_taumeasure;                 // number of imag time where we do measurement 
//  const frequency_t n_self;                        //number of self energy (W) binning points
  const boost::uint64_t mc_steps;                        
  const unsigned long therm_steps;                
  
  const itime_t temperature;                               
  const itime_t beta;  
  const itime_t timestepinv;// n_tau * temperature 
  const itime_t timestep;   // 1/(n_tau * temperature) = beta/ n_tau  
  const double V;                        
  //const double delta; 
  
  const unsigned int recalc_period;                
  const unsigned int measurement_period;        
  //const unsigned int convergence_check_period;        
  
  //M matrix  
  std::vector<m_matrix> M;
  m_matrix Msuper; 

  green_function bare_green_itime;
  super_green_function super_bare_green_itime;
    
  unsigned long sweeps;        

  const double eta; //coef before W
  double logweight;
  double sign;

  const double Add; 
  const double Remove; 
  const double ZtoW;
  const double WtoZ;  //(Zupdate) *2 + Z2W + W2Z + Wupdate = 1
  std::vector<double> probs; 

  unsigned sector; 

  std::vector<std::pair<unsigned, unsigned> > table; // index -> (icopy, vert)

  double S2;

  unsigned int randomint(const unsigned int i) {return random() * i;}//random int [0, i) 

};

#endif
