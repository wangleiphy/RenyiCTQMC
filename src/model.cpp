#include "interaction_expansion.hpp"
/// add or remove vertices from either Z or W sector 

void InteractionExpansion::add()
{
  // add vertices

  int pert_order = Msuper.num_vertices(); 

  if(pert_order+1 > max_order) 
    return; 

  itime_t tau = 2.*beta*random(); 

  std::vector<site_t> sites; 
  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));

  // true means compute_only_weight
  std::vector<double> wratios = add_impl(tau, sites, true);    // wratios in the Z and W sectors 
  double metropolis_weight = (Remove/Add)*(-2.*beta*n_bond*V)/(pert_order+1)* wratios[sector]; // sector =0 (Z) =1 (W)

  //if (metropolis_weight<0.){
  //  std::cout << metropolis_weight << " < 0 in add" << std::endl; 
  //  abort();  
  //}

  if(fabs(metropolis_weight) > random()){

    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<sector;
    measurements[obs_name.str().c_str()] << 1.;
    add_impl(tau, sites, false);

    assert(Msuper.creators().size() == Msuper.matrix().rows()); 
    assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 

    sign*=metropolis_weight<0.?-1.:1.;
    
    logweight += log(fabs(wratios[1])) - log(fabs(wratios[0]));  // wratios individually may be negative 

  }else{

    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<sector;
    measurements[obs_name.str().c_str()] << 0.;


    assert(Msuper.creators().size() == Msuper.matrix().rows()); 
    assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 
  }
}


void InteractionExpansion::remove()
{

    unsigned int pert_order = Msuper.num_vertices(); 

    if(pert_order < 1)
      return;    
    
    unsigned int vertex = randomint(pert_order);// pickup a random vertex 

    std::vector<double> wratios = remove_impl(vertex, true); 
    double metropolis_weight = (Add/Remove)*pert_order/(-2.*beta*n_bond*V)* wratios[sector];

    //std::cout << "after remove_impl" << std::endl; 
    //if (metropolis_weight<0.){
    //  std::cout << metropolis_weight << " < 0 in remove" << std::endl; 
    // abort();  
    //}

    if(fabs(metropolis_weight) > random()){ //do the actual update

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<sector;
      measurements[obs_name.str().c_str()] << 1.;

      remove_impl(vertex, false);  // false means really perform, not only compute weight

      assert(Msuper.creators().size() == Msuper.matrix().rows()); 
      assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 

      sign*=metropolis_weight<0.?-1.:1.;

      logweight += log(fabs(wratios[1])) - log(fabs(wratios[0]));  // wratios individually may be negative 

    }else{

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<sector;
      measurements[obs_name.str().c_str()] << 0.;

      //do nothing
      assert(Msuper.creators().size() == Msuper.matrix().rows()); 
      assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 

    }
}

