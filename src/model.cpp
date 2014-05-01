#include "interaction_expansion.hpp"
/// add or remove vertices from either Z or W sector 

void InteractionExpansion::add()
{
  // add vertices

  int pert_order = Msuper.num_vertices(); 

  if(pert_order+1 > max_order) 
    return; 

  std::vector<site_t> sites; 

  itime_t tau = 2.*beta*random(); 

  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));

  // true means compute_only_weight
  std::vector<double> wratios = -V*add_impl(tau, sites, true);    // wratios in the Z and W sectors 
  double metropolis_weight = (Remove/Add)*2.*beta*n_bond/(pert_order+1)* wratios[sector]; // sector =0 (Z) =1 (W)

  //if (metropolis_weight<0.){
  //  std::cout << metropolis_weight << " < 0 in add" << std::endl; 
  //  abort();  
  //}

  if(fabs(metropolis_weight) > random()){

    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<sector;
    measurements[obs_name.str().c_str()] << 1.;

    add_impl(taus, sites, false);
 
    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 

    sign*=metropolis_weight<0.?-1.:1.;
    
    weight *= wratios[1]/wratios[0]; 

  }else{

    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<sector;
    measurements[obs_name.str().c_str()] << 0.;


    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 
  }
}


void InteractionExpansion::remove()
{

    unsigned int pert_order = Msuper.num_vertices(); 

    if(pert_order < 1)
      return;    
    
    unsigned int vertex_nr=randomint( pert_order);// pickup a random vertex 

    std::vector<double> wratios = remove_impl(vertices, true)/(-V); 
    double metropolis_weight = (Add/Remove)*pert_order/(beta*n_bond)* wratio[sector];

    //std::cout << "after remove_impl" << std::endl; 
    //if (metropolis_weight<0.){
    //  std::cout << metropolis_weight << " < 0 in remove" << std::endl; 
    // abort();  
    //}

    if(fabs(metropolis_weight) > random()){ //do the actual update

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<sector;
      measurements[obs_name.str().c_str()] << 1.;

      remove_impl(vertices, false);  // false means really perform, not only compute weight

      assert(Msuper.creators().size() == Msuper.matrix().rows()); 
      assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 

      sign*=metropolis_weight<0.?-1.:1.;
      weight *= wratios[1]/wratios[0]; 

    }else{

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<sector;
      measurements[obs_name.str().c_str()] << 0.;

      //do nothing
      assert(Msuper.creators().size() == Msuper.matrix().rows()); 
      assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 

    }
}

