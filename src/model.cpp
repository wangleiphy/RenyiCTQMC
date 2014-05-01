#include "interaction_expansion.hpp"

double npickk( unsigned int n, unsigned int k );  

void InteractionExpansion::add()
{
  //std::cout << "##add##" << std::endl; 

  if (not Zflag) return; 

  int pert_order = M.num_vertices(); 
  //we want to add n vertices  
  unsigned int n = randomint(n_max)+1;

  if(pert_order+n > max_order) 
    return; 

  std::vector<site_t> sites; 
  std::vector<double> taus; 

  for (unsigned int i=0; i< n; ++i){
    itime_t tau = beta*random(); 
    taus.push_back(tau); 
    taus.push_back(tau); 

    alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
    sites.push_back(lattice.source(b));
    sites.push_back(lattice.target(b));
  }

  // true means compute_only_weight
  double metropolis_weight = pow(-beta*V*n_bond, n)/npickk(M.num_vertices()+n, n)*add_impl(taus, sites, true);

  //if (metropolis_weight<0.){
  //  std::cout << metropolis_weight << " < 0 in add" << std::endl; 
  //  abort();  
  //}

  if(fabs(metropolis_weight) > random()){

    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<n;
    measurements[obs_name.str().c_str()] << 1.;

    add_impl(taus, sites, false);
 
    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 

    sign*=metropolis_weight<0.?-1.:1.;
    
    weight *= metropolis_weight; 

  }else{

    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<n;
    measurements[obs_name.str().c_str()] << 0.;



    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 
  }
}


void InteractionExpansion::remove()
{
    //std::cout << "##remove##" << std::endl; 
    if (not Zflag) return; 

    unsigned int pert_order = M.num_vertices(); 
    //we want to remove n vertices  
    unsigned int n = randomint(n_max)+1;

    if(pert_order < n)
      return;    
    
   //pick up n distinguishable vertices  
    std::vector<unsigned int> vertices;
    while (vertices.size()<n){
        unsigned int vertex_nr=randomint( pert_order);// pickup a random vertex 
        if (std::find(vertices.begin(), vertices.end(), vertex_nr) == vertices.end())
            vertices.push_back(vertex_nr); //accept only when it is different from the exsisting ones 
    }
    std::sort (vertices.begin(), vertices.end());  

    //std::cout << "before remove_impl,n:"<< vertices.size() << " "<< n << std::endl; 
    //std::cout << "vertices:\n";  
    //std::ostream_iterator<unsigned int> out_it (std::cout," ");
    //std::copy(vertices.begin(), vertices.end(), out_it );
    //std::cout << std::endl; 

    double metropolis_weight = npickk(pert_order, n)/pow(-beta*V*n_bond, n) * remove_impl(vertices, true);
    //std::cout << "after remove_impl" << std::endl; 
    //if (metropolis_weight<0.){
    //  std::cout << metropolis_weight << " < 0 in remove" << std::endl; 
    // abort();  
    //}

    if(fabs(metropolis_weight) > random()){ //do the actual update

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<n;
      measurements[obs_name.str().c_str()] << 1.;

      remove_impl(vertices, false);  // false means really perform, not only compute weight

      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()); 

      sign*=metropolis_weight<0.?-1.:1.;

      weight *= metropolis_weight; 

    }else{

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<n;
      measurements[obs_name.str().c_str()] << 0.;

      //do nothing
      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()); 
    }
}

