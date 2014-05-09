#include "interaction_expansion.hpp"

void InteractionExpansion::wanglandau_run(const unsigned kc){

  sweeps =0; 
  //accumute hisogram 
  while (sweeps  < wanglandau_steps){
    sweeps++;
    wanglandau_step();                
    
    unsigned pert_order; 
    if (sector==0)
        pert_order = M[0].num_vertices() +  M[1].num_vertices(); 
    else
        pert_order = Msuper.num_vertices(); 

    pertorder_hist[sector][pert_order] += 1.; 

    if (pert_order < kc) //only modify it upto kc 
    {
       lng[sector][pert_order] += lnf[sector]; 
       //substract a common background  such that lng[kc-1] always =0
       if (pert_order ==kc-1){
        for (unsigned i=0; i< kc; ++i)
            lng[sector][i] -= lng[sector][kc-1]; 
       }
    }
 
    if(sweeps % recalc_period ==0)
       wanglandau_rebuildmatrix(); 
   }; 
}

//wanglandau before actual simuation 
//here we also combine the simulation of two ensembles (in principle we can seperate them)
void  InteractionExpansion::wanglandau(const int node) 
{

 //estimate kc using a unmodified run 
 sector = 0; 
 wanglandau_run(0); // kc = 0 means we donot modify the dynamics at all 
 sector = 1; 
 wanglandau_run(0);  
 
 //pick the larger one as the cutoff 
 unsigned kc = std::max ( pertorder_hist[0].top_index(),  pertorder_hist[1].top_index()); 

 for (sector =0; sector < 2; ++sector) {

    if (node ==0){
        std::cout << "kc: " << kc << " initial histogram of sector: " << sector << std::endl;  
        print_histogram(); 
    }
   
    pertorder_hist.clear(); 
   
    //start wang-landau iteration 
    unsigned iter = 0; 
    while (lnf[sector] > 1E-8 && iter < 100){
       wanglandau_run(kc); 
       ++iter; 
   
       //if it is flat enough, refine it  
       if (pertorder_hist[sector].is_flat(kc)){ 
   
           if (node ==0){
               std::cout << "#iteration: " << iter << ", lnf: " << lnf[sector] << std::endl; 
               print_histogram(); 
           }
           pertorder_hist[sector].clear(); 
           lnf[sector] /= 2.; 
       }
    }
 }
}


//after this are add, remove and rebuild matrix code only for wang-landau purpose 
//we have huge code dupulication 
/////////////////////////////////////////////////////
void InteractionExpansion::wanglandau_step()
{
    if(random()< 0.5){     
        wanglandau_add();
    }else{
        wanglandau_remove(); 
    }
}


void InteractionExpansion::wanglandau_add()
{
  // add vertices

  int pert_order = Msuper.num_vertices(); 

  if(pert_order+1 >= max_order) 
    return; 

  itime_t tau = 2.*beta*random(); 

  std::vector<site_t> sites; 
  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));

  double wratio; 
  if (sector==0)
      wratio =  Zadd_impl(tau, sites, true);
  else
      wratio =  Wadd_impl(tau, sites, true);

  double metropolis_weight = (Remove/Add)*(-2.*beta*n_bond*V)/(pert_order+1)* wratio* exp(lng[sector][pert_order]-lng[sector][pert_order+1])    ; // sector =0 (Z) =1 (W)


  if(fabs(metropolis_weight) > random()){

    if (sector==0)
        Zadd_impl(tau, sites, false);
    else
        Wadd_impl(tau, sites, false);

  }else{

  }
}

void InteractionExpansion::wanglandau_remove()
{

    unsigned int pert_order = Msuper.num_vertices(); 

    if(pert_order < 1)
      return;    
    
    unsigned int vertex = randomint(pert_order);// pickup a random vertex 

    double wratio; 
    if (sector==0)
        wratio =  Zremove_impl(vertex, true);
    else
        wratio =  Wremove_impl(vertex, true);

    double metropolis_weight = (Add/Remove)*pert_order/(-2.*beta*n_bond*V)* wratio * exp(lng[sector][pert_order]-lng[sector][pert_order-1]);

    if(fabs(metropolis_weight) > random()){ //do the actual update

      std::stringstream obs_name;
      obs_name<<"VertexRemoval_"<<sector;
      measurements[obs_name.str().c_str()] << 1.;

    if (sector==0)
        Zremove_impl(vertex, false);
    else
        Wremove_impl(vertex, false);

    }else{

    }
}

void InteractionExpansion::wanglandau_rebuildmatrix(){

    if (sector==0){
         for (unsigned icopy=0; icopy<2; ++icopy) {
        
         assert(M[icopy].creators().size() == 2*M[icopy].num_vertices()); 
        
         M[icopy].matrix() = Eigen::MatrixXd::Zero(M[icopy].creators().size(), M[icopy].creators().size());  
         for (unsigned int i=0; i< M[icopy].creators().size(); ++i){
             for (unsigned int j=i+1; j< M[icopy].creators().size(); ++j){ //do not fill diagonal 
                 M[icopy].matrix()(i,j) = green0_spline(M[icopy].creators()[i], M[icopy].creators()[j]); 
                 M[icopy].matrix()(j,i) = -M[icopy].creators()[i].parity()*M[icopy].creators()[j].parity()*M[icopy].matrix()(i,j);//anti-symmetrization 
             }
         }
        
         M[icopy].matrix() = M[icopy].matrix().inverse().eval();     
        }
    }else{
    
         //rebuild super Matrix 
         assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 
         assert(Msuper.creators().size() == M[0].creators().size() + M[1].creators().size()); 

          Msuper.matrix() = Eigen::MatrixXd::Zero(Msuper.creators().size(), Msuper.creators().size());  
          for (unsigned int i=0; i< Msuper.creators().size(); ++i){
              for (unsigned int j=i+1; j< Msuper.creators().size(); ++j){ //do not fill diagonal 

                  Msuper.matrix()(i,j) = super_green0_spline(Msuper.creators()[i], Msuper.creators()[j]); 
                  Msuper.matrix()(j,i) = -Msuper.creators()[i].parity()*Msuper.creators()[j].parity()*Msuper.matrix()(i,j);//anti-symmetrization 

              }
          }

          Msuper.matrix() = Msuper.matrix().inverse().eval();     
         }
}


