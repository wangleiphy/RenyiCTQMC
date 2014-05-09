#include "interaction_expansion.hpp"

//wanglandau before actual simuation 
//here we also combine the simulation of two ensembles (in principle we can seperate them)
void  InteractionExpansion::wanglandau() 
{

 for (unsigned iter =0; iter< 10; ++iter) {

    //reset histogram 
    for (unsigned i=0; i<2; ++i) {
    }
 
    sweeps =0; 
    //thermolization 
    do{
      sweeps++;
      interaction_expansion_step();                
 
      if(sweeps % recalc_period ==0)
           reset_perturbation_series();
     }while (sweeps  <  therm_steps) ; 
 
    //accumute hisogram 
    do{
      sweeps++;
      interaction_expansion_step();                
 
      if(sweeps % recalc_period ==0)
           reset_perturbation_series();
 
      if(sweeps % measurement_period ==0)
         {//measure 
          unsigned pert_order = Msuper.num_vertices(); 
          pertorder_hist[sector][pert_order]++; 
          lng[sector][pert_order] += lnf[sector]; 
         }
     }while (sweeps  < therm_steps + mc_steps*measurement_period) ; 

    //if it is not flat enough, refine it  
    for (unsigned i=0; i<2; ++i) 
        if (not pertorder_hist[i].is_flat(pertorder_hist[i].top_index())){
        
            lnf[i] /= 2.; 
            pertorder_hist[i].clear(); 
        }

 }
}
