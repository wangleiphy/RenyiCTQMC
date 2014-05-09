#include "interaction_expansion.hpp"

//void InteractionExpansion::thermalization(){
//}

//wanglandau before actual simuation 
//here we also combine the simulation of two ensembles (in principle we can seperate them)
void  InteractionExpansion::wanglandau() 
{

  //estimate kc 
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
        pertorder_hist[pert_order]++; 
       }
  }while (sweeps  < therm_steps + mc_steps*measurement_period) ; 

 unsigned kc = unsigned(1.2*pertorder_hist.average()); 
 pertorder_hist.clear(); 

 std::cout << "kc: " << kc << std::endl;  

 for (unsigned iter =0; iter< 10; ++iter) {

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
          pertorder_hist[pert_order]++; 
          if (pert_order < kc) //only modify it upto kc 
            lng[pert_order] += lnf; 
         }
     }while (sweeps  < therm_steps + mc_steps*measurement_period) ; 

    std::cout << "#iteration: " << iter << std::endl; 
    print_histogram(); 

    //if it is not flat enough, refine it  
    if (not pertorder_hist.is_flat(kc)){ 
        lnf /= 2.; 
        pertorder_hist.clear(); 
     }
 }

}
