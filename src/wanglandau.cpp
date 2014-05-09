#include "interaction_expansion.hpp"

void InteractionExpansion::wanglandau_run(unsigned kc){

  sweeps =0; 
  //accumute hisogram 
  do{
    sweeps++;
    interaction_expansion_step();                

    unsigned pert_order = Msuper.num_vertices(); 
    //std::cout << "pert_order " <<  pert_order << std::endl; 
    pertorder_hist[pert_order] += 1.; 

    if (pert_order < kc) //only modify it upto kc 
    {
       lng[pert_order] += lnf; 
       //substract a common background  such that lng[kc-1] always =0
       if (pert_order ==kc-1){
        for (unsigned i=0; i< kc; ++i)
            lng[i] -= lng[kc-1]; 
       }
    }
 
    if(sweeps % recalc_period ==0)
         reset_perturbation_series();

   }while (sweeps  < wanglandau_steps) ; 
}

//wanglandau before actual simuation 
//here we also combine the simulation of two ensembles (in principle we can seperate them)
void  InteractionExpansion::wanglandau() 
{

 //estimate kc 
 wanglandau_run(0); // kc = 0 means we donot modify the dynamics at all 
 unsigned kc = pertorder_hist.top_index(); //unsigned(2.*pertorder_hist.mean()); 
 std::cout << "kc: " << kc << std::endl;  
 print_histogram(); 

 pertorder_hist.clear(); 

 for (unsigned iter =0; iter< 50; ++iter) {
    wanglandau_run(kc); 

    std::cout << "#iteration: " << iter << ", flat: " << std::boolalpha <<  pertorder_hist.is_flat(kc)  << " " << lnf << std::endl; 
    print_histogram(); 

    //if it is flat enough, refine it  
    if (pertorder_hist.is_flat(kc)){ 
        pertorder_hist.clear(); 
        lnf /= 2.; 
    }

 }

}
