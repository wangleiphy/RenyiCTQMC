#include "interaction_expansion.hpp"

void InteractionExpansion::wanglandau_run(unsigned kc){

  sweeps =0; 
  //accumute hisogram 
  while (sweeps  < wanglandau_steps){
    sweeps++;
    interaction_expansion_step();                

    unsigned pert_order = Msuper.num_vertices(); 
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

   }; 
}

//wanglandau before actual simuation 
//here we also combine the simulation of two ensembles (in principle we can seperate them)
void  InteractionExpansion::wanglandau(const int node) 
{

 //estimate kc using a unmodified run 
 wanglandau_run(0); // kc = 0 means we donot modify the dynamics at all 
 unsigned kc = pertorder_hist.top_index(); //unsigned(2.*pertorder_hist.mean()); 

 if (node ==0){
    std::cout << "kc: " << kc << " initial histogram" << std::endl;  
     print_histogram(); 
 }

 pertorder_hist.clear(); 

 //start wang-landau iteration 
 unsigned iter = 0; 
 while (lnf > 1E-8 && iter < 100){
    wanglandau_run(kc); 
    ++iter; 

    //if it is flat enough, refine it  
    if (pertorder_hist.is_flat(kc)){ 

        if (node ==0){
            std::cout << "#iteration: " << iter << ", flat: " << std::boolalpha <<  pertorder_hist.is_flat(kc)  << ", lnf: " << lnf << std::endl; 
            print_histogram(); 
        }
        pertorder_hist.clear(); 
        lnf /= 2.; 
    }
 }
}
