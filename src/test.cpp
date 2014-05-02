#include <alps/ngs.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include "interaction_expansion.hpp"

void InteractionExpansion::test(){


   {//this block adds ONE vertex
      std::vector<site_t> sites;  
      sites.push_back(1); 
      sites.push_back(4); 

      double tau = 1.90143; 

      std::vector<double> wratio = add_impl(tau, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 

      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << Msuper.num_vertices() << std::endl; 
      weight *= wratio[1]/wratio[0]; 
      std::cout << "weight: " << weight << std::endl; 
   }
   
   {//this block adds vertex
      std::vector<site_t> sites;  
      sites.push_back(3); 
      sites.push_back(4); 

      double tau = 0.3; 

      std::vector<double> wratio = add_impl(tau, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << Msuper.num_vertices() << std::endl; 
      weight *= wratio[1]/wratio[0]; 
      std::cout << "weight: " << weight << std::endl; 
 
   }

   {//this block adds vertex
      std::vector<site_t> sites;  

      sites.push_back(12); 
      sites.push_back(13); 

      double tau = 1.824; 

      std::vector<double> wratio = add_impl(tau, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << Msuper.num_vertices() << std::endl; 
      weight *= wratio[1]/wratio[0]; 
      std::cout << "weight: " << weight << std::endl; 
 
   }

   {//this block adds vertex
      std::vector<site_t> sites;  

      sites.push_back(1); 
      sites.push_back(12); 

      double tau = 0.824; 

      std::vector<double> wratio = add_impl(tau, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << Msuper.num_vertices() << std::endl; 
      weight *= wratio[1]/wratio[0]; 
      std::cout << "weight: " << weight << std::endl; 
 
   }


   {//this block removes vertex 
      unsigned vertex = 0; 

      std::vector<double> wratio = remove_impl(vertex, false);
      std::cout << "remove vertex" << std::endl; 
      std::cout << "number of vertices: " << Msuper.num_vertices() << std::endl; 
      weight *= wratio[1]/wratio[0]; 
      std::cout << "weight: " << weight << std::endl; 
 
   }


   {//this block adds vertex
      std::vector<site_t> sites;  

      sites.push_back(12); 
      sites.push_back(3); 

      double tau = 0.4231; 

      std::vector<double> wratio = add_impl(tau, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << Msuper.num_vertices() << std::endl; 
      weight *= wratio[1]/wratio[0]; 
      std::cout << "weight: " << weight << std::endl; 

   }

   std::cerr << "cdag:" << std::endl ;
   for (int i=0; i< Msuper.creators().size(); ++i){
       std::cerr << "(" << Msuper.creators()[i].s() << "," << Msuper.creators()[i].t() << ") ";
   }
   std::cerr << std::endl;
   std::cout << "Msuperinv.matrix() from fast update:\n" << Msuper.matrix().inverse() << std::endl; 
   //std::cout << "Msuper.matrix() from fast update:\n" << Msuper.matrix() << std::endl; 
   std::cout << "det(Msuper) from fast update= " << Msuper.matrix().determinant() << std::endl; 
    
   for (unsigned icopy=0; icopy<2; ++icopy){
    std::cout << "M" << icopy << "from fast update:\n" << M[icopy].matrix() << std::endl; 
    std::cout << "det(M"<< icopy <<")= " << M[icopy].matrix().determinant() << std::endl; 
   }

   std::cout << "###############################################" << std::endl; 

}

int main(int argc, char** argv){
   alps::parseargs options(argc, argv); 
   alps::params params(options.input_file);

   InteractionExpansion sim(params, 0); 
   std::cout << "initialization done" << std::endl; 
   sim.test(); 
   sim.build_matrix(); 
   return 0; 
}

