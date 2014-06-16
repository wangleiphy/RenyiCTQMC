#include "interaction_expansion.hpp"

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
void InteractionExpansion::interaction_expansion_step()
{
    double update_type=random();

    if(update_type < probs[0]){     
//        std::cout << "before add" << std::endl; 
      add();
//        std::cout << "after add" << std::endl; 
  }else if(update_type < probs[1]){
//        std::cout << "before remove" << std::endl; 
      remove(); 
//        std::cout << "after remove" << std::endl; 
  }else if(update_type < probs[2]){
//        std::cout << "before Z2W" << std::endl; 
      Z_to_W(); 
//        std::cout << "after Z2W" << std::endl; 
  }else if(update_type < probs[3]){
//        std::cout << "before W2Z" << std::endl; 
      W_to_Z(); 
//        std::cout << "after W2Z" << std::endl; 
    }
}

/*
void InteractionExpansion::build_matrix(){
  //rebuild matrix by adding vertices back one by one 
    
  //copy creators out as clear will destroy it  
  std::vector<creator> creators = Msuper[0].creators(); 

  Msuper[0].clear(); 
  Msuper[1].clear();  

  logweight = 0.;  
  for (unsigned i=0; i< creators.size()/2; ++i){
     double tau = creators[2*i].t(); 

     std::vector<site_t> sites; 
     sites.push_back(creators[2*i].s());
     sites.push_back(creators[2*i+1].s());

     std::vector<double> wratios = add_impl(tau, sites, true);  
     add_impl(tau, sites, false);  

     logweight += log(fabs(wratios[1])) - log(fabs(wratios[0]));  
    }
   
   //std::cout << "Msuper from scratch:\n" << Msuper.matrix() << std::endl; 
   //std::cout << "det(Msuper)= " << Msuper.matrix().determinant() << std::endl; 
    
   //for (unsigned sec = 0; sec< 2; ++sec){
   // std::cout << "M" << sec << " from scratch:\n" << M[sec].matrix() << std::endl; 
   // std::cout << "det(M"<< sec <<") = " << M[sec].matrix().determinant() << std::endl; 
   //}

   //std::cout << "weight from scratch " <<  M[0].matrix().determinant()*M[1].matrix().determinant()/Msuper.matrix().determinant() << std::endl;  
}
*/


void InteractionExpansion::build_matrix(){

    //rebuild  Matrix 
   for (unsigned sec=0; sec<2; ++sec) {

    Msuper[sec].matrix() = Eigen::MatrixXd::Zero(Msuper[sec].creators().size(), Msuper[sec].creators().size());  
    for (unsigned int i=0; i< Msuper[sec].creators().size(); ++i){
        for (unsigned int j=i+1; j< Msuper[sec].creators().size(); ++j){ //do not fill diagonal 
            Msuper[sec].matrix()(i,j) = super_green0_spline(sec, Msuper[sec].creators()[i], Msuper[sec].creators()[j]); 
            Msuper[sec].matrix()(j,i) = -Msuper[sec].creators()[i].parity()*Msuper[sec].creators()[j].parity()*Msuper[sec].matrix()(i,j);//anti-symmetrization 
        }
    }

    Msuper[sec].matrix() = Msuper[sec].matrix().inverse().eval();     

    //std::cout << "M" << sec << " from scratch:\n" << M[sec].matrix() << std::endl; 
    //std::cout << "det(M"<< sec <<") = " << M[sec].matrix().determinant() << std::endl; 
  }
    //std::cout << "weight from scratch " <<  M[0].matrix().determinant()*M[1].matrix().determinant()/Msuper.matrix().determinant() << std::endl;  
}


///Every now and then we have to recreate M from scratch to avoid roundoff error.
void InteractionExpansion::reset_perturbation_series()
{
  sign=1.;
  
  if (Msuper[0].creators().size()<1) return; //do not rebuilt for empty matrix 
  
  std::vector<m_matrix::matrix_t> Mdiff(2); 
  for (unsigned i=0; i<2; ++i)
    Mdiff[i] = Msuper[i].matrix(); //make a copy of M.matrix()

  build_matrix(); // rebuild Msuper matrix 
  double new_logweight = log(fabs(Msuper[0].matrix().determinant())) - log(fabs(Msuper[1].matrix().determinant()));
  
  //check logweight 
  if ( fabs(exp(logweight-new_logweight)-1.) >1E-6)
      std::cout<<"WARNING: roundoff errors in weight " <<  fabs(exp(logweight-new_logweight)-1.) <<  " " << logweight << " " << logweight << std::endl;

  logweight = new_logweight; 

  for (unsigned i=0; i<2; ++i) {
      //check the difference of M matrix 
      Mdiff[i] -= Msuper[i].matrix(); //subtract the new one 
      Mdiff[i] = Mdiff[i].cwiseAbs(); //and take absolute value 
   }

   Eigen::MatrixXd::Index r0, c0, r1, c1;//keep track the location of the max difference
   //double max_diff = Mdiff[0].maxCoeff(&r0, &c0) + Mdiff[1].maxCoeff(&r1, &c1) ; 

   if( Mdiff[0].maxCoeff(&r0, &c0) + Mdiff[1].maxCoeff(&r1, &c1)  > 1.e-6){
     std::cout<<"WARNING: roundoff errors " << Mdiff[0].maxCoeff(&r0, &c0) << " " <<  Mdiff[1].maxCoeff(&r1, &c1)  << std::endl;  //<<  " " << Msuper[0].matrix()(r0, c0) << " " << Msuper[1].matrix()(r1, c1) <<  std::endl;

    //std::cout << Mdiff << std::endl; 
    //std::cout << "creators: ";  
    //for (unsigned  i=0; i< Msuper.creators().size(); ++i) {
    //  std::cout << Msuper.creators()[i].s()<< "("<< Msuper.creators()[i].t() << ")"  << ","; 
    //}
    //std::cout << std::endl; 
    //abort(); 
  }
}

