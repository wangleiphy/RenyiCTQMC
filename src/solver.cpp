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
   }else{                       
//        std::cout << "before W2W" << std::endl; 
        W_to_W(); 
//        std::cout << "after W2W" << std::endl; 
    }
}

void InteractionExpansion::build_matrix(){
    assert(M.creators().size() == 2*M.num_vertices()); 


    M.matrix() = Eigen::MatrixXd::Zero(M.creators().size(), M.creators().size());  
    for (unsigned int i=0; i< M.creators().size(); ++i){
        //M.matrix()(i,i) = lattice.parity(M.creators()[i].s())*delta; //staggered on site potential  
        for (unsigned int j=i+1; j< M.creators().size(); ++j){ //do not fill diagonal 
            //std::cout << "i,j " <<  i << " " << j << std::endl; 
            M.matrix()(i,j) = green0_spline(M.creators()[i], M.creators()[j]); 
            M.matrix()(j,i) = -M.creators()[i].parity()*M.creators()[j].parity()*M.matrix()(i,j);//anti-symmetrization 
        }
    }

   //std::cout << "Minv from scratch:\n" << M.matrix() << std::endl; 
   M.matrix() = M.matrix().inverse().eval();     

   //std::cout << "M from scratch:\n" << M.matrix() << std::endl; 
   //std::cout << "det(M)= " << M.matrix().determinant() << std::endl; 
}

///Every now and then we have to recreate M from scratch to avoid roundoff error.
void InteractionExpansion::reset_perturbation_series()
{
  sign=1.;
  
  if (M.creators().size()<1) return; //do not rebuilt for empty matrix 

  m_matrix::matrix_t Mdiff(M.matrix()); //make a copy of M.matrix()
  build_matrix(); 

  int pert_order = M.num_vertices(); 
  double new_weight = std::pow(-V,pert_order)/M.matrix().determinant()/eta;  
  if (fabs(new_weight/weight-1.)>1E-8) {
    std::cout<<"WARNING: roundoff errors " << weight << " " << new_weight << std::endl;
  }

  weight = new_weight; 


  //std::cout << "weight:" << weight << std::endl; 

  Mdiff -= M.matrix(); //subtract the new one 
  Mdiff = Mdiff.cwiseAbs(); //and take absolute value 
  double max_diff = Mdiff.maxCoeff(); 


  if(max_diff > 1.e-8){
    std::cout<<"WARNING: roundoff errors " <<max_diff << std::endl;

    //std::cout << Mdiff << std::endl; 
    //std::cout << "creators: ";  
    //for (unsigned int  i=0; i< M.creators().size(); ++i) {
    //  std::cout << M.creators()[i].s()<< "("<< M.creators()[i].t() << ")"  << ","; 
    //}
    //std::cout << std::endl; 
    //abort(); 
  }
}

