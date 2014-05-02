#include "interaction_expansion.hpp"

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
void InteractionExpansion::interaction_expansion_step()
{
    double update_type=random();

    if(update_type < probs[0]){     
        std::cout << "before add" << std::endl; 
        add();
        std::cout << "after add" << std::endl; 
    }else if(update_type < probs[1]){
        std::cout << "before remove" << std::endl; 
        remove(); 
        std::cout << "after remove" << std::endl; 
    }else if(update_type < probs[2]){
        std::cout << "before Z2W" << std::endl; 
        Z_to_W(); 
        std::cout << "after Z2W" << std::endl; 
    }else if(update_type < probs[3]){
        std::cout << "before W2Z" << std::endl; 
        W_to_Z(); 
        std::cout << "after W2Z" << std::endl; 
    }
}


void InteractionExpansion::build_matrix(){

    //rebuild super Matrix 
    assert(Msuper.creators().size() == 2*Msuper.num_vertices()); 
    assert(Msuper.creators().size() == M[0].creators().size() + M[1].creators().size()); 

    Msuper.matrix() = Eigen::MatrixXd::Zero(Msuper.creators().size(), Msuper.creators().size());  
    for (unsigned int i=0; i< Msuper.creators().size(); ++i){
        for (unsigned int j=i+1; j< Msuper.creators().size(); ++j){ //do not fill diagonal 
            Msuper.matrix()(i,j) = super_bare_green_itime.gf(Msuper.creators()[i].t(), Msuper.creators()[j].t(), Msuper.creators()[i].s(), Msuper.creators()[j].s());  //super_green0_spline(Msuper.creators()[i], Msuper.creators()[j]); 
            //Msuper.matrix()(j,i) = -Msuper.creators()[i].parity()*Msuper.creators()[j].parity()*Msuper.matrix()(i,j);//anti-symmetrization 
            //Msuper.matrix()(j,i) = super_green0_spline(Msuper.creators()[j], Msuper.creators()[i]); 
            Msuper.matrix()(j,i) = super_bare_green_itime.gf(Msuper.creators()[j].t(), Msuper.creators()[i].t(), Msuper.creators()[j].s(), Msuper.creators()[i].s());  

        }
    }

   std::cout << "Msuperinv from scratch:\n" << Msuper.matrix() << std::endl; 
   Msuper.matrix() = Msuper.matrix().inverse().eval();     

   //std::cout << "Msuper from scratch:\n" << Msuper.matrix() << std::endl; 
   std::cout << "det(Msuper)= " << Msuper.matrix().determinant() << std::endl; 

    
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

    std::cout << "M" << icopy << " from scratch:\n" << M[icopy].matrix() << std::endl; 
    std::cout << "det(M"<< icopy <<") = " << M[icopy].matrix().determinant() << std::endl; 

  }

}

///Every now and then we have to recreate M from scratch to avoid roundoff error.
void InteractionExpansion::reset_perturbation_series()
{
  sign=1.;
  
  if (Msuper.creators().size()<1) return; //do not rebuilt for empty matrix 

  m_matrix::matrix_t Mdiff(Msuper.matrix()); //make a copy of M.matrix()
  build_matrix(); // rebuild Msuper and M 
    
  //reset the weight ratio 
  double new_weight = Msuper.matrix().determinant()/M[0].matrix().determinant()/M[1].matrix().determinant();  

  if (fabs(new_weight/weight-1.)>1E-8) {
    std::cout<<"WARNING: roundoff errors in weight " << weight << " " << new_weight << std::endl;
  }

  weight = new_weight; 

  Mdiff -= Msuper.matrix(); //subtract the new one 
  Mdiff = Mdiff.cwiseAbs(); //and take absolute value 
  double max_diff = Mdiff.maxCoeff(); 


  if(max_diff > 1.e-8){
    std::cout<<"WARNING: roundoff errors " <<max_diff << std::endl;

    std::cout << Mdiff << std::endl; 
    std::cout << "creators: ";  
    for (unsigned  i=0; i< Msuper.creators().size(); ++i) {
      std::cout << Msuper.creators()[i].s()<< "("<< Msuper.creators()[i].t() << ")"  << ","; 
    }
    std::cout << std::endl; 
    abort(); 
  }
}

