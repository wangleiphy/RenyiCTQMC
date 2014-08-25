#include "interaction_expansion.hpp"

//wrapper functions add and remove 

std::vector<double> InteractionExpansion::add_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight) {
    
    std::vector<double> wratios(2); 
    for (unsigned sec=0; sec<2; ++sec)
         wratios[sec] = Wadd_impl(tau, sites, sec, compute_only_weight); 

    return wratios; 
}


std::vector<double> InteractionExpansion::remove_impl(const unsigned vertex, const bool compute_only_weight)
{
    std::vector<double> wratios(2); 
    for (unsigned sec=0; sec<2; ++sec)
        wratios[sec] = Wremove_impl(vertex, sec, compute_only_weight); 

    return wratios; 
}

/*implement add vertex in W sector*/
double InteractionExpansion::Wadd_impl(const double tau, const std::vector<site_t>& sites, const unsigned sec, const bool compute_only_weight)
{
 
    assert(Msuper[sec].matrix().rows() == Msuper[sec].matrix().cols()); 
    assert(Msuper[sec].matrix().rows() == Msuper[sec].creators().size()); 

    unsigned Msize = Msuper[sec].matrix().rows();
 
   if(compute_only_weight){

        Eigen::MatrixX2d Q(Msize, 2);

        Msuper[sec].Stilde = Eigen::Matrix2d::Zero(); // diagonal term is always zero 
        Msuper[sec].R.resize(2, Msize);
        Msuper[sec].MQ.resize(Msize, 2);
       
        Msuper[sec].Stilde(0, 1) = super_green0_spline(sec, tau, tau, sites[0], sites[1]);  
        Msuper[sec].Stilde(1, 0) = -Msuper[sec].Stilde(0, 1)* parity(sec, tau, sites[0])* parity(sec, tau, sites[1]);   
           
        for(unsigned i=0; i<2; ++i){
         for(unsigned j=0; j< Msize; ++j){
              Q(j,i) = super_green0_spline(sec, Msuper[sec].creators()[j].t(), tau, Msuper[sec].creators()[j].s(), sites[i]);
              //Q(j,i) = super_bare_green_itime.gf(Msuper.creators()[j].t(), tau, Msuper.creators()[j].s(), sites[i]);
              Msuper[sec].R(i,j) = -parity(sec, tau, sites[i]) *parity(sec, Msuper[sec].creators()[j].t(), Msuper[sec].creators()[j].s())* Q(j,i);//anti-symmetrization 
          }
        }
        
        if(Msize>0){
          Msuper[sec].MQ.noalias() = Msuper[sec].matrix() * Q; 
          Msuper[sec].Stilde.noalias() -= Msuper[sec].R * Msuper[sec].MQ; 
        }

        
        return  Msuper[sec].Stilde.determinant();// we have not yet perform the inverse, so it is actually 1./det(Stilde)
  }
 
  Msuper[sec].Stilde = Msuper[sec].Stilde.inverse().eval(); 

  Eigen::Matrix2Xd RM(2, Msize);
  if (Msize>0)
      RM.noalias() = Msuper[sec].R*Msuper[sec].matrix(); 

  Msuper[sec].matrix().conservativeResize(Msize+2, Msize+2);//conservativeResize keep the content 

  //perform the actual update  
  if(Msize>0){
    Msuper[sec].matrix().topRightCorner(Msize, 2).noalias()   = -Msuper[sec].MQ * Msuper[sec].Stilde ;       
    Msuper[sec].matrix().topLeftCorner(Msize, Msize).noalias() -=  Msuper[sec].matrix().topRightCorner(Msize, 2) * RM; 
    Msuper[sec].matrix().bottomLeftCorner(2, Msize).noalias() = -Msuper[sec].Stilde * RM; 
  }

  Msuper[sec].matrix().bottomRightCorner(2, 2) = Msuper[sec].Stilde; 

  for(unsigned i=0; i<2; ++i){
     Msuper[sec].creators().push_back(creator(sites[i], tau)); 
  }

  Msuper[sec].num_vertices() += 1; 

  return 1./Msuper[sec].Stilde.determinant();
}


/*implement remove vertex in W sector*/
double InteractionExpansion::Wremove_impl(const unsigned vertex, const unsigned sec, const bool compute_only_weight)
{

  unsigned Msize = Msuper[sec].matrix().rows();

  if(compute_only_weight){
    // the block we want to remove
    Msuper[sec].Stilde = Msuper[sec].matrix().block<2,2>(2*vertex, 2*vertex);
    return Msuper[sec].Stilde.determinant();
  }

  //swap with the back 
  unsigned pos = 2*vertex; 
  unsigned pos_j = Msize-2; 

  assert(pos <= pos_j); 

  Msuper[sec].matrix().row(pos).swap(Msuper[sec].matrix().row(pos_j));
  Msuper[sec].matrix().row(pos+1).swap(Msuper[sec].matrix().row(pos_j+1));

  Msuper[sec].matrix().col(pos).swap(Msuper[sec].matrix().col(pos_j));
  Msuper[sec].matrix().col(pos+1).swap(Msuper[sec].matrix().col(pos_j+1));

  std::swap(Msuper[sec].creators()[pos], Msuper[sec].creators()[pos_j]);
  std::swap(Msuper[sec].creators()[pos+1], Msuper[sec].creators()[pos_j+1]);
 

  //now perform fastupdate of M
  Msize -= 2; 
  Eigen::MatrixX2d Qtilde(Msize, 2); 
  Eigen::Matrix2Xd Rtilde(2, Msize);

  if(Msize>0){
    Qtilde = Msuper[sec].matrix().topRightCorner(Msize, 2) ; //block(i,j,rows,cols)
    Rtilde = Msuper[sec].matrix().bottomLeftCorner(2, Msize) ; 
  }

  //std::cout << "remove:5" << std::endl; 

  Msuper[sec].matrix().conservativeResize(Msize, Msize); 

  if(Msize>0)
    Msuper[sec].matrix().noalias() -= Qtilde * Msuper[sec].Stilde.inverse() * Rtilde; 

  //std::cout << "remove:6" << std::endl; 
  //get rid of operators
  Msuper[sec].creators().pop_back();
  Msuper[sec].creators().pop_back();

  Msuper[sec].num_vertices() -= 1; 

  return Msuper[sec].Stilde.determinant(); 
}
