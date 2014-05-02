#include "interaction_expansion.hpp"

//wrapper functions add and remove 

std::vector<double> InteractionExpansion::add_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight) {
    
    std::vector<double> wratios(2); 
    wratios[0] = Zadd_impl(tau, sites, compute_only_weight); 
    wratios[1] = Wadd_impl(tau, sites, compute_only_weight); 

    if (not compute_only_weight) {
        unsigned icopy = tau < beta ? 0 : 1;  
        table.push_back(std::make_pair(icopy, M[icopy].num_vertices())); 
    }

    return wratios; 
}


std::vector<double> InteractionExpansion::remove_impl(const unsigned vertex, const bool compute_only_weight)
{
    std::vector<double> wratios(2); 
    wratios[0] = Zremove_impl(vertex, compute_only_weight); 
    wratios[1] = Wremove_impl(vertex, compute_only_weight); 
    
    if (not compute_only_weight) {
        unsigned icopy = table[vertex].first; 
        unsigned vert = table[vertex].second; 
        
        std::swap(table[vertex], table.back());
        table.pop_back(); 

        //for that icopy decrease all the numbering after vert by one 
        for (std::vector<std::pair<unsigned, unsigned> >::iterator it=table.begin(); it!=table.end(); ++it){
            if (it->first == icopy && it->second>vert)
                --(it->second) ; 
        }
    }

    return wratios; 
}

//////////////////////////////////////////////////////////////////////////
/*implement add vertex in Z sector*/
double InteractionExpansion::Zadd_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight)
{
 
    //ratio of Z, it depends on whether tau is < beta or not 
    unsigned int icopy = tau < beta ? 0: 1;  
    unsigned int Msize = M[icopy].matrix().rows();
 
    Eigen::MatrixXd Stilde = Eigen::MatrixXd::Zero(2, 2); // diagonal term is always zero 
    Eigen::MatrixXd RM(2, Msize), R(2, Msize);
    Eigen::MatrixXd Q(Msize, 2), MQ(Msize, 2);
 
    Stilde(0,1) = green0_spline(0., sites[0], sites[1]);
    Stilde(1,0) = -Stilde(0,1)* lattice.parity(sites[0])* lattice.parity(sites[1]);   
 
    //std::cout << "Stilde:\n"<< Stilde << std::endl; 
       
    for(unsigned int i=0; i<2; ++i){
     for(unsigned int j=0; j< Msize; ++j){
          Q(j,i) = green0_spline(M[icopy].creators()[j].t()-tau, M[icopy].creators()[j].s(), sites[i]);
          R(i,j) = -lattice.parity(sites[i]) * M[icopy].creators()[j].parity()* Q(j,i);//anti-symmetrization 
      }
    }
 
    if(Msize>0){
      MQ.noalias() = M[icopy].matrix() * Q; 
      Stilde.noalias() -= R * MQ; 
    }

  //return weight if we have nothing else to do
  if(compute_only_weight){
    return Stilde.determinant();// we have not yet perform the inverse, so it is actually 1./det(Stilde)
  }

 
  Stilde = Stilde.inverse().eval(); 

  if (Msize>0)
      RM.noalias() = R*M[icopy].matrix(); 

  M[icopy].matrix().conservativeResize(Msize+2, Msize+2);//conservativeResize keep the content 

  //perform the actual update  
  if(Msize>0){
    M[icopy].matrix().topRightCorner(Msize, 2).noalias()   = -MQ * Stilde ;       
    M[icopy].matrix().topLeftCorner(Msize, Msize).noalias() -=  M[icopy].matrix().topRightCorner(Msize, 2) * RM; 
    M[icopy].matrix().bottomLeftCorner(2, Msize).noalias() = -Stilde * RM; 
  }

  M[icopy].matrix().bottomRightCorner(2, 2) = Stilde; 

  for(unsigned int i=0; i<2; ++i)
     M[icopy].creators().push_back(creator(sites[i], lattice.parity(sites[i]), tau)); 

  M[icopy].num_vertices() += 1;

  return 1./Stilde.determinant();
}

/*implement add vertex in W sector*/
double InteractionExpansion::Wadd_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight)
{
 
    assert(Msuper.matrix().rows() == Msuper.matrix().cols()); 
    assert(Msuper.matrix().rows() == Msuper.creators().size()); 

    unsigned int Msize = Msuper.matrix().rows();
 
    Eigen::MatrixXd Stilde = Eigen::MatrixXd::Zero(2, 2); // diagonal term is always zero 
    Eigen::MatrixXd RM(2, Msize), R(2, Msize);
    Eigen::MatrixXd Q(Msize, 2), MQ(Msize, 2);
 
    Stilde(0, 1) = super_green0_spline(tau, tau, sites[0], sites[1]);  
    Stilde(1, 0) = -Stilde(0, 1)* lattice.parity(sites[0])* lattice.parity(sites[1]);   
       
    for(unsigned int i=0; i<2; ++i){
     for(unsigned int j=0; j< Msize; ++j){
          Q(j,i) = super_green0_spline(Msuper.creators()[j].t(), tau, Msuper.creators()[j].s(), sites[i]);
          R(i,j) = -lattice.parity(sites[i]) * Msuper.creators()[j].parity()* Q(j,i);//anti-symmetrization 
      }
    }
 
    if(Msize>0){
      MQ.noalias() = Msuper.matrix() * Q; 
      Stilde.noalias() -= R * MQ; 
    }

 
  //return weight if we have nothing else to do
  if(compute_only_weight){
    return  Stilde.determinant();// we have not yet perform the inverse, so it is actually 1./det(Stilde)
  }

 
  Stilde = Stilde.inverse().eval(); 

  if (Msize>0)
      RM.noalias() = R*Msuper.matrix(); 

  Msuper.matrix().conservativeResize(Msize+2, Msize+2);//conservativeResize keep the content 

  //perform the actual update  
  if(Msize>0){
    Msuper.matrix().topRightCorner(Msize, 2).noalias()   = -MQ * Stilde ;       
    Msuper.matrix().topLeftCorner(Msize, Msize).noalias() -=  Msuper.matrix().topRightCorner(Msize, 2) * RM; 
    Msuper.matrix().bottomLeftCorner(2, Msize).noalias() = -Stilde * RM; 
  }

  Msuper.matrix().bottomRightCorner(2, 2) = Stilde; 

  for(unsigned int i=0; i<2; ++i){
     Msuper.creators().push_back(creator(sites[i], lattice.parity(sites[i]), tau)); 
  }

  Msuper.num_vertices() += 1; 

  return 1./Stilde.determinant();
}


/*implement remove vertex in Z sector*/
double InteractionExpansion::Zremove_impl(const unsigned vertex, const bool compute_only_weight)
{

  //copy and vert number in that copy 
  unsigned icopy = table[vertex].first; 
  unsigned vert = table[vertex].second; 

  unsigned Msize = M[icopy].matrix().rows();

  // the block we want to remove
  Eigen::MatrixXd Stilde = M[icopy].matrix().block<2,2>(2*vert, 2*vert);

  if(compute_only_weight){
    return Stilde.determinant();
  }

  //swap with the back 
  unsigned pos = 2*vert; 
  unsigned pos_j = Msize-2; 

  assert(pos <= pos_j); 

  M[icopy].matrix().row(pos).swap(M[icopy].matrix().row(pos_j));
  M[icopy].matrix().row(pos+1).swap(M[icopy].matrix().row(pos_j+1));

  M[icopy].matrix().col(pos).swap(M[icopy].matrix().col(pos_j));
  M[icopy].matrix().col(pos+1).swap(M[icopy].matrix().col(pos_j+1));

  std::swap(M[icopy].creators()[pos], M[icopy].creators()[pos_j]);
  std::swap(M[icopy].creators()[pos+1], M[icopy].creators()[pos_j+1]);
 

  //now perform fastupdate of M
  Msize -= 2; 
  Eigen::MatrixXd Qtilde(Msize, 2), Rtilde(2, Msize);

  if(Msize>0){
    Qtilde = M[icopy].matrix().topRightCorner(Msize, 2) ; //block(i,j,rows,cols)
    Rtilde = M[icopy].matrix().bottomLeftCorner(2, Msize) ; 
  }

  //std::cout << "remove:5" << std::endl; 

  M[icopy].matrix().conservativeResize(Msize, Msize); 

  if(Msize>0)
    M[icopy].matrix().noalias() -= Qtilde * Stilde.inverse() * Rtilde; 

  //std::cout << "remove:6" << std::endl; 
  //get rid of operators
  M[icopy].creators().pop_back();
  M[icopy].creators().pop_back();

  M[icopy].num_vertices() -= 1; 

  return Stilde.determinant(); 
}


/*implement remove vertex in W sector*/
double InteractionExpansion::Wremove_impl(const unsigned vertex, const bool compute_only_weight)
{

  unsigned int Msize = Msuper.matrix().rows();

  // the block we want to remove
  Eigen::MatrixXd Stilde = Msuper.matrix().block<2,2>(2*vertex, 2*vertex);

  if(compute_only_weight){
    return Stilde.determinant();
  }

  //swap with the back 
  unsigned pos = 2*vertex; 
  unsigned pos_j = Msize-2; 

  assert(pos <= pos_j); 

  Msuper.matrix().row(pos).swap(Msuper.matrix().row(pos_j));
  Msuper.matrix().row(pos+1).swap(Msuper.matrix().row(pos_j+1));

  Msuper.matrix().col(pos).swap(Msuper.matrix().col(pos_j));
  Msuper.matrix().col(pos+1).swap(Msuper.matrix().col(pos_j+1));

  std::swap(Msuper.creators()[pos], Msuper.creators()[pos_j]);
  std::swap(Msuper.creators()[pos+1], Msuper.creators()[pos_j+1]);
 

  //now perform fastupdate of M
  Msize -= 2; 
  Eigen::MatrixXd Qtilde(Msize, 2), Rtilde(2, Msize);

  if(Msize>0){
    Qtilde = Msuper.matrix().topRightCorner(Msize, 2) ; //block(i,j,rows,cols)
    Rtilde = Msuper.matrix().bottomLeftCorner(2, Msize) ; 
  }

  //std::cout << "remove:5" << std::endl; 

  Msuper.matrix().conservativeResize(Msize, Msize); 

  if(Msize>0)
    Msuper.matrix().noalias() -= Qtilde * Stilde.inverse() * Rtilde; 

  //std::cout << "remove:6" << std::endl; 
  //get rid of operators
  Msuper.creators().pop_back();
  Msuper.creators().pop_back();

  Msuper.num_vertices() -= 1; 

  return Stilde.determinant(); 
}
