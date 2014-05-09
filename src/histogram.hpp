#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

class histogram
{
public:

  histogram():hist_(){}
  histogram(unsigned int N):hist_(N, 0){}

  void resize(unsigned int N, const double val){
      hist_.resize(N, val); 
  }
  
  double &operator[](unsigned int n){return hist_[n];}
  const double &operator[](unsigned int n) const{return hist_[n];}
  unsigned size() const{return hist_.size();}

  /*
  double mean() const
  {
      double up= 0.0; 
      double down = 0.0; 
      for(unsigned i=0;i<hist_.size();++i){
          up += i * hist_[i]; 
          down += hist_[i]; 
      }
      return up/down; 
  } 
  
  unsigned max_index() const
  { 
    unsigned max_index=0; 
    double max=0; 
    for(unsigned i=0;i<hist_.size();++i){
      if(max<hist_[i]){
        max=hist_[i];
        max_index=i;
      }
    }
    return max_index;
  }
  */
  
  unsigned top_index() const// the last nonvanishing element past one past one 
  { 
    unsigned top_index=0;
    for(unsigned i=0;i<hist_.size();++i){
      if(hist_[i]>0.){
        top_index=i;
      }
    }
    return top_index+1;
  }
  
  double max(const unsigned index) const
  { 
    double max=0; 
    for(unsigned int i=0;i<index;++i){
      if(max<hist_[i]){
        max=hist_[i];
      }
    }
    return max;
  }
  
  double average(const unsigned index) const{ 
    double average=0; 
    for(unsigned i=0;i<index;++i){
      average+=hist_[i];
    }
    return average/index;
  }
  
  bool is_flat(const unsigned index)const{return max(index)*0.8<average(index);}
  
  void clear()
  {
    for(unsigned i=0;i<hist_.size();++i){
      hist_[i]=0.;
    }
  }

private:
  
  std::vector<double> hist_;
};



#endif 
