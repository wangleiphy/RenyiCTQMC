#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

class histogram
{
public:

  histogram():hist_(){}
  histogram(unsigned int N):hist_(N, 0){}

  void resize(unsigned int N){
      hist_.resize(N, 0); 
  }
  
  unsigned long &operator[](unsigned int n){return hist_[n];}
  const unsigned long &operator[](unsigned int n) const{return hist_[n];}
  unsigned int size() const{return hist_.size();}
  
  unsigned int max_index() const
  { 
    unsigned int max_index=0; 
    double max=0; 
    for(unsigned int i=0;i<hist_.size();++i){
      if(max<hist_[i]){
        max=hist_[i];
        max_index=i;
      }
    }
    return max_index;
  }
  
  unsigned int top_index() const
  { 
    unsigned int top_index=0;
    for(unsigned int i=0;i<hist_.size();++i){
      if(hist_[i]!=0){
        top_index=i;
      }
    }
    return top_index;
  }
  
  double max(const unsigned int index) const
  { 
    double max=0; 
    for(unsigned int i=0;i<index;++i){
      if(max<hist_[i]){
        max=hist_[i];
      }
    }
    return max;
  }
  
  double average(const unsigned int index) const{ 
    double average=0; 
    for(unsigned int i=0;i<index;++i){
      average+=hist_[i];
    }
    return average/index;
  }
  
  bool is_flat(const unsigned int index)const{return max(index)*0.8<average(index);}
  
  void clear()
  {
    for(unsigned int i=0;i<hist_.size();++i){
      hist_[i]=0;
    }
  }

private:
  
  std::vector<unsigned long> hist_;
};



#endif 
