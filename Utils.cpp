#include "Utils.h"

using namespace PCA;
  // Useful stats functions
float PCA::percentile(std::vector<float>& v, float pct) {
    if (v.size()==0) return 0;

    // copy to temp vector because order is changed
    std::vector<float> tmp(v.size());
    for(int i=0;i<v.size();++i) tmp[i]=v[i];
    
    int pctileIndex= static_cast<int> ((v.size()-1)*pct);
    if (pctileIndex>=v.size()) pctileIndex=v.size()-1;

    std::nth_element(tmp.begin(), tmp.begin()+pctileIndex, tmp.end());
    return tmp[pctileIndex];
  }
  
float PCA::median(std::vector<float>& v) {return percentile(v,0.5);}
  
  // compute median and median absolute deviation
float PCA::median_mad(std::vector<float> & v,float &mad) {
    float med=median(v);
    std::vector<float> medr(v.size());
    for(int i=0;i<v.size();++i) medr[i]=std::abs(v[i]-med);
    mad=1.4826*median(medr);
    return med;
  }

  // flag objects that are greater than n sigma from median
  // compute mad, median and mean of good objects
float PCA::median_mad_flag(std::vector<float> & v,float &mad,
			    float &mean,float sigma,
			 std::vector<int> &ind) {
    float med=median(v);
    std::vector<float> medr(v.size());
    for(int i=0;i<v.size();++i) medr[i]=std::abs(v[i]-med);
    mad=1.4826*median(medr);
    mean=0;
    int n_used=0;
    for(int i=0;i<v.size();++i) {
      if( std::abs(v[i]-med)/mad > sigma) {
	ind.push_back(i);
	continue;
      }
      mean+=v[i];
      n_used++;
    }
    mean/=n_used;
    return med;
}


void PCA::Tokenize(const string& str,
              std::vector<string>& tokens,
              const string& delimiters)
{
  // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

// Useful functions for fitting functions Legendre polynomials
DVector PCA::definePXY(int order, float x, float xmin, float xmax)
{
  DVector temp(order+1);
  float newx = (2.*x-xmin-xmax)/(xmax-xmin);
  temp[0] = 1.;
  if(order>0) temp[1] = newx;
  for(int i=2;i<=order;++i) {
    temp[i] = ((2.*i-1.)*newx*temp[i-1] - (i-1.)*temp[i-2])/i;
  }
  return temp;
}

// return matrix of x,y points 
void PCA::setPRow(int fitorder, Position<float> pos, 
	     const Bounds<float>& bounds, DVectorView prow)
{
  DVector px =
    definePXY(fitorder,pos.x,bounds.getXMin(),bounds.getXMax());
  DVector py =
    definePXY(fitorder,pos.y,bounds.getYMin(),bounds.getYMax());
  int pq = 0;
  for(int n=0;n<=fitorder;++n) {
    for(int p=n,q=n-p;q<=n;--p,++q) {
      prow(pq) = px[p]*py[q];
      ++pq;
    }
  }
  
}


void  PCA::meanRemove(FMatrix &m, FVector &mean,
			const std::vector<std::vector<bool> >  &missing)
{

  for(int i=0;i<m.ncols();++i) {
    float sum=0;
    
    int ntot=0;
    for(int j=0;j<m.nrows();++j) {
      if(missing[j][i]) continue;
      sum+=m(j,i);
      ntot++;
    }
    mean(i)=sum/ntot;
  }
  
  for(int i=0;i<m.ncols();++i) {
    m.col(i).addToAll(-mean(i));
  }

}
