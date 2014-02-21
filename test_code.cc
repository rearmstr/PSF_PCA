#include "Detection.h"
#include "Data.h"
#include "Random.h"
#include <vector>
#include <iostream>
using namespace std;
using namespace PCA;
std::ostream* dbgout = 0;//&cout;
bool XDEBUG = true;

int main() 
{

  int n=10;
  int ndata=30;
  ran::UniformDeviate u;
  ran::GaussianDeviate gauss;
  vector<Detection<>*> dets;
  float max=10;
  Bounds<float> b(0,max,0,max);
  Cell c(ndata,-5,5,-5,5);
  int nexp=10;

  BinnedData *bd=new BinnedData(2,4,ndata,2);
  for(int iexp=0;iexp<nexp;++iexp) {
    Chip<> *chip2=new Chip<>(1,b);
    Chip<> *chip1=new Chip<>(2,b);
    for(int i=0;i<n;++i) {
      
      Detection<> *d=new Detection<>(u()*max,u()*max,ndata);
      
      for(int j=0;j<ndata;++j) {
	(*d)(j)=gauss();
      }
      
      dets.push_back(d);
      c.addDet(d);
      if(i>n/2)     chip1->addDet(d);
      else chip2->addDet(d);
    }
     
    // test that Cells work
    vector<float> param(1,3);
    c.calcVals("MeanClipped",param);
    
    stringstream ss;
    ss<<"exp"<<iexp;
    Exposure<> *exp=new Exposure<>(ss.str(),ndata,b);
    exp->addChip(0,chip1);
    exp->addChip(1,chip2);
    bd->addExp(exp);
  }
  
   
  vector<float> p(1,3);
  bd->calcData("Mean",p); 

  EMPcaMethod pca(3,bd,10,100,1e-4);
  pca.solve();
  pca.writeFits("test.fits");

  PcaMethod pca2(bd);
  pca2.solve();
  pca.writeFits("test2.fits");

//    vector<Cell> g;
//   g.push_back(c);
//   Cell g2=c;
//   Cell g3(c);
//   cout<<c.getData()<<" "<<g3.getData()<<endl;
//   vector<vector<Cell> > vg;
//   vg.push_back(g);
//   cout<<*vg[0][0].getDet(0)<<endl;

} 
 
