#include <vector>
#include "PCAObjects.h"
#include "myTypeDef.h"
#include "myClass.h"
#include "myIO.h"
#include "ConfigFile.h"
using namespace std;
using namespace PCA;
std::ostream* dbgout = 0;
bool XDEBUG = false;

int main(int argc,char*argv[])
{

  std::vector<Exposure> exps;

  
  ConfigFile params;
  params.setDelimiter("=");
  params.setInclude("+");
  params.setComment("#");
  params.load(argv[1]);
  for(int k=2;k<argc;k++) params.append(argv[k]);

  std::string filename= params.read<std::string>("file");
  int ccd= params.read<int>("ccd");
  int nvar= params.read<int>("nvar");
  int nx= params.read<int>("nx");
  int ny= params.read<int>("ny");
  float xmax= params.read<float>("xmax",2048.0);
  float ymax= params.read<float>("ymax",4096.0);
  std::string dir= params.read<std::string>("dir");
  int max_exp= params.read<int>("max_exp",-1);
  bool skip61= params.read<bool>("skip61",true);
  std::string outname= params.read<std::string>("outname");
  std::string grid_file= params.read<std::string>("grid_file");
  bool subtract_mean=params.read<bool>("subtract_mean",true);
  //bool skip61= params.read<bool>("skip61",true);

  cout<<"Settings..."<<endl;
  cout<<params<<endl;

  ifstream file(filename.c_str());  
  string name;
  while(file>>name) {
    
    Exposure exp(name,ccd,nvar);
    exp.setChipDivide(nx,ny);
    exp.setChipMax(xmax,ymax);
    if(skip61) exp.addSkip(61);
    bool suc=exp.readShapelet(dir+name+"/");
    if(suc) exps.push_back(exp);
    if(exps.size()>(max_exp-1) && max_exp>0) break;
  }
  
  // take the cell boundaries from the first chip of the first exposure
  std::vector<Bounds<float> > cb=exps[0][1]->getCellBounds();
  std::ofstream ogrid(grid_file.c_str());
  for(int j=0;j<ccd-1;++j) {
    for(int i=0;i<cb.size();++i) {
      ogrid<<cb[i].getXMin()<<" "<<cb[i].getYMin()<<" "<<cb[i].getXMax()<<" "<<cb[i].getYMax()<<endl;
    }
  }

  int nexp=exps.size();
  nvar*=nx*ny*(ccd-1);
  FMatrix dataM(nexp,nvar);
  for(int i=0;i<nexp;++i) {

    FVector med=exps[i].getMeanVals();
    //FVector med=exps[i].getMedianVals();
    dataM.row(i)=med;
  }
  
  FMatrix dataM_mr(nexp,nvar);
  dataM_mr.setZero();
  // Remove mean from variables
  if(subtract_mean) {
    for(int i=0;i<nvar;++i) {
      
      FVector col=dataM.col(i);
      
      double sum=0.;
      for(int j=0;j<nexp;++j) {
        sum+=dataM(j,i);
      }
      for(int j=0;j<nexp;++j) {
        dataM_mr(j,i)=dataM(j,i)-sum/nexp;
      }
      
    }
  }
  else {
    dataM_mr=dataM;
  }
  
  
  FDiagMatrix Svec(1);
  FMatrix U(1,1),V(1,1);
  
  if(nexp > nvar) {
    V.resize(nvar,nvar);
    Svec.resize(nvar);
    U.resize(nexp,nvar);
    
    U=dataM_mr;

    SV_Decompose(U,Svec,V,true);


  }
  else {
    cout<<"Resizing"<<endl;
    V.resize(nexp,nvar);
    Svec.resize(nexp);
    U.resize(nexp,nexp);
    cout<<"Assigining"<<endl;
    V = dataM_mr;
    cout<<"Decomposing"<<endl;
    SV_Decompose(V.transpose(),Svec,U.transpose());
  }


  outputToFileF (dataM.transpose(), outname+"_data");
  outputToFileF (V.transpose(),  outname+"_vec");
  outputToFileF (U.transpose(), outname+"_coeff");
  outputToFileF (Svec.diag(), outname+"_singular"); 
  

  //cout<<med<<endl;

}
  
