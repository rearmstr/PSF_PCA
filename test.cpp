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
  std::string type=params.read<std::string>("type","mean");
  float exp_cut= params.read<float>("exp_cut",0.15);
  //bool skip61= params.read<bool>("skip61",true);

  cout<<"Settings..."<<endl;
  cout<<params<<endl;

  ifstream file(filename.c_str());  
  string name;
  while(file>>name) {
    
    Exposure exp(name,ccd);
    exp.setChipDivide(nx,ny);
    exp.setChipMax(xmax,ymax);
    if(skip61) exp.addSkip(61);
    bool suc=exp.readShapelet(dir+name+"/",nvar);
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
  if(type=="plin") nvar*=3;
  if(type=="pquad") nvar*=6;

  // Build the data matrix
  FMatrix dataM(nexp,nvar);
  for(int i=0;i<nexp;++i) {
    FVector med=exps[i].getVals(type);
    dataM.row(i)=med;
  }
  outputToFileF (dataM.transpose(), outname+"_data");
  
  // Remove mean from variables
  if(subtract_mean) {
    for(int i=0;i<nvar;++i) {
      
      FVector col=dataM.col(i);
      
      double sum=0.;
      for(int j=0;j<nexp;++j) {
        sum+=dataM(j,i);
      }
      for(int j=0;j<nexp;++j) {
        dataM(j,i)=dataM(j,i)-sum/nexp;
      }
      
    }
  }

  
  
  FDiagMatrix Svec(1);
  FMatrix U(1,1),V(1,1);
  
  if(nexp > nvar) {
    V.resize(nvar,nvar);
    Svec.resize(nvar);
    U.resize(nexp,nvar);
    U=dataM;
    SV_Decompose(U,Svec,V,true);


  }
  else {

    V.resize(nexp,nvar);
    Svec.resize(nexp);
    U.resize(nexp,nexp);
    V = dataM;
    SV_Decompose(V.transpose(),Svec,U.transpose());
  }

  // Check for outliers at the exposure level
  // if a single pca contributes more than exp_cut% then remove it and do the fit again
  int noutlier=0;
  for(int iexp=0;iexp<U.nrows();++iexp) {
    
    int sum=0;
    for(int ipca=0;ipca<U.nrows();++ipca) {
      double var=U(iexp,ipca)*U(iexp,ipca);
      sum+=var;

      if(var>exp_cut) {
        cout<<"Removing Exposure "<<exps[iexp].getLabel()
            <<" , has PC with fractional variance "<<var<<endl;;
        exps[iexp].setOutlier(1);
        noutlier++;
        break;
      }
    }
  }

  if(noutlier>0) {
    // Build the data matrix
    int nexp_cut=nexp-noutlier;
    std::ofstream oexp((outname+"_exp").c_str());
    FMatrix dataM(nexp_cut,nvar);
   
    int cur_exp=0;
    for(int i=0;i<nexp;++i) {
      if(exps[i].isOutlier()) continue;
      FVector med=exps[i].getVals(type);
      dataM.row(cur_exp)=med;
      oexp<<exps[i].getLabel()<<endl;
      cur_exp++;
    }
    outputToFileF (dataM.transpose(), outname+"_data");
  
    // Remove mean from variables
    if(subtract_mean) {
      for(int i=0;i<nvar;++i) {
        
        FVector col=dataM.col(i);
        
        double sum=0.;
        for(int j=0;j<nexp_cut;++j) {
          sum+=dataM(j,i);
        }
        for(int j=0;j<nexp_cut;++j) {
          dataM(j,i)=dataM(j,i)-sum/nexp_cut;
        }
        
      }
    }
    
    if(nexp_cut > nvar) {
      V.resize(nvar,nvar);
      Svec.resize(nvar);
      U.resize(nexp_cut,nvar);
      U=dataM;
      SV_Decompose(U,Svec,V,true);
    }
    else {
      
      V.resize(nexp_cut,nvar);
      Svec.resize(nexp_cut);
      U.resize(nexp_cut,nexp_cut);
      V = dataM;
      SV_Decompose(V.transpose(),Svec,U.transpose());
    }
  }
    
    
    
    outputToFileF (V.transpose(),  outname+"_vec");
    outputToFileF (U.transpose(), outname+"_coeff");
    outputToFileF (Svec.diag(), outname+"_singular"); 
    

  //cout<<med<<endl;

}
  
