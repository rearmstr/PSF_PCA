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
  bool use_dash=params.read<bool>("use_dash",false);
  std::string prefix=params.read<std::string>("prefix","");
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
    bool suc=exp.readShapelet(dir+name+"/",nvar,use_dash,prefix+name);
    if(suc) exps.push_back(exp);
    if(exps.size()>(max_exp-1) && max_exp>0) break;
  }
  
  // take the cell boundaries from the first chip of the first exposure
  std::vector<Bounds<float> > cb=exps[0][1]->getCellBounds();
  std::ofstream ogrid(grid_file.c_str());
  
  int nccd=ccd;
  if(skip61) nccd-=1;
  for(int j=0;j<nccd;++j) {
    for(int i=0;i<cb.size();++i) {
      ogrid<<cb[i].getXMin()<<" "<<cb[i].getYMin()<<" "<<cb[i].getXMax()<<" "<<cb[i].getYMax()<<endl;
    }
  }

  int nexp=exps.size();
  // scale the number of variables to include the total number per exposure
  nvar*=nx*ny*nccd;

  if(type=="plin") nvar*=3; // fit a linear term to each cell
  if(type=="pquad") nvar*=6;// fit a quadratic function to each cell

  // Build the data matrix
  FMatrix dataM(nexp,nvar);
  for(int i=0;i<nexp;++i) {
    FVector med=exps[i].getVals(type);
    dataM.row(i)=med;
  }

  // output raw data file now before it is altered
  outputToFileF (dataM.transpose(), outname+"_data");
  
  // Remove mean from the data
  // probably can bemore efficient by using tmv operations
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

  
  // matrices for svd
  FDiagMatrix Svec(1);
  FMatrix U(1,1),Vt(1,1);
  
  if(nexp > nvar) {
    Vt.resize(nvar,nvar);
    Svec.resize(nvar);
    U.resize(nexp,nvar);
    U=dataM;
    SV_Decompose(U,Svec,Vt,true);


  }
  else {

    Vt.resize(nexp,nvar);
    Svec.resize(nexp);
    U.resize(nexp,nexp);
    Vt = dataM;
    SV_Decompose(Vt.transpose(),Svec,U.transpose());
  }

  // Check for outliers at the exposure level
  // if a single pca contributes more than exp_cut% then remove it and do the fit again
  // maybe we want to iterate this
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
    // go throuth the process again
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
      Vt.resize(nvar,nvar);
      Svec.resize(nvar);
      U.resize(nexp_cut,nvar);
      U=dataM;
      SV_Decompose(U,Svec,Vt,true);
    }
    else {
      
      Vt.resize(nexp_cut,nvar);
      Svec.resize(nexp_cut);
      U.resize(nexp_cut,nexp_cut);
      Vt = dataM;
      SV_Decompose(Vt.transpose(),Svec,U.transpose());
    }
  }
    
    
    
  outputToFileF (Vt.transpose(),  outname+"_vec");
  outputToFileF (U.transpose(), outname+"_coeff");
  outputToFileF (Svec.diag(), outname+"_singular"); 
  
  
  }
  
