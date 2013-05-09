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


void doSVD(FMatrix &data,int nvar,int nexp,FMatrix &U,
           FDiagMatrix &Svec,FMatrix &Vt)
{
  if(nexp > nvar) {
    Vt.resize(nvar,nvar);
    Svec.resize(nvar);
    U.resize(nexp,nvar);
    U=data;
    SV_Decompose(U,Svec,Vt,true);

  }
  else {

    Vt.resize(nexp,nvar);
    Svec.resize(nexp);
    U.resize(nexp,nexp);
    Vt = data;
    SV_Decompose(Vt.transpose(),Svec,U.transpose());
  }
}

int identifyOutliers(FMatrix &m,vector<bool> &outliers,float cut)
{
  int noutlier=0;
  outliers.resize(m.nrows());
  for(int iexp=0;iexp<m.nrows();++iexp) {
    
    int sum=0;
    outliers[iexp]=false;
    for(int ipca=0;ipca<m.ncols();++ipca) {
      double var=m(iexp,ipca)*m(iexp,ipca);
      sum+=var;

      if(var>cut)  {
        cout<<" Found outlier with % contribution: "<<var<<endl;
        outliers[iexp]=true;
        noutlier++;
        break;
      }
    }
  }
  return noutlier;
}


void meanRemove(FMatrix &m)
{
  for(int i=0;i<m.ncols();++i) {

    FVectorView col=m.col(i);
    double sum=m.col(i).sumElements();
    col.addToAll(-sum/m.nrows());
  }
}




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
  int max_outlier_iter=params.read<int>("max_outlier_iter",100);
  bool do_exp_rej=params.read<bool>("do_exp_rej",true);
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
  if(skip61 && ccd>61) nccd-=1;
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
  if(subtract_mean) meanRemove(dataM);

  
  // matrices for svd
  FDiagMatrix Svec(1);
  FMatrix U(1,1),Vt(1,1);
  doSVD(dataM,nvar,nexp,U,Svec,Vt);

  if(do_exp_rej) {
    // Check for outliers at the exposure level
    // if a single pca contributes more than exp_cut to the total
    // remove it and do the fit again
    int noutlier=0;
    int outlier_iter=0;
    do {
      cout<<"\nOutlier rejection iter "<<outlier_iter<<endl;
      //cout<<"Exposures remaining: "<<U.nrows()<<endl;
      vector<bool> outliers;
      noutlier=identifyOutliers(U,outliers,exp_cut);
      int nexp_cur=U.nrows();
      int iexp=0;
      cout<<"Found "<<noutlier<<" outliers"<<endl;

      for(int i=0;i<nexp;++i) {
        if(exps[i].isOutlier()) continue;
        
        if(outliers[iexp]) {
          cout<<"Removing Exposure "<<exps[iexp].getLabel()<<" outlier"<<endl;
          exps[i].setOutlier(1);
        }
        iexp++;
      }
      //cout<<"Found "<<iexp<<" that were not rejected "<<endl;
      if(noutlier>0) {
      int nexp_cut=nexp_cur-noutlier;
      //cout<<"Reducing size to "<<nexp_cut<<endl;
      dataM.setZero();
      dataM.resize(nexp_cut,nvar);
      
      int cur_exp=0;
      for(int i=0;i<nexp;++i) {
        if(exps[i].isOutlier()) continue;

        FVector med=exps[i].getVals(type);
        dataM.row(cur_exp)=med;
        cur_exp++;
      }

      // do I really want to write/overwrite at each stage of the rejection
      // iteration?  Other option is to keep dataM and not modify it
      outputToFileF (dataM.transpose(), outname+"_data");
      
      // Remove mean from variables
      if(subtract_mean) meanRemove(dataM);
      doSVD(dataM,nvar,nexp_cut,U,Svec,Vt);    
      }
      outlier_iter++;
    } while (noutlier>0 && outlier_iter-1<max_outlier_iter);
  }

  std::ofstream oexp((outname+"_exp").c_str());           
  for(int i=0;i<nexp;++i) {
    if(exps[i].isOutlier()) continue;
    oexp<<exps[i].getLabel()<<endl;
  }
    
    
  outputToFileF (Vt.transpose(),  outname+"_vec");
  outputToFileF (U.transpose(), outname+"_coeff");
  outputToFileF (Svec.diag(), outname+"_singular"); 
  
  
  }
  
