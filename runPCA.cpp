#include <vector>
#include "PCAObjects.h"
#include "myTypeDef.h"
#include "myClass.h"
#include "myIO.h"
#include "ConfigFile.h"
#include <cassert>
#include "Log.h"
#include "NR.h"
using namespace std;
using namespace PCA;
std::ostream* dbgout = 0;
bool XDEBUG = false;

template<class T1,class T2>
void doSVD(T1 &data,int nvar,int nexp,T1 &U,
           T2 &Svec,T1 &Vt,std::vector<std::vector<bool> > &missing,
	   bool use_em=false,int npc=20,
	   int max_iter=1000,double tol=1e-6,bool do_missing=false)
  
{
  if(!use_em) {
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
  else {

    assert(npc<=nexp);

    // select a random set of exposures to initialize the solution
    vector<int> rands;
    while(rands.size()<npc) {
      int n=ran01()*nexp;
      if(find(rands.begin(),rands.end(),n)==rands.end()) {
	rands.push_back(n);
      }
    }

    // assume that we will always have more variables than exposures
    assert(nexp<nvar);
    Vt.resize(npc,nvar);
    Svec.resize(npc);
    U.resize(npc,npc);

    // use a random subset of data equal to the number of pcs to
    // give an initial solution via svd
    for(int i=0;i<npc;++i) {
      Vt.row(i)=data.row(i);
    }
    FILE_LOG(logINFO)<<"Initial Decomposition "<<endl;
    SV_Decompose(Vt.transpose(),Svec,U.transpose());
    T1 C(nvar,npc);
    C=Vt.transpose();
    T1 x(npc,nexp);
    for(int i=0;i<rands.size();++i) {
      x.col(i)=Svec(i)*U.col(i);
    }
    
    
    for(int iter=0;iter<max_iter;++iter) {

      // for no missing data can solve
      if(!do_missing) {
	T1 tmp=C.transpose()*C;
	x=C.transpose()*data.transpose()/tmp;
	
	T1 Cnew=data.transpose()*x.transpose()%(x*x.transpose());
	double diff=(Cnew-C).norm();
	diff/=(Cnew.nrows()*Cnew.ncols());
	FILE_LOG(logDEBUG1)<<"diff "<<diff<<" "<<tol<<endl;
	if(diff<tol) break;
	
	C=Cnew;
      
      }
      else {
	
	assert(missing.size()>1);
	int nvar_o=missing[0].size();
	x.setZero();
	// Solve for each exposure independently.  Could add openmp here later
	for(int iexp=0;iexp<nexp;iexp++) {
	  
	  bool cell_miss=false;
	  // loop over cells to see if any data are missing for this exposure
	  int nmiss=0;
	  for(int ivar=0;ivar<nvar;++ivar) {
	    int icell=ivar%nvar_o;

	    if (missing[iexp][icell]) {
	      cell_miss=true;
	      nmiss++;
	    }
	  }

	  FILE_LOG(logDEBUG2)<<"Exposure "<<iexp<<" missing "<<nmiss<<" variables"<<endl;

	  if(!cell_miss) {
	    
	    T1 tmp=C.transpose()*C;
	    x.subMatrix(0,npc,iexp,iexp+1)=C.transpose()*data.transpose().subMatrix(0,nvar,iexp,iexp+1)/tmp;
	  }
	  else {
	    // resshuffle C into Cnew so that the cells with missing data are in the lowest rows
	    // and cells with data are in the highest rows
	    T1 Cnew(nvar,npc);
	    int cur_missing=0;
	    int cur_here=0;
	    int nhere=nvar-nmiss;
	    T1 Y(nhere,1); // actual data values
	    T1 Dm(nmiss,1);// missing data
	    for(int ivar=0;ivar<nvar;++ivar) {
	      int icell=ivar%nvar_o;
	      	      
	      if (missing[iexp][icell]) {
		Cnew.subMatrix(cur_missing,cur_missing+1,0,npc)=C.subMatrix(ivar,ivar+1,0,npc);
		cur_missing++;
	      }
	      else {
		Cnew.subMatrix(nmiss+cur_here,nmiss+cur_here+1,0,npc)=C.subMatrix(ivar,ivar+1,0,npc);
		Y(cur_here,0)=data(iexp,ivar);
		
		cur_here++;
	      }
	    }

	    // solve for x with known data points
	    x.subMatrix(0,npc,iexp,iexp+1)=Y/Cnew.subMatrix(nmiss,nvar,0,npc);
	    
	    // calculate missing points from C and x
	    Dm=Cnew.subMatrix(0,nmiss,0,npc)*x.subMatrix(0,npc,iexp,iexp+1);

	    // fill in missing points into data matrix
	    cur_missing=0;
	    for(int ivar=0;ivar<nvar;++ivar) {
	      int icell=ivar%nvar_o;
	      if (missing[iexp][icell]) { 
		data(iexp,ivar)=Dm(cur_missing,0);
		cur_missing++;
	      }
	    }
	    
	  }
	}
	
	
	T1 Cnew=data.transpose()*x.transpose()%(x*x.transpose());

	double diff=(Cnew-C).norm();
	diff/=(Cnew.nrows()*Cnew.ncols());
	FILE_LOG(logINFO)<<"EM iteration "<<iter<<" diff "<<diff<<" "<<tol<<endl;
	if(diff<tol) break;
	C=Cnew;
      }
    }
      
      
    Svec.setZero();
    Svec.diag()=(x*x.transpose()).diag();
    U=(x.transpose()%Svec).subMatrix(0,npc,0,npc);
    Vt=C.transpose();

	//  // normalize the rows of U to one
	//     for(int i=0;i<U.nrows();++i) {
	//       double norm=U.row(i).normSq();
	//       Svec(i)*=norm;
	//       U.row(i)/=norm;
	//     }
      
  
  }
}



template<class T>
int identifyOutliers(T &m,vector<bool> &outliers,float cut)
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

template<class T>
void meanRemove(T &m)
{
  for(int i=0;i<m.ncols();++i) {

    double sum=m.col(i).sumElements();
    m.col(i).addToAll(-sum/m.nrows());
  }
}




int main(int argc,char*argv[])
{

  std::vector<Exposure<double> > exps;

  
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
  std::string image_dir= params.read<std::string>("image_dir","");
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
  int fit_order=params.read<int>("fit_order",-1);
  float sigma_clip=params.read<float>("sigma_clip",-1.);
  int logging=params.read<int>("logging",3);
  int npix=params.read<int>("npix",10);
  bool shapelet=params.read<bool>("shapelet",true);
  bool rm_zero=params.read<bool>("rm_zero",true);
  bool do_em=params.read<bool>("do_em",true);
  int max_iter=params.read<int>("max_iter",1000);
  int em_pc=params.read<int>("em_pc",20);
  float tol=params.read<float>("tol",1e-6);
  float add_missing=params.read<float>("add_missing",-1);

  FILELog::ReportingLevel() = FILELog::FromInt(logging);
  FILE_LOG(logINFO)<<"Settings...\n"<<params<<endl;

  ifstream file(filename.c_str());  
  string name;
  while(file>>name) {
    
    Exposure<double> exp(name,ccd);
    exp.setChipDivide(nx,ny);
    exp.setChipMax(xmax,ymax);
    if(skip61) exp.addSkip(61);
    bool suc;
    if(shapelet) {
      suc=exp.readShapelet(dir+name+"/",nvar,use_dash,prefix+name);
    }
    else {
      string fitsname=name;
      // erase the leading zeros in the exposure number
      // this is to fix some issues
      if(rm_zero) fitsname.erase(6,2);
      nvar=4*(npix+1)*(npix+1);
      suc=exp.readPixels(image_dir+fitsname+"/",npix,nvar,
			 dir+name+"/",use_dash,prefix+name);
    }
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
  
  int nvar_tot=nvar*nx*ny*nccd;
  std::vector<float> vparams(1,0.);
  
  if(type=="fit") {
    assert(fit_order>0);
    vparams[0]=fit_order; 
    nvar_tot*=(fit_order+1)*(fit_order+2)/2;
  }
  if(type=="mean_clip") {
    assert(sigma_clip>0);
    vparams[0]=sigma_clip; 
  }

  // artificially remove data from each exposure
  if(add_missing>0 && add_missing <1) {
    for(int i=0;i<nexp;++i) {
      exps[i].setMissing(add_missing);
    }
  }  
   

  

  // Build the data matrix
  DMatrix dataM(nexp,nvar_tot);
  std::vector<std::vector<bool> > missing(nexp,std::vector<bool>(exps[0].getCells(),false));
  bool hasMissing=false;
  for(int i=0;i<nexp;++i) {
    DVector med=exps[i].getVals(type,vparams);
    dataM.row(i)=med;
    missing[i]=exps[i].getMissing();
    FILE_LOG(logDEBUG)<<<<"exposure "<<i<<" missing cells"<<endl;
    for(int j=0;j<missing[i].size();++j) {
       FILE_LOG(logDEBUG)<<missing[i][j]<<" ";
    }
    FILE_LOG(logDEBUG)<<endl;
  }


  // output raw data file now before it is altered
  writeMatrix(dataM,outname+"_data");
  
  // Remove mean from the data
  // probably can bemore efficient by using tmv operations
  if(subtract_mean) meanRemove(dataM);

  
  // matrices for svd
  DDiagMatrix Svec(1);
  DMatrix U(1,1),Vt(1,1);
   hasMissing=true;
  cout<<"doing svd"<<endl;
  doSVD<DMatrix,DDiagMatrix>(dataM,nvar_tot,nexp,U,Svec,Vt,missing,do_em,em_pc,max_iter,tol,hasMissing);

  if(do_exp_rej) {
    // Check for outliers at the exposure level
    // if a single pca contributes more than exp_cut to the total
    // remove it and do the fit again
    int noutlier=0;
    int outlier_iter=0;
    do {
      FILE_LOG(logINFO)<<"\nOutlier rejection iter "<<outlier_iter<<endl;
      //FILE_LOG(logINFO)<<"Exposures remaining: "<<U.nrows()<<endl;
      vector<bool> outliers;
      noutlier=identifyOutliers<DMatrix>(U,outliers,exp_cut);
      int nexp_cur=U.nrows();
      int iexp=0;
      FILE_LOG(logINFO)<<"Found "<<noutlier<<" outliers"<<endl;

      for(int i=0;i<nexp;++i) {
        if(exps[i].isOutlier()) continue;
        
        if(outliers[iexp]) {
          FILE_LOG(logINFO)<<"Removing Exposure "<<exps[iexp].getLabel()<<" outlier"<<endl;
          exps[i].setOutlier(1);
        }
        iexp++;
      }
      FILE_LOG(logDEBUG)<<"Found "<<iexp<<" that were not rejected "<<endl;
      if(noutlier>0) {
      int nexp_cut=nexp_cur-noutlier;
      FILE_LOG(logDEBUG)<<"Reducing size to "<<nexp_cut<<endl;
      dataM.setZero();
      dataM.resize(nexp_cut,nvar_tot);
      
      int cur_exp=0;
      for(int i=0;i<nexp;++i) {
        if(exps[i].isOutlier()) continue;

        DVector med=exps[i].getVals(type,vparams);
        dataM.row(cur_exp)=med;
        cur_exp++;
      }

      // do I really want to write/overwrite at each stage of the rejection
      // iteration?  Other option is to keep dataM and not modify it
      writeMatrix(dataM,outname+"_data");

      // Remove mean from variables
      if(subtract_mean) meanRemove<DMatrix>(dataM);
      doSVD<DMatrix,DDiagMatrix>(dataM,nvar_tot,nexp_cut,U,Svec,Vt,missing,do_em,em_pc,max_iter,tol,hasMissing);    
      }
      outlier_iter++;
    } while (noutlier>0 && outlier_iter-1<max_outlier_iter);
  }

  std::ofstream oexp((outname+"_exp").c_str());           
  for(int i=0;i<nexp;++i) {
    if(exps[i].isOutlier()) continue;
    oexp<<exps[i].getLabel()<<endl;
  }
    
  writeMatrix(Vt,outname+"_vec");
  writeMatrix(U,outname+"_coeff");
  writeVector(Svec.diag(),outname+"_singular");
  
}

