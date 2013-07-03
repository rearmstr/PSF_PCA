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

// principal component decomposition, can do svd or em
template<class T1,class T2>
void doPCD(T1 &data,int nvar,int nexp,T1 &U,
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
      SV_Decompose(Vt.transpose(),Svec,U.transpose(),true);
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
	FILE_LOG(logDEBUG)<<"Doing EM for all exposures "<<endl;
	T1 tmp=C.transpose()*C;
	x=C.transpose()*data.transpose()/tmp;
	
	T1 Cnew=data.transpose()*x.transpose()%(x*x.transpose());
	double diff=(Cnew-C).norm();
	diff/=(Cnew.nrows()*Cnew.ncols());
	FILE_LOG(logDEBUG1)<<"diff "<<diff<<" "<<tol<<endl;
	if(diff<tol) break;
	FILE_LOG(logINFO)<<"EM iteration "<<iter<<" diff "
			 <<diff<<" "<<tol<<endl;
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
	  if(iter==0) {
	    FILE_LOG(logDEBUG)<<"Exposure "<<iexp<<" missing "
			      <<nmiss<<" variables"<<endl;
	  }

	  if(!cell_miss) {
	    
	    T1 tmp=C.transpose()*C;
	    x.subMatrix(0,npc,iexp,iexp+1)=C.transpose()*
	      data.transpose().subMatrix(0,nvar,iexp,iexp+1)/tmp;
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
		Cnew.subMatrix(cur_missing,cur_missing+1,
			       0,npc)=C.subMatrix(ivar,ivar+1,0,npc);
		cur_missing++;
	      }
	      else {
		Cnew.subMatrix(nmiss+cur_here,nmiss+cur_here+1,
			       0,npc)=C.subMatrix(ivar,ivar+1,0,npc);
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
      
       
    U.resize(nexp,npc);
    Svec.setZero();

    Svec.diag()=(C.transpose()*C).diag();
    for(int i=0;i<npc;++i) Svec(i)=std::sqrt(Svec(i));
    Vt=C.transpose()/Svec;
    U=x.transpose();

    //normalize the cols of U to one and scale Svec accordingly
    for(int i=0;i<U.ncols();++i) {
      double norm=U.col(i).normSq();
      Svec(i)*=std::sqrt(norm);
      U.col(i)/=std::sqrt(norm);
    }

    
    // reorder matrices in terms of largest eigenvalue
    tmv::Permutation p;
    DVector diag=Svec.diag();
    diag.sort(p,tmv::Descend);
    Svec.diag()=p*Svec.diag();
    Vt=p*Vt;
    U=U*p;
      
  }
}


// identify outliers where they are found by examining the
// components of a single exposure.  If any of the components
// are above cut label the exposure as bad.  It does not
// assume the vectors of an exposure are orthogonal
template<class T>
int identifyOutliers(T &m,vector<bool> &outliers,float cut)
{
  int noutlier=0;
  outliers.resize(m.nrows());
  for(int iexp=0;iexp<m.nrows();++iexp) {
    
    double sum=0;
    outliers[iexp]=false;
    for(int ipca=0;ipca<m.ncols();++ipca) {
      sum+=m(iexp,ipca)*m(iexp,ipca);
    }

    for(int ipca=0;ipca<m.ncols();++ipca) {
      double var=m(iexp,ipca)*m(iexp,ipca)/sum;
      
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


// remove the mean of each column
template<class T>
DVector meanRemove(T &m)
{
  DVector mean(m.ncols());
  for(int i=0;i<m.ncols();++i) {

    double sum=m.col(i).sumElements();
    m.col(i).addToAll(-sum/m.nrows());
    mean(i)=sum/m.nrows();
  }
  return mean;
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
  bool do_obj_rej=params.read<bool>("do_obj_rej",false);
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
  bool write_fits=params.read<bool>("write_fits",true);
  string read_fits=params.read<string>("read_fits","");

  FILELog::ReportingLevel() = FILELog::FromInt(logging);
  FILE_LOG(logINFO)<<"Settings...\n"<<params<<endl;

  ifstream file(filename.c_str());  
  string name;

  vector<string> exp_names;
  while(file>>name) {
    
    Exposure<double> exp(name,ccd);
    exp_names.push_back(name);
    exp.setChipDivide(nx,ny);
    exp.setChipMax(xmax,ymax);
    if(skip61) exp.addSkip(61);
    bool suc;
    if(shapelet) {
      suc=exp.readShapelet(dir+name+"/",nvar,do_em,use_dash,prefix+name);
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
  int nccd=ccd;
  if(skip61 && ccd>61) nccd-=1;

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

  std::vector<std::vector<bool> > missing(nexp,
					  std::vector<bool>(exps[0].getCells(),false));
  bool hasMissing=false;
  for(int i=0;i<nexp;++i) {
    DVector med=exps[i].getVals(type,vparams);
    dataM.row(i)=med;
    missing[i]=exps[i].getMissing();
    //FILE_LOG(logDEBUG)<<"exposure "<<i<<" missing cells"<<endl;
    
    for(int j=0;j<missing[i].size();++j) {
      //cout<<missing[i][j]<<" ";
       if(missing[i][j])hasMissing=true;
    }
    //cout<<endl;
  }

  // keep this around to write out full data matrix
  DMatrix original_data(nexp,nvar_tot);
  if(subtract_mean) original_data=dataM;

  FITS *fitfile=0;
  if(write_fits) {
    long naxis    =   2;      
    long naxes1[2] = { 1, 1 }; 
    fitfile=new FITS("!"+outname+".fits",USHORT_IMG , naxis , naxes1 );
  }
  
  // Remove mean from the data
  // probably can bemore efficient by using tmv operations
  DVector mean(dataM.ncols());
  if(subtract_mean) mean=meanRemove(dataM);

  
  // matrices for svd
  DDiagMatrix Svec(1);
  DMatrix U(1,1),Vt(1,1);
  
  doPCD<DMatrix,DDiagMatrix>(dataM,nvar_tot,nexp,U,Svec,Vt,missing,
			     do_em,em_pc,max_iter,tol,hasMissing);


  
  // identify outliers using the scores of the pcs.  Iterate until
  // no more exposures are removed
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

      exp_names.clear();

      for(int i=0;i<nexp;++i) {
	if(exps[i].isOutlier()) continue;
        
        if(outliers[iexp]) {
          FILE_LOG(logINFO)<<"Removing Exposure "<<exps[iexp].getLabel()<<" outlier"<<endl;
          exps[i].setOutlier(1);
        }
	exp_names.push_back(exps[i].getLabel());
        iexp++;
      }
      FILE_LOG(logDEBUG)<<"Found "<<iexp<<" that were not rejected "<<endl;
      if(noutlier>0) {
	int nexp_cut=nexp_cur-noutlier;
	FILE_LOG(logDEBUG)<<"Reducing size to "<<nexp_cut<<endl;
	dataM.setZero();
	dataM.resize(nexp_cut,nvar_tot);
	missing.resize(nexp_cut);
	int cur_exp=0;
	
	hasMissing=false;
	for(int i=0;i<nexp;++i) {
	  if(exps[i].isOutlier()) continue;
	  
	  DVector med=exps[i].getVals(type,vparams);
	  missing[i]=exps[i].getMissing();
	  for(int j=0;j<missing[i].size();++j) {
	    if(missing[i][j])hasMissing=true;
	  }
	  dataM.row(cur_exp)=med;
	  cur_exp++;
	  
	}
	
	
	// Remove mean from variables
	if(subtract_mean) {
	  
	  //  keep original
	  original_data.resize(nexp_cut,nvar_tot);
	  original_data=dataM;
	  mean.resize(dataM.ncols());
	  mean=meanRemove<DMatrix>(dataM);
	}
	
	doPCD<DMatrix,DDiagMatrix>(dataM,nvar_tot,nexp_cut,U,Svec,Vt,missing,
				   do_em,em_pc,max_iter,tol,hasMissing);    
	
	outlier_iter++;
      }
    } while (noutlier>0 && outlier_iter-1<max_outlier_iter);
  }
  
  
  if(do_obj_rej) {
    
    FILE_LOG(logINFO)<<"Doing star rejection on each exposure"<<endl;
    DMatrix dataR=U*Svec*Vt;
    if(subtract_mean) {
      for(int i=0;i<dataR.ncols();i++) dataR.col(i).addToAll(mean(i));
    }
    
    
    for(int i=0;i<nexp;++i) {
      if(exps[i].isOutlier()) continue;
      DVector data_exp=dataR.row(i);
      double obj_sigma_clip=3;

      exps[i].outlierReject(data_exp,obj_sigma_clip,type,vparams);
    }

    bool hasMissing=false;
    int cur_exp=0;
    for(int i=0;i<nexp;++i) {
      if(exps[i].isOutlier()) continue;
      
      DVector med=exps[i].getVals(type,vparams);
      dataM.row(cur_exp)=med;
      missing[i]=exps[i].getMissing();
      for(int j=0;j<missing[i].size();++j) {
	if(missing[i][j])hasMissing=true;
      }
      cur_exp++;
    }
      
    
    // Remove mean from variables
    if(subtract_mean) {
      original_data.resize(cur_exp,nvar_tot);
      original_data=dataM;
      mean.resize(dataM.ncols());
      mean=meanRemove<DMatrix>(dataM);

    }
    FILE_LOG(logINFO)<<"Redoing Decomposition"<<endl;
    doPCD<DMatrix,DDiagMatrix>(dataM,nvar_tot,cur_exp,U,Svec,Vt,missing,
			       do_em,em_pc,max_iter,tol,hasMissing); 
   
  }
  


  if(!write_fits) {
    writeMatrix(dataM,outname+"_data");
    writeMatrix(Vt,outname+"_vec");
    writeMatrix(U,outname+"_coeff");
    writeVector(Svec.diag(),outname+"_singular");
    std::ofstream oexp((outname+"_exp").c_str());           
   
    for(int i=0;i<exp_names.size();++i) {
      oexp<<exp_names[i]<<endl;
    }

    std::vector<Bounds<float> > cb=exps[0][1]->getCellBounds();
    std::ofstream ogrid(grid_file.c_str());
  
    
    for(int j=0;j<nccd;++j) {
      for(int i=0;i<cb.size();++i) {
	ogrid<<cb[i].getXMin()<<" "<<cb[i].getYMin()<<" "
	     <<cb[i].getXMax()<<" "<<cb[i].getYMax()<<endl;
      }
    }
  }
  else {

   
    

    // write the exposure information
    int nwvar=1;
    int nrows=exp_names.size();
    std::vector<string> colName(nwvar,"");
    std::vector<string> colForm(nwvar,"");
    std::vector<string> colUnit(nwvar,"");
    colName[0] = "exposure";
    colForm[0] = "16A";
    colUnit[0] = "";
    Table* newTable = fitfile->addTable("exps",nrows,colName,colForm,colUnit);
    newTable->column(colName[0]).write(exp_names,1);
 // write the header information
    newTable->addKey("ccd",ccd,"number of ccds");
    newTable->addKey("nvar",nvar,"number of variabls");
    newTable->addKey("nx",nx,"cells in x direction");
    newTable->addKey("ny",ny,"cells in y direction");
    newTable->addKey("xmax",xmax,"maximum ccd x");
    newTable->addKey("ymax",ymax,"maximum ccd y");
    newTable->addKey("rm_mean",subtract_mean,"remove mean");
    newTable->addKey("type",type,"cell estimation");
    newTable->addKey("exp_cut",exp_cut,"exposure cut");
    newTable->addKey("order",fit_order,"fit order within cell");
    newTable->addKey("clip",sigma_clip,"sigma clip within cell");
    if(do_em) {
      newTable->addKey("em_pc",em_pc,"EM PCs");
      newTable->addKey("em_tol",tol,"EM tol erance");
      newTable->addKey("em_iter",max_iter,"EM maximum iterations");
    }
    


    // write the grid information
    vector<double> lx,ux,ly,uy;
    std::vector<Bounds<float> > cb=exps[0][1]->getCellBounds();
    for(int i=0;i<cb.size();++i) {
      lx.push_back(cb[i].getXMin());
      ly.push_back(cb[i].getYMin());
      ux.push_back(cb[i].getXMax());
      uy.push_back(cb[i].getYMax());
    }
    
    nwvar=4;
    nrows=cb.size();
    colName.resize(nwvar);
    colForm.resize(nwvar);
    colUnit.resize(nwvar);

    colName[0] = "lower_x";
    colName[1] = "lower_y";
    colName[2] = "upper_x";
    colName[3] = "upper_y";
    colForm[0] = "1E";
    colForm[1] = "1E";
    colForm[2] = "1E";
    colForm[3] = "1E";
    colUnit[0] = "pixels";
    colUnit[1] = "pixels";
    colUnit[2] = "pixels";
    colUnit[3] = "pixels";
    Table* newTable2 = fitfile->addTable("grid",nrows,colName,colForm,colUnit);
    newTable2->column(colName[0]).write(lx,1);
    newTable2->column(colName[1]).write(ly,1);
    newTable2->column(colName[2]).write(ux,1);
    newTable2->column(colName[3]).write(uy,1);

    // write the matrices
    if(subtract_mean) {
      writeMatrixToFits<DMatrix>(fitfile,original_data,"data");
      writeMatrixToFits<DMatrix>(fitfile,dataM,"data_mr");
    }
    else {
      writeMatrixToFits<DMatrix>(fitfile,dataM,"data");
    }
    writeMatrixToFits<DMatrix>(fitfile,Vt,"vec");
    writeMatrixToFits<DMatrix>(fitfile,U,"coeff");
    writeVectorToFits<DVector>(fitfile,Svec.diag(),"singular");
        
  }
  
 
}

