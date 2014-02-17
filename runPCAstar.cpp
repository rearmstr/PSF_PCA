#include <vector>
#include "PCAObjects.h"
#include "myTypeDef.h"
#include "myClass.h"
#include "myIO.h"
#include "ConfigFile.h"
#include <cassert>
#include "Log.h"
#include "NR.h"
#include "omp.h"
#include "Random.h"
#include <complex>
using namespace std;
using namespace PCA;
using namespace ran;
std::ostream* dbgout = 0;
bool XDEBUG = false;



void MakePsi(DMatrix& psi, CDVectorView z, int order)
{
    // For p>=q:
    //
    // psi_pq = (pi p! q!)^-1/2 z^m exp(-r^2/2) K_pq(r^2)
    //
    // where K_pq(r^2) satisfies the recursion relation:
    // K_pq = (r^2 - (N-1)) K_p-1,q-1 - (p-1)(q-1) K_p-2,q-2
    // with K_p0 = 1
    //
    // The recursion for psi_pq can then be derived as:
    //
    // psi_pq = (pq)^-1/2 (r^2-(N-1)) psi_p-1,q-1
    //          - sqrt( (p-1)(q-1)/(pq) ) psi_p-2,q-2
    //
    // psi_p0 = (z/sqrt(p)) psi_p-1,0
    // 
    // psi_00 = 1/sqrt(pi) exp(-r^2/2)

    Assert(int(psi.rowsize()) >= (order+1)*(order+2)/2);
    Assert(psi.colsize() == z.size());
    Assert(psi.iscm());
    Assert(!psi.isconj());

    const double invsqrtpi = 1./sqrt(PI);

    // Setup rsq, z vectors and set psi_00
    DVector rsq(z.size());
    double* rsqit = rsq.ptr();
    double* psi00it = psi.ptr();
    const std::complex<double>* zit = z.cptr();
    const int zsize = z.size();
    for(int i=0;i<zsize;++i) {
        rsqit[i] = std::norm(zit[i]);
        psi00it[i] = invsqrtpi * exp(-(rsqit[i])/2.);
    }


    DVector zr = z.realPart();
    DVector zi = z.imagPart();
    if (order >= 1) {
        // Set psi_10
        // All m > 0 elements are intrinsically complex.
        // However, we are fitting to a real intensity pattern
        // with complex coefficients of the complex shapelets.
        // Since psi_pq = psi_qp* (* = complex conjugate),
        // we know that b_pq must be b_qp*
        // b_pq psi_pq + b_pq* psi_pq* = 2 b_pqr psi_pqr - 2 b_pqi psi_pqi
        // So the values we want for the real fitter are
        // really 2 real(psi_pq) and -2 imag(psi_pq)
        // Putting the 2's here carries through to the rest of the 
        // elements via the recursion.
        psi.col(1) = 2. * DiagMatrixViewOf(zr) * psi.col(0);
        psi.col(2) = -2. * DiagMatrixViewOf(zi) * psi.col(0);
    }
    for(int N=2,k=3;N<=order;++N) {
        // Set psi_N0
        // The signs of these are not what you naively think due to 
        // the +2, -2 discussed above.  You just have to follow through
        // what the complex psi_N0 is, and what value is stored in the
        // psi_N-1,0 location, and what needs to get stored here.
        psi.col(k) = sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N);
        psi.col(k) += sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N+1);
        psi.col(k+1) = -sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N);
        psi.col(k+1) += sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N+1);
        k+=2;

        // Set psi_pq with q>0
        // The rsq part of this calculation can be done in batch, which 
        // speeds things up a bit.
        psi.colRange(k,k+N-1) =
            DiagMatrixViewOf(rsq) * psi.colRange(k-2*N-1,k-N-2);
        psi.colRange(k,k+N-1) -= (N-1.) * psi.colRange(k-2*N-1,k-N-2);
        // The other calculation steps are different for each component:
        for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
            double pq = p*q;
            if (m==0) {
                psi.col(k) /= sqrt(pq);
                if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
                ++k;
            } else {
                psi.colRange(k,k+2) /= sqrt(pq);
                if (q > 1)
                    psi.colRange(k,k+2) -= 
                        sqrt(1.-(N-1.)/pq)*psi.colRange(k+2-4*N,k+4-4*N);
                k+=2;
            }
        }
    }
}

// principal component decomposition, can do svd or em
template<class T1,class T2,class T3>
void doPCD(T1 &data,int nvar,int nvar_single,int nexp,T3 &U,
           T2 &Svec,T3 &Vt,const std::vector<std::vector<bool> > &missing,
	   bool use_em,int npc,
	   int max_iter,int min_iter,double tol,bool do_missing,T3 &C,T3 &x)
  
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


    if(npc>nexp) {
      FILE_LOG(logINFO)<<"Not enough exposures "<<nexp<<" to have "
		       <<npc<<" principal components"<<endl;
      exit(1);
    }
      

    // assume that we will always have more variables than exposures
    if(nexp>nvar) {
      FILE_LOG(logINFO)<<"Not enough variables for "
		       <<npc<<" principal components increase the variables"<<endl;
      exit(1);
    }

    Vt.resize(npc,nvar);
    Svec.resize(npc);
    U.resize(npc,npc);

    //T3 C(nvar,npc);
    //T3 x(npc,nexp);
    bool use_old=true;
//     if(C.ncols()!=npc && C.nrows()!=nvar ) {
//       && x.ncols()!=nexp && x.nrows()!=npc) {
//       C.resize(nvar,npc);
//       x.resize(npc,nexp);
//       C.setZero();
//       x.setZero();
//       use_old=false;
//     }
    if(C.ncols()==1 && C.nrows()==1 
       && x.ncols()==1 && x.nrows()==1) {
      
      C.resize(nvar,npc);
      x.resize(npc,nexp);
      C.setZero();
      x.setZero();

      use_old=false;
    }
    else {
      if( nexp<x.ncols() ) {
	x.resize(npc,nexp);
      }
    }

      


    int nall=0;
    for(int i=0;i<missing.size();++i) {
      bool miss=false;
      for(int j=0;j<missing[i].size();++j) {
	if(missing[i][j]) miss=true;
      }
      if (!miss) nall++;
    }
    if(!do_missing) nall=npc;
    
    if(nall>=npc && !use_old) {

      // use a random subset of data equal to the number of pcs to
      // give an initial solution via svd
      // select a random set of exposures to initialize the solution


      vector<int> rands;
      while(rands.size()<npc) {
	int n=ran01()*nexp;
	if(find(rands.begin(),rands.end(),n)==rands.end()) {
	  bool miss=false;
	  for(int i=0;i<missing[n].size();++i) {
	    if(missing[n][i]) miss=true;
	  }
	  
	  if(miss) continue;
	  rands.push_back(n);
	}
      }
    
      
      for(int i=0;i<npc;++i) {
	FILE_LOG(logDEBUG)<<"Selected "<<rands[i]<<" for initial svd"<<endl;
	Vt.row(i)=data.row(rands[i]);
      }
      FILE_LOG(logDEBUG)<<"Initial Decomposition "<<endl;
      FILE_LOG(logDEBUG)<<"Vt: "<<Vt<<endl;
      SV_Decompose(Vt.transpose(),Svec,U.transpose(),true);
      
      
      C=Vt.transpose();
      
      for(int i=0;i<rands.size();++i) {
	x.col(i)=Svec(i)*U.col(i);
      }
    }
    else if (!use_old) {
      FILE_LOG(logINFO)<<"Not enough full exposures.  Using random starting matrix"<<endl;
      GaussianDeviate gaus;
      for(int i=0;i<nvar;++i) {
	for(int j=0;j<npc;++j) {
	  C(i,j)=gaus()*0.01;
	}
      }
      
      for(int i=0;i<npc;++i) {
	for(int j=0;j<npc;++j) {
	  x(i,j)=gaus()*0.01;
	}
      }
    }
      
    
    FILE_LOG(logDEBUG1)<<"Initial C:"<<C<<endl;
    FILE_LOG(logDEBUG1)<<"Initial x:"<<x<<endl;
    double prev_diff=1e10;
    for(int iter=0;iter<max_iter;++iter) {

      // for no missing data can solve
      if(!do_missing ) {
	FILE_LOG(logDEBUG)<<"Doing EM for all exposures "<<endl;
	T3 tmp=C.transpose()*C;
	x=C.transpose()*data.transpose()/tmp;
	
	T3 Cnew=data.transpose()*x.transpose()%(x*x.transpose());
	double diff=(Cnew-C).norm();
	diff/=(Cnew.nrows()*Cnew.ncols());
	FILE_LOG(logDEBUG1)<<"diff "<<diff<<" "<<tol<<endl;

	FILE_LOG(logINFO)<<"EM iteration "<<iter<<" diff "
			 <<diff<<" "<<std::abs(diff-prev_diff)/prev_diff
			 <<" "<<tol<<endl;

	C=Cnew;
	if(iter>min_iter && 
	   (diff<tol || std::abs(diff-prev_diff)/(prev_diff)<tol) ) break;
	prev_diff=diff;


      
      }
      else {
	
	//assert(missing.size()>1);
	//x.setZero();
	// Solve for each exposure independently.  Could add openmp here later

#pragma omp parallel for shared(x,C,data)

	for(int iexp=0;iexp<nexp;iexp++) {

	  T3 xt(npc,nexp);
	  bool cell_miss=false;
	  // loop over cells to see if any data are missing for this exposure
	  int nmiss=0;
	  for(int ivar=0;ivar<nvar;++ivar) {
	    int icell=ivar/nvar_single;

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
	    
	    T3 tmp=C.transpose()*C;
	    x.subMatrix(0,npc,iexp,iexp+1)=C.transpose()*
	      data.transpose().subMatrix(0,nvar,iexp,iexp+1)/tmp;
	  }
	  else {
	    
	    // resshuffle C into Cnew so that the cells with missing data are in the lowest rows
	    // and cells with data are in the highest rows
	    T3 Cnew(nvar,npc);
	    Cnew.setZero();
	    int cur_missing=0;
	    int cur_here=0;
	    int nhere=nvar-nmiss;
	    T3 Y(nhere,1); // actual data values
	    T3 Dm(nmiss,1);// missing data
	    Y.setZero();
	    Dm.setZero();
	    //FILE_LOG(logDEBUG1)<<"Missing "<<nmiss<<" here "<<nhere<<endl;
	    for(int ivar=0;ivar<nvar;++ivar) {
	      int icell=ivar/nvar_single;

	      //FILE_LOG(logDEBUG1)<<"Variable "<<ivar<<" cell: "<<icell<<endl;
	      if (missing[iexp][icell]) {
		//FILE_LOG(logDEBUG1)<<" is missing "<<data(iexp,ivar)<<endl;
		
		Cnew.subMatrix(cur_missing,cur_missing+1,
			       0,npc)=C.subMatrix(ivar,ivar+1,0,npc);
		cur_missing++;
	      }
	      else {
		//FILE_LOG(logDEBUG1)<<" is not missing "<<data(iexp,ivar)<<endl;
	
		Cnew.subMatrix(nmiss+cur_here,nmiss+cur_here+1,
			       0,npc)=C.subMatrix(ivar,ivar+1,0,npc);
		
		Y(cur_here,0)=data(iexp,ivar);
		
		cur_here++;
	      }
	    }
	    
	    // solve for x with known data points
	    FILE_LOG(logDEBUG1)<<"data_good "<<Y<<endl;
	    FILE_LOG(logDEBUG1)<<"Cnew_good "<<Cnew.subMatrix(nmiss,nvar,0,npc)<<endl;

	    // This needs to be flagged
	      x.subMatrix(0,npc,iexp,iexp+1)=Y/Cnew.subMatrix(nmiss,nvar,0,npc);


	    FILE_LOG(logDEBUG1)<<" x_good:= "<<x.subMatrix(0,npc,iexp,iexp+1)<<endl;
	    
	    // calculate missing points from C and x
	    Dm=Cnew.subMatrix(0,nmiss,0,npc)*x.subMatrix(0,npc,iexp,iexp+1);
	    FILE_LOG(logDEBUG1)<<" Missing Data "<<Dm<<endl;

	    // fill in missing points into data matrix
	    cur_missing=0;
	    for(int ivar=0;ivar<nvar;++ivar) {
	      int icell=ivar/nvar_single;
	      if (missing[iexp][icell]) { 
		FILE_LOG(logDEBUG1)<<" Setting missing data  "<<ivar<<" "<<icell<<" "<<Dm(cur_missing,0)<<endl;

		data(iexp,ivar)=Dm(cur_missing,0);
		cur_missing++;
	      }
	    }
	    
	  }
	}


	FILE_LOG(logDEBUG1)<<"Update data "<<data<<endl;
	FILE_LOG(logDEBUG1)<<"Update x "<<x<<endl;
	T3 Cnew=data.transpose()*x.transpose()%(x*x.transpose());
	FILE_LOG(logDEBUG1)<<"Update C "<<Cnew<<endl;

	double diff=((Cnew*x).transpose()-data).norm();
	diff/=(data.nrows()*data.ncols());
	
	FILE_LOG(logINFO)<<"EM iteration "<<iter<<" diff "
			 <<diff<<" "<<std::abs(diff-prev_diff)/prev_diff
			 <<" "<<tol<<endl;
	
	C=Cnew;
	if(iter>min_iter && 
	   (diff<tol || std::abs(diff-prev_diff)/(prev_diff)<tol) ) break;
	prev_diff=diff;

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
  outliers.clear();
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
         FILE_LOG(logINFO)<<" Found outlier with % contribution: "<<var<<endl;
        outliers[iexp]=true;
        noutlier++;
        break;
      }
    }
  }
  return noutlier;
}


// remove the mean of each column
template<class T,class T2>
DVector meanRemove(T &m,int nvar,const T2 &missing)
{
  DVector mean(m.ncols());

  for(int i=0;i<m.ncols();++i) {
    int icell=i/nvar;
    double sum=0;
    int ntot=0;
    for(int j=0;j<m.nrows();++j) {
      if(missing[j][icell] || m(j,i)<-100) continue;
      sum+=m(j,i);
      ntot++;
    }
    mean(i)=sum/ntot;
  }
      
  for(int i=0;i<m.ncols();++i) {
    m.col(i).addToAll(-mean(i));
  }

  return mean;
}



double mean(DVector& v) {
  if (v.size()==0) return 0;
  return v.sumElements()/v.size();
}

double covariance(DVector& v1, DVector& v2) {
  int N = v1.size();
  if (N != v2.size()) {
    cerr << "covariance vector size incompatible";
    return 0;
  }
  if (N == 0) return 0;
  double cov = 0, mean1 = mean(v1), mean2 = mean(v2);
  for (int i=0; i<N; i++) cov += (v1[i]-mean1) * (v2[i]-mean2);
  cov = (cov / N);
  return cov;
}


double variance(DVector& v) {
  return covariance(v, v);
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
  std::string pcname= params.read<std::string>("pcfile","pcdata");
  int ccd= params.read<int>("ccd");
  int nvar= params.read<int>("nvar");
  int nx= params.read<int>("nx",100);
  int ny= params.read<int>("ny",100);
  std::string dir= params.read<std::string>("dir");
  std::string image_dir= params.read<std::string>("image_dir","");
  int max_exp= params.read<int>("max_exp",-1);
  bool skip61= params.read<bool>("skip61",true);
  std::string outname= params.read<std::string>("outname");
  bool subtract_mean=params.read<bool>("subtract_mean",true);
  std::string type=params.read<std::string>("type","mean");
  float exp_cut= params.read<float>("exp_cut",0.15);
  bool use_dash=params.read<bool>("use_dash",false);
  std::string prefix=params.read<std::string>("prefix","");
  int max_outlier_iter=params.read<int>("max_outlier_iter",100);
  bool do_exp_rej=params.read<bool>("do_exp_rej",false);
  bool do_obj_rej=params.read<bool>("do_obj_rej",false);
  int fit_order=params.read<int>("fit_order",-1);
  int shape_order=params.read<int>("shape_order",10);
  float sigma_clip=params.read<float>("sigma_clip",-1.);
  float obj_sigma_clip=params.read<float>("obj_sigma_clip",3);
  int logging=params.read<int>("logging",3);
  int npix=params.read<int>("npix",10);
  bool shapelet=params.read<bool>("shapelet",true);
  bool rm_zero=params.read<bool>("rm_zero",true);
  bool do_em=params.read<bool>("do_em",true);
  int max_iter=params.read<int>("max_iter",1000);
  int min_iter=params.read<int>("min_iter",100);
  int em_pc=params.read<int>("em_pc",20);
  float tol=params.read<float>("tol",1e-6);
  float add_missing=params.read<float>("add_missing",-1);
  bool use_missing=params.read<bool>("use_missing",false);
  bool write_fits=params.read<bool>("write_fits",true);
  bool write_obj=params.read<bool>("write_obj",false);
  string read_fits=params.read<string>("read_fits","");
  string suffix=params.read<string>("suffix","psf.fits");
  bool add_size=params.read<bool>("add_size",false);
  bool read2=params.read<bool>("read2",false);
  float readmax=params.read<float>("readmax",1);
  string use_dir=params.read<string>("use_dir",".");
  int threads=params.read<int>("threads",1);
  int seed=params.read<int>("seed",11111);

  double rmax=params.read<double>("rmax",5);
  int shapestart=3;

  omp_set_num_threads(threads);

  if ( seed==0) seed=unsigned (std::time(0));
  std::srand(seed);


  if(add_size) {
    nvar+=1;
    //shapestart--;
  }
  FILELog::ReportingLevel() = FILELog::FromInt(logging);
  FILE_LOG(logINFO)<<"Settings...\n"<<params<<endl;

  ifstream file(filename.c_str());  
  string name;

  // take the cell boundaries from the first chip of the first exposure

  vector<string> exp_names;
  while(file>>name) {
    
    Exposure<double> exp(name,ccd,shapestart);
    exp.setChipDivide(1,1);
    exp.setChipMax(2048,4096);

    if(skip61) exp.addSkip(61);
    bool suc;
    if(shapelet) {
      suc=exp.readShapelet(dir+name+"/",nvar,add_size,
			 do_em,use_dash,suffix, prefix+name,readmax
			 ,use_dir);
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
  
  
  
  // scale the number of variables to include the total number per exposure
  
  int nexp=0;
  for(int i=0;i<exps.size();++i) {
    nexp+=exps[i].getNGood();
  }
  FILE_LOG(logINFO)<<"Found "<<nexp<<" stars for "<<exps.size()<<" exposures"<<endl;
  // Build the data matrix
  DMatrix data(nexp,nvar);
  std::vector<std::vector<bool> > 
    missing(nexp,std::vector<bool>(1,false));

  bool hasMissing=false;
  int notMissing=0;
  int current_star=0;
  map<int,map<int,vector<int> > > exp_ccd;
  int nused_ccd=0;
  for(int i=0;i<exps.size();++i) {

    std::map<int,Chip<double>*>::const_iterator iter=exps[i].chips.begin();
    int ichip=0;
    for(; iter!=exps[i].chips.end();++iter,++ichip) {
      // There is only one cell here
      for(int istar=0;istar < iter->second->getCell(0)->getNDet();++istar) {
	DVector star(nvar);
	for(int ivar=0;ivar<nvar;++ivar) {

	  star(ivar)=iter->second->getCell(0)->getDet(istar)->getVal(ivar);
	}
	
	data.row(current_star)=star;
	exp_ccd[i][ichip].push_back(current_star);
	current_star++;
      }
    }
    if(i==0) nused_ccd=ichip;
  }

  DMatrixView dataM=data.rowRange(0,nexp);
 
  DMatrix original_data(nexp,nvar);
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
  if(subtract_mean) mean=meanRemove(dataM,nvar,missing);


  // Try to reject bad

  
  if(hasMissing && !use_missing) hasMissing=false;
  // matrices for svd
  DDiagMatrix Svec(1);
  DMatrix U(1,1),Vt(1,1);
  DMatrix C(1,1),x(1,1);
  doPCD<DMatrixView,DDiagMatrix,DMatrix>(dataM,nvar,1,nexp,U,Svec,Vt,missing,
					 do_em,em_pc,max_iter,min_iter,tol,hasMissing,C,x);

  
  if(hasMissing && do_em && use_missing) {
    for(int i=0;i<dataM.ncols();++i) {
      for(int j=0;j<dataM.nrows();++j) {
	if(subtract_mean) {
	  original_data(j,i)=dataM(j,i)+mean(i);
	}
	else {
	  original_data(j,i)=dataM(j,i);
	}
	
      }
    }
  }


  DMatrix dataR=U*Svec*Vt;
  if(subtract_mean) {
    for(int i=0;i<dataR.ncols();i++) dataR.col(i).addToAll(mean(i));
  }

    

  // write the matrices
  if(subtract_mean) {
    writeMatrixToFits<DMatrix>(fitfile,original_data,"data");
    writeMatrixToFits<DMatrixView>(fitfile,dataM,"data_mr");
    writeVectorToFits<DVector>(fitfile,mean,"mean");
  }
  else {
    writeMatrixToFits<DMatrixView>(fitfile,dataM,"data");
  }
  writeMatrixToFits<DMatrix>(fitfile,dataR,"dataR");
  writeMatrixToFits<DMatrix>(fitfile,Vt,"vec");
    writeMatrixToFits<DMatrix>(fitfile,U,"coeff");
    writeVectorToFits<DVector>(fitfile,Svec.diag(),"singular");
    //writeVVectorToFits<bool>(fitfile,missing,"missing");
    double var=0;
    DVector cumvar(Svec.diag().size());
    for(int i=0;i<Svec.diag().size();++i) {
      var+=Svec(i,i)*Svec(i,i);
    }
    double tot=0;
    for(int i=0;i<Svec.diag().size();++i) {
      tot+=Svec(i,i)*Svec(i,i)/var;
      cumvar(i)=tot;
    }
    writeVectorToFits<DVector>(fitfile,cumvar,"cumvar");

    
    // Check CCD to CCD variance for each PC coefficient on each exposure

    for(int ipc=0;ipc<nvar;++ipc) {
      DVector ave_exp(exps.size());
      for(int iexp=0;iexp<exps.size();++iexp) {
	DVector ave(nused_ccd);      
	int exp_ndet=0;


	// store info of all detections on this exposure
	for(int iccd=0;iccd<nused_ccd;++iccd) exp_ndet+=exp_ccd[iexp][iccd].size();
	DVector exps_det(exp_ndet);


	int cur_det=0;
	for(int iccd=0;iccd<nused_ccd;++iccd) {
	  double sum_ccd=0;
	  
	  int ndet=exp_ccd[iexp][iccd].size();
	  for(int idet=0;idet<ndet;++idet) {

	    double pc=U(exp_ccd[iexp][iccd][idet],ipc);
	    sum_ccd+=pc;
	    
	    exps_det(cur_det)=pc;
	    cur_det++;
	  }
	  ave(iccd)=sum_ccd/ndet;

	  //cout<<"PC "<<ipc<<" ccd "<<iccd<<" "<<ave(iccd)<<" "<<det_var<<endl;
	}
	double var=variance(ave);
	double det_var=variance(exps_det);
	cout<<"PC "<<ipc<<" "<<var<<" "<<det_var/sqrt(exp_ndet)<<endl;
	ave_exp(iexp)=var;
      }

      double vmean=ave_exp.sumElements()/exps.size();
      double var=variance(ave_exp);
      //cout<<"PC "<<ipc<<" "<<vmean<<" +- "<<sqrt(var)<<endl;
    }
    
    
    
    ofstream pcdata(pcname.c_str());
    if(shapelet) {
      int npix=0;
      CDVector rpos(nx*ny),pos(nx*ny);
      for(int ix=0;ix<nx;++ix) {
	for(int iy=0;iy<ny;++iy) {
	  double x=rmax/nx*(ix+0.5);
	  double y=rmax/ny*(iy+0.5);
	  complex<double> z1(x,y);
	  pos(npix)=z1;
	  x=x-rmax/2;
	  y=y-rmax/2;
	  complex<double> z(x,y);
	  rpos(npix)=z;
	  npix++;
	}
      }
      
      

      int shape_var=(shape_order+1)*(shape_order+2)/2;
      int start=0;
      int end=shape_var-3;
      if(add_size) {
	start=1;
	end++;
      }
      for(int ipc=0;ipc<Vt.nrows();++ipc) {
	DVector b(shape_var);
	
	// always set first three to default
	b(0)=1;
	b(1)=0;
	b(2)=0;
	
	b.subVector(3,shape_var)=Vt.row(ipc).subVector(start,end);
	
	DMatrix psi(npix,shape_var);
	MakePsi(psi, rpos.view(), shape_order);
	DVector data=psi*b;
	for(int ipix=0;ipix<npix;++ipix) {
	  pcdata<<ipc<<" "<<pos(ipix).real()<<" "<<pos(ipix).imag()<<" "<<data(ipix)<<endl;
	}   
      }
    }
  
  
  if(fitfile) delete fitfile;
}

