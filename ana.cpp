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

int main(int argc,char*argv[])
{


  // use for now to get focal plane coor.
  map<int,vector<double> > focal;
  focal[1].push_back(-201.348000);
  focal[1].push_back(-94.610000);
  focal[2].push_back(-201.348000);
  focal[2].push_back(-30.720000);
  focal[3].push_back(-201.348000);
  focal[3].push_back(33.170000);
  focal[4].push_back(-167.532000);
  focal[4].push_back(-126.555000);
  focal[5].push_back(-167.532000);
  focal[5].push_back(-62.665000);
  focal[6].push_back(-167.532000);
  focal[6].push_back(1.225000);
  focal[7].push_back(-167.532000);
  focal[7].push_back(65.115000);
  focal[8].push_back(-133.716000);
  focal[8].push_back(-158.500000);
  focal[9].push_back(-133.716000);
  focal[9].push_back(-94.610000);
  focal[10].push_back(-133.716000);
  focal[10].push_back(-30.720000);
  focal[11].push_back(-133.716000);
  focal[11].push_back(33.170000);
  focal[12].push_back(-133.716000);
  focal[12].push_back(97.060000);
  focal[13].push_back(-99.900000);
  focal[13].push_back(-190.445000);
  focal[14].push_back(-99.900000);
  focal[14].push_back(-126.555000);
  focal[15].push_back(-99.900000);
  focal[15].push_back(-62.665000);
  focal[16].push_back(-99.900000);
  focal[16].push_back(1.225000);
  focal[17].push_back(-99.900000);
  focal[17].push_back(65.115000);
  focal[18].push_back(-99.900000);
  focal[18].push_back(129.005000);
  focal[19].push_back(-66.084000);
  focal[19].push_back(-190.445000);
  focal[20].push_back(-66.084000);
  focal[20].push_back(-126.555000);
  focal[21].push_back(-66.084000);
  focal[21].push_back(-62.665000);
  focal[22].push_back(-66.084000);
  focal[22].push_back(1.225000);
  focal[23].push_back(-66.084000);
  focal[23].push_back(65.115000);
  focal[24].push_back(-66.084000);
  focal[24].push_back(129.005000);
  focal[25].push_back(-32.268000);
  focal[25].push_back(-222.390000);
  focal[26].push_back(-32.268000);
  focal[26].push_back(-158.500000);
  focal[27].push_back(-32.268000);
  focal[27].push_back(-94.610000);
  focal[28].push_back(-32.268000);
  focal[28].push_back(-30.720000);
  focal[29].push_back(-32.268000);
  focal[29].push_back(33.170000);
  focal[30].push_back(-32.268000);
  focal[30].push_back(97.060000);
  focal[31].push_back(-32.268000);
  focal[31].push_back(160.950000);
  focal[32].push_back(1.548000);
  focal[32].push_back(-222.390000);
  focal[33].push_back(1.548000);
  focal[33].push_back(-158.500000);
  focal[34].push_back(1.548000);
  focal[34].push_back(-94.610000);
  focal[35].push_back(1.548000);
  focal[35].push_back(-30.720000);
  focal[36].push_back(1.548000);
  focal[36].push_back(33.170000);
  focal[37].push_back(1.548000);
  focal[37].push_back(97.060000);
  focal[38].push_back(1.548000);
  focal[38].push_back(160.950000);
  focal[39].push_back(35.364000);
  focal[39].push_back(-190.445000);
  focal[40].push_back(35.364000);
  focal[40].push_back(-126.555000);
  focal[41].push_back(35.364000);
  focal[41].push_back(-62.665000);
  focal[42].push_back(35.364000);
  focal[42].push_back(1.225000);
  focal[43].push_back(35.364000);
  focal[43].push_back(65.115000);
  focal[44].push_back(35.364000);
  focal[44].push_back(129.005000);
  focal[45].push_back(69.180000);
  focal[45].push_back(-190.445000);
  focal[46].push_back(69.180000);
  focal[46].push_back(-126.555000);
  focal[47].push_back(69.180000);
  focal[47].push_back(-62.665000);
  focal[48].push_back(69.180000);
  focal[48].push_back(1.225000);
  focal[49].push_back(69.180000);
  focal[49].push_back(65.115000);
  focal[50].push_back(69.180000);
  focal[50].push_back(129.005000);
  focal[51].push_back(102.996000);
  focal[51].push_back(-158.500000);
  focal[52].push_back(102.996000);
  focal[52].push_back(-94.610000);
  focal[53].push_back(102.996000);
  focal[53].push_back(-30.720000);
  focal[54].push_back(102.996000);
  focal[54].push_back(33.170000);
  focal[55].push_back(102.996000);
  focal[55].push_back(97.060000);
  focal[56].push_back(136.812000);
  focal[56].push_back(-126.555000);
  focal[57].push_back(136.812000);
  focal[57].push_back(-62.665000);
  focal[58].push_back(136.812000);
  focal[58].push_back(1.225000);
  focal[59].push_back(136.812000);
  focal[59].push_back(65.115000);
  focal[60].push_back(170.628000);
  focal[60].push_back(-94.610000);
  focal[61].push_back(170.628000);
  focal[61].push_back(-30.720000);
  focal[62].push_back(170.628000);
  focal[62].push_back(33.170000);

  double xsize=2048*15e-6*1000;
  double ysize=4096*15e-6*1000;


  
  ConfigFile params;
  params.setDelimiter("=");
  params.setInclude("+");
  params.setComment("#");
  params.load(argv[1]);
  for(int k=2;k<argc;k++) params.append(argv[k]);

  std::string filename= params.read<std::string>("file","");

  std::string pcname= params.read<std::string>("pcfile");
  int pc= params.read<int>("pc",-1);
  std::string dir= params.read<std::string>("dir","./");
  std::string image_dir= params.read<std::string>("image_dir","");
  int max_exp= params.read<int>("max_exp",-1);
  std::string outname= params.read<std::string>("outname","out");
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
  bool read2=params.read<bool>("read2",false);
  float readmax=params.read<float>("readmax",1);
  string use_dir=params.read<string>("use_dir",".");
  string cdir=params.read<string>("cdir",".");
  string out_dir=params.read<string>("out_dir",".");
  int threads=params.read<int>("threads",1);
  int seed=params.read<int>("seed",11111);

  FILELog::ReportingLevel() = FILELog::FromInt(logging);
  FILE_LOG(logINFO)<<"Settings...\n"<<params<<endl;


  // read in pc file
  CCfits::FITS *pcfits=new FITS(pcname, CCfits::Read);
    
  CCfits::ExtHDU& table=pcfits->extension("exps");
  long nexp=table.rows();
  
  vector<string> pc_exposures;
  table.column("exposure").read(pc_exposures, 1, nexp);

  int ccd,nvar,nx,ny;
  bool add_size;
  double xmax,ymax;
  table.readKey("ccd",ccd);
  table.readKey("nvar",nvar);
  table.readKey("nx",nx);
  table.readKey("ny",ny);
  table.readKey("xmax",xmax);
  table.readKey("ymax",ymax);
  table.readKey("add_size",add_size);

  // I didn't always write this out, but I 
  // have never used anything different
  int shapestart=3;

  //  table.readKey("add_size",skip61);
  // assume it skipped for now
  bool skip61=true;
  //if(skip61) ccd-=1;
    
  int ntot=nx*ny*nvar*ccd;

  // read in the reconstruction
  DMatrix coeff, vec,data,dataR;
  DDiagMatrix S;
  DVector mean,Sdiag;

  // read in the singluar vector to get the number of PCs in the file
  readVecFromFits(pcfits,"singular",Sdiag);
  int file_pc=Sdiag.size();
  if(pc<0) pc=file_pc;

  readMatrixFromFits(pcfits,"dataR",data);
  // we don't need to reconstruct the matrix just read it in
  if(file_pc==pc) {
    readMatrixFromFits(pcfits,"dataR",dataR);
  }
  else {
    readVecFromFits(pcfits,"mean",mean);
    readMatrixFromFits(pcfits,"coeff",coeff);
    readMatrixFromFits(pcfits,"vec",vec);

    Assert(nexp==coeff.nrows());
    Assert(ntot==vec.ncols());
    S.resize(Sdiag.size());
    S.diag()=Sdiag;
    dataR.resize(nexp,ntot);
    dataR=coeff.subMatrix(0,nexp,0,pc)*S.subDiagMatrix(0,pc)*vec.subMatrix(0,pc,0,ntot);

    // add in mean
    for(int icol=0;icol<dataR.ncols();++icol) {
      dataR.col(icol).addToAll(mean(icol));
    }
  }


  // Read in shapelet file
  ifstream file(filename.c_str());  
  string name;
  std::vector<Exposure<double> > exps;

  // map exposure index to index in pc file
  map<int,int> exp_map;
  int iexp=0;
  while(file>>name) {

    // assume shapestart
    // shapestart=3

    // check that the exposure exists in the file and find the index
    bool found=false;
    for(int i=0;i<nexp;++i) {
      if(pc_exposures[i]==name) {
	exp_map[iexp]=i;
	found=true;
	break;
      }
    }

    if(!found) {
      cout<<"Cannot find the exposure name "<<name<<" skipping "<<endl;
      continue;
    }
    
    Exposure<double> exp(name,ccd,shapestart);
    exp.setChipDivide(1,1);// put all objects into one cell
    exp.setChipMax(xmax,ymax);

    if(skip61) exp.addSkip(61);
    bool suc;
    suc=exp.readShapelet(dir+name+"/",nvar,add_size,
			 do_em,use_dash,suffix, prefix+name,readmax
			 ,use_dir,cdir);

    if(suc) exps.push_back(exp);
    if(exps.size()>(max_exp-1) && max_exp>0) break;
    iexp++;
  }

  // Now iterate through all the exposures/ccds


  int current_star=0;
#pragma omp parallel for shared(exps,dataR)
  for(int iexp=0;iexp<exps.size();++iexp) {


    std::vector<double> ie1_vec,ie2_vec,isize_vec,e1_diff,
      e2_diff,size_diff,x_vec,y_vec,xf_vec,yf_vec,ra_vec,dec_vec;
    std::vector<int> ccd_vec;
    // go through the chips
    std::map<int,Chip<double>*>::const_iterator iter=exps[iexp].chips.begin();
    int ichip=0;
    for(; iter!=exps[iexp].chips.end();++iter,++ichip) {

      int use_chip;
      if(ichip<60) use_chip=ichip+1;
      else use_chip=ichip+2;

      // There is only one cell here
      for(int istar=0;istar < iter->second->getCell(0)->getNDet();++istar) {
	DVector star(nvar);
	double xstar,ystar;
	Position<float> pos=iter->second->getCell(0)->getDet(istar)->getPos();
	Position<double> sky=iter->second->getCell(0)->getDet(istar)->getSky();
	for(int ivar=0;ivar<nvar;++ivar) {
	  
	  star(ivar)=iter->second->getCell(0)->getDet(istar)->getVal(ivar);
	}
	
	double x=pos.x;
	double y=pos.y;
	double ra=sky.x;
	double dec=sky.y;
	// rough focal plane for now, update with WCS coords later
	// could read it into exposure class
	double xf=x*15e-6*1000+focal[use_chip][0];
	double yf=y*15e-6*1000+focal[use_chip][1];


	int e1_index,e2_index,size_index;
	if(add_size) e1_index=1;
	else e1_index=0;
	if(add_size) e2_index=2;
	else e2_index=1;
	if(add_size) size_index=3;
	else size_index=2;


	double e1=sqrt(2)*star(e1_index)/(1+star(size_index));
	double e2=sqrt(2)*star(e2_index)/(1+star(size_index));
	double size=-1;
	if(add_size) {
	  size=2.35*(sqrt(star(0)*star(0)*(1+star(size_index))));
	}

	// Now we need the reconstructed values
	// find the bin that this belongs to
	// add interpolation smoothing here across ccd
	int binx=static_cast<int>(x/(xmax/nx));
	int biny=static_cast<int>(y/(ymax/ny));
	int bin=binx*ny+biny;
	int icell=(ichip)*nx*ny*nvar+bin*nvar;
	DVector istar=dataR.subVector(exp_map[iexp],icell,0,1,nvar);

	double ie1=sqrt(2)*istar(e1_index)/(1+istar(size_index));
	double ie2=sqrt(2)*istar(e2_index)/(1+istar(size_index));
	double isize=-1;
	if(add_size) {
	  isize=2.35*(sqrt(istar(0)*istar(0)*(1+istar(size_index))));
	}
	  
	//cout<<x<<" "<<y<<" "<<e1<<" "<<ie1<<" "<<ra<<" "<<dec<<endl;
	//cout<<x<<" "<<y<<" "<<binx<<" "<<biny<<" "<<icell<<" "<<star(0)<<" "<<istar(0)<<" "<<star(1)<<" "<<istar(1)<<endl;
	ie1_vec.push_back(ie1);
	ie2_vec.push_back(ie2);
	isize_vec.push_back(isize);
	e1_diff.push_back(ie1-e1);
	e2_diff.push_back(ie2-e2);
	size_diff.push_back(isize-size);
	x_vec.push_back(x);
	y_vec.push_back(y);
	xf_vec.push_back(xf);
	yf_vec.push_back(yf);
	ra_vec.push_back(ra);
	dec_vec.push_back(dec);
	ccd_vec.push_back(use_chip);
	current_star++;
      }
    }

    // Write out a fits file with residual info in it
    FITS *rfitfile=0;
    long naxis    =   0;      
    long naxes1[2] = { 0, 0 }; 
    rfitfile=new FITS("!r_"+outname+"_"+exps[iexp].getLabel()+".fits",USHORT_IMG , naxis , naxes1 );
    // write the exposure information
    int nwvar=13;
    int nrows=current_star;
    std::vector<string> colName(nwvar,"");
    std::vector<string> colForm(nwvar,"");
    std::vector<string> colUnit(nwvar,"");
    colName[0] = "ccd";
    colName[1] = "x";
    colName[2] = "y";
    colName[3] = "xf";
    colName[4] = "yf";
    colName[5] = "e1_model";
    colName[6] = "e2_model";
    colName[7] = "size_model";
    colName[8] = "e1_diff";
    colName[9] = "e2_diff";
    colName[10]= "size_diff";
    colName[11] = "ra";
    colName[12] = "dec";
        
    colForm[0] = "1J";
    colForm[1] = "1E";
    colForm[2] = "1E";
    colForm[3] = "1E";
    colForm[4] = "1E";
    colForm[5] = "1E";
    colForm[6] = "1E";
    colForm[7] = "1E";
    colForm[8] = "1E";
    colForm[9] = "1E";
    colForm[10]= "1E";
    colForm[11] = "1E";
    colForm[12]= "1E";
    
      
    colUnit[0] = "";
    colUnit[1] = "";
    colUnit[2] = "";
    colUnit[3] = "";
    colUnit[4] = "";
    colUnit[5] = "";
    colUnit[6] = "";
    colUnit[7] = "";
    colUnit[8] = "";
    colUnit[9] = "";
    colUnit[10] = "";
    colUnit[11] = "";
    colUnit[12] = "";
    Table* newTable = rfitfile->addTable("residuals",nrows,colName,colForm,colUnit);
    newTable->column(colName[0]).write(ccd_vec,1);
    newTable->column(colName[1]).write(x_vec,1);
    newTable->column(colName[2]).write(y_vec,1);
    newTable->column(colName[3]).write(xf_vec,1);
    newTable->column(colName[4]).write(yf_vec,1);
    newTable->column(colName[5]).write(ie1_vec,1);
    newTable->column(colName[6]).write(ie2_vec,1);
    newTable->column(colName[7]).write(isize_vec,1);
    newTable->column(colName[8]).write(e1_diff,1);
    newTable->column(colName[9]).write(e2_diff,1);
    newTable->column(colName[10]).write(size_diff,1);
    newTable->column(colName[11]).write(ra_vec,1);
    newTable->column(colName[12]).write(dec_vec,1);
    
  }


}

  

