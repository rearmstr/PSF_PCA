#include "Data.h"
#include "dbg.h"
#include "Utils.h"
#include <CCfits/CCfits>
#include "Random.h"
using namespace CCfits;
using namespace PCA;


// calculate cell values
void Cell::calcVals(std::string type,const std::vector<float> &params) 
{
  
  if(type=="MeanClipped") {
    
    Assert(params.size()>=1);
    double sigma=params[0];
    for(int ivar=0;ivar<nvar;++ivar) {
      
      xdbg<<"Getting variable: "<<ivar<<std::endl;
      std::vector<float> tmp;     
      
      int cur=0;
      for(int idet=0;idet<dets.size();++idet) {
	if(dets[idet]->isClipped()) continue;
	xdbg<<"From det: "<<idet<<" "<<(*dets[idet])(ivar)<<std::endl;
	tmp.push_back((*dets[idet])(ivar));
	cur++;
      }
      
      // if number of objects too small mad will not work 
      // just return the mean
      if(cur<3) {
	xdbg<<"Too few objects to do robust estimate: "
	    <<cur<<" using mean"<<std::endl;
	float sum=0;
	for(int i=0;i<cur;++i) sum+=tmp[i];
	data(ivar)=sum/cur;
      }
      else {

	float med,mean,mad;
	std::vector<int> clipped;
	med=median_mad_flag(tmp,mad,mean,sigma,clipped);
	data(ivar)=mean;
      }
    }
  }
  else if(type=="Mean") {

    for(int ivar=0;ivar<nvar;++ivar) {
      
      float sum=0;
      int cur=0;
      for(int idet=0;idet<dets.size();++idet) {
	if(dets[idet]->isClipped()) continue;
	sum+=(*dets[idet])(ivar);
	cur++;
      } 
      data(ivar)=sum/cur;
    }
  }
}



bool BinnedData::addExp(Exposure<float> *exp)
{
  int cur_size=exp_list.size();
  exp_list.push_back(exp);
  exp_map[exp->getLabel()]=cur_size;

  if(exp->getNChip()!=nchips) {
    std::cout<<"Skipping exposure "<<exp->getLabel()
	     <<" not the correct number of ccds"<<std::endl;
    return false;
  }
    
  std::vector<std::vector<Cell> > exp_cells;
  
  int nchip=0;
  Exposure<float>::ChipIter chip_iter=exp->ChipBegin();
  for(; chip_iter!=exp->ChipEnd(); chip_iter++,nchip++) {

    // assume area of each ccd is the same
    if(chip_iter==exp->ChipBegin()) {
      cbounds=chip_iter->second->getBounds().divide(nx,ny);
    }
    // store boundaries for future use

    std::vector<Cell> ccd_cells(nx*ny,Cell(nvar));
    float xmax=chip_iter->second->getBounds().getXMax();
    float ymax=chip_iter->second->getBounds().getYMax();
    float xmin=chip_iter->second->getBounds().getXMin();
    float ymin=chip_iter->second->getBounds().getYMin();
    for(int i=0;i<cbounds.size();++i) {
      xdbg<<"creating bounds for cell "<<i<<" on chip "
	  <<chip_iter->second->getLabel()<<" "<<cbounds[i]<<std::endl;
      ccd_cells[i].setBounds(cbounds[i]);
    }
    
    Chip<float>::DetIter det_iter=chip_iter->second->DetBegin();
    int idet=0;
    for(; det_iter!=chip_iter->second->DetEnd(); ++det_iter,++idet) {
      
      if(!chip_iter->second->getBounds().includes((*det_iter)->getPos())) {
	xdbg<<"invalid object at "<<(*det_iter)->getPos()
	    <<" skipping object "<<std::endl;
	continue;
      }
	
      int bin_x=static_cast<int>(((*det_iter)->getX()-xmin)/((xmax-xmin)/nx));
      int bin_y=static_cast<int>(((*det_iter)->getY()-ymin)/((ymax-ymin)/ny));
      int bin=bin_x*ny+bin_y;


      Assert(bin>=0 && bin<cbounds.size());
      ccd_cells[bin].addDet(*det_iter);

      xdbg<<"adding det "<<idet<<" at "<<(*det_iter)->getPos()	
	  <<" from chip "
	  <<chip_iter->second->getLabel()<<" to bin "
	  <<bin<<std::endl;
    }
    
    exp_cells.push_back(ccd_cells);
  }
  
  cells.push_back(exp_cells);
  return true;
}

void BinnedData::calcData(std::string type,const std::vector<float> &param) 
{
  cur_exp=0;

  if(!initialized) {
    data.resize(exp_list.size(),nchips*nx*ny*nvar);
    for(int i=0;i<exp_list.size();++i) {
      std::vector<bool> tmp(nchips*nx*ny*nvar,false);
      missing.push_back(tmp);
    }
    data.setAllTo(defaultVal);
    initialized=true;
  }
  
  for(int iexp=0;iexp<exp_list.size();++iexp) {

    if(exp_list[iexp]->isOutlier()) continue;

    for(int iccd=0;iccd<cells[iexp].size();++iccd) {

      for(int icell=0;icell<cells[iexp][iccd].size();++icell) {

	int firstvar=indx.getFirstVar(iccd,icell);
	int lastvar=firstvar+nvar;
	
	// Check if cell is missing or clipped.  If so reset
	// data values to defaultVal
	if(cells[iexp][iccd][icell].isMissing() ||
	   cells[iexp][iccd][icell].isClipped()) {
	
	  for(int imiss=firstvar;imiss<firstvar;++imiss) {
	    missing[iexp][imiss]=true;
	  }
	 

	  xdbg<<"setting "<<iccd<<" "<<icell
	      <<" to missing and defualt"<<std::endl;
	  data.subMatrix(cur_exp,cur_exp+1,
			 firstvar,
			 lastvar).setAllTo(defaultVal);
	  continue;
	}

	
	cells[iexp][iccd][icell].calcVals(type,param);	
	
	xdbg<<"setting "<<iccd<<" "<<icell
	    <<" data to "<<cells[iexp][iccd][icell].getData()<<std::endl;
	
	data.subVector(cur_exp,firstvar,0,1,
		       lastvar-firstvar)=cells[iexp][iccd][icell].getData();
      }
      
    }
    
    cur_exp++;
    
  }

  mean.resize(data.ncols());
  meanRemove(data,mean,missing);
  
}



void PcaMethod::solve() {

  FMatrixView mdata=data->getData();
  int nexp=mdata.nrows();
  int nvar=mdata.ncols();
  if(nexp > nvar) {
    Vt.resize(nvar,nvar);
    S.resize(nvar);
    U.resize(nexp,nvar);
    U=mdata;
    SV_Decompose(U,S,Vt,true);
      
  }
  else {
    
    Vt.resize(nexp,nvar);
    S.resize(nexp);
    U.resize(nexp,nexp);
    Vt = mdata;
    SV_Decompose(Vt.transpose(),S,U.transpose(),true);
    
  }
}

void PcaMethod::writeFits(string file)
{
  
  try {
    FITS *fitfile=0;
    long naxis    =   0;      
    long naxes1[2] = { 1, 1 }; 
    fitfile=new FITS("!"+file,USHORT_IMG , naxis , naxes1 );

    // Write exposure names to output file
    // make two lists, those that were rejected and those that
    // succeeded.

    std::vector<std::string> exps;
    std::vector<std::string> miss_exps;

    int nexp=data->getNExp();
    for(int i=0;i<nexp; ++i) {
      Exposure<> *exp=data->getExp(i);
      if(exp->isOutlier()) miss_exps.push_back(exp->getLabel());
      else exps.push_back(exp->getLabel());
    }
   
    int nwvar=1;
    int nrows=exps.size();
    std::vector<string> colName(nwvar,"");
    std::vector<string> colForm(nwvar,"");
    std::vector<string> colUnit(nwvar,"");
    colName[0] = "exposure";
    colForm[0] = "16A";
    colUnit[0] = "";
    Table* newTable = fitfile->addTable("exps",nrows,colName,colForm,colUnit);
    newTable->column(colName[0]).write(exps,1);

    newTable->addKey("nx",data->getNx(),"bins in x");
    newTable->addKey("ny",data->getNy(),"bins in x");
    newTable->addKey("nccd",data->getNCCD(),"number of ccds");
    newTable->addKey("nvar",data->getNVar(),"number of variables per cell");

    // include a table of failed exposures
    nrows=miss_exps.size();
    if(nrows>0) {
      colName[0] = "missing_exposure";
    
      Table* newTable2 = fitfile->addTable("miss_exps",nrows,colName,
					   colForm,colUnit);
      newTable2->column(colName[0]).write(miss_exps,1);
    }
  

    // include table of grid positions
    std::vector<double> lx,ux,ly,uy;
    std::vector<Bounds<float> > cb=data->getCellBounds();
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
    Table* newTable3 = fitfile->addTable("grid",nrows,colName,colForm,colUnit);
    newTable3->column(colName[0]).write(lx,1);
    newTable3->column(colName[1]).write(ly,1);
    newTable3->column(colName[2]).write(ux,1);
    newTable3->column(colName[3]).write(uy,1);


    // need to write info from Binned Data
    writeVectorToFits(fitfile,data->getMean(),"mean");
    writeMatrixToFits(fitfile,U,"U");
    writeMatrixToFits(fitfile,Vt,"Vt");
    writeVectorToFits(fitfile,S.diag(),"S");
  }
  catch (CCfits::FitsException& e) {
    std::cout<<e.message()<<std::endl;
    
  }
}






void EMPcaMethod::solve() {

  FMatrixView mdata=data->getData();
  int nexp=mdata.nrows();
  int nvar=mdata.ncols();
  if(nexp > nvar) {
   
    if(npc>nvar) {
      std::cerr<<"Not enough variables for npc"<<std::endl;
      exit(1);
    }
    C.resize(nexp,npc);
    x.resize(npc,nvar);
  }
  else {
    if(npc>nexp) {
      std::cerr<<"Not enough exposures for npc"<<std::endl;
      exit(1);
    }
    
    C.resize(nvar,npc);
    x.resize(npc,nexp);
  }
  
  C.setZero();
  x.setZero();
  
  // for now fill with random values
  // COULD fill with something else
  ran::GaussianDeviate gauss;
  for(int i=0;i<C.nrows();++i) {
    for(int j=0;j<C.ncols();++j) {
      C(i,j)=gauss()*sigma;
    }
  }
  
  for(int i=0;i<x.nrows();++i) {
    for(int j=0;j<x.ncols();++j) {
      x(i,j)=gauss()*sigma;
    }
  }
  
  bool do_missing=false;
  double prev_diff=1e10;
  for(int iter=0;iter<max_iter;++iter) {
    
    if(!do_missing) {
      
      FMatrix tmp=C.transpose()*C;
      x=C.transpose()*mdata.transpose()/tmp;
      
      FMatrix Cnew=mdata.transpose()*x.transpose()%(x*x.transpose());

      double diff=((Cnew*x).transpose()-mdata).norm();
      diff/=(Cnew.nrows()*Cnew.ncols());

      double rel_diff=std::abs(diff-prev_diff)/prev_diff;
      C=Cnew;
      if(iter>min_iter && 
	 (diff<tol || rel_diff<tol) ) break;
      prev_diff=diff;
      
    }
    
    
  }
  
}


void EMPcaMethod::writeFits(string file)
{
  
  try {
    FITS *fitfile=0;
    long naxis    =   0;      
    long naxes1[2] = { 1, 1 }; 
    fitfile=new FITS("!"+file,USHORT_IMG , naxis , naxes1 );

    // Write exposure names to output file
    // make two lists, those that were rejected and those that
    // succeeded.

    std::vector<std::string> exps;
    std::vector<std::string> miss_exps;

    int nexp=data->getNExp();
    for(int i=0;i<nexp; ++i) {
      Exposure<> *exp=data->getExp(i);
      if(exp->isOutlier()) miss_exps.push_back(exp->getLabel());
      else exps.push_back(exp->getLabel());
    }
   
    int nwvar=1;
    int nrows=exps.size();
    std::vector<string> colName(nwvar,"");
    std::vector<string> colForm(nwvar,"");
    std::vector<string> colUnit(nwvar,"");
    colName[0] = "exposure";
    colForm[0] = "16A";
    colUnit[0] = "";
    Table* newTable = fitfile->addTable("exps",nrows,colName,colForm,colUnit);
    newTable->column(colName[0]).write(exps,1);

    newTable->addKey("nx",data->getNx(),"bins in x");
    newTable->addKey("ny",data->getNy(),"bins in x");
    newTable->addKey("nccd",data->getNCCD(),"number of ccds");
    newTable->addKey("nvar",data->getNVar(),"number of variables per cell");

    // include a table of failed exposures
    nrows=miss_exps.size();
    if(nrows>0) {
      colName[0] = "missing_exposure";
    
      Table* newTable2 = fitfile->addTable("miss_exps",nrows,colName,
					   colForm,colUnit);
      newTable2->column(colName[0]).write(miss_exps,1);
    }
  

    // include table of grid positions
    std::vector<double> lx,ux,ly,uy;
    std::vector<Bounds<float> > cb=data->getCellBounds();
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
    Table* newTable3 = fitfile->addTable("grid",nrows,colName,colForm,colUnit);
    newTable3->column(colName[0]).write(lx,1);
    newTable3->column(colName[1]).write(ly,1);
    newTable3->column(colName[2]).write(ux,1);
    newTable3->column(colName[3]).write(uy,1);


    // need to write info from Binned Data
    writeVectorToFits(fitfile,data->getMean(),"mean");
    writeMatrixToFits(fitfile,C,"C");
    writeMatrixToFits(fitfile,x,"x");

  }
  catch (CCfits::FitsException& e) {
    std::cout<<e.message()<<std::endl;
    
  }
}
