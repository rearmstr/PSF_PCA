#include "Detection.h"
#include <CCfits/CCfits>
using namespace PCA;
using namespace CCfits;




template<class T>
Exposure<T>::Exposure(std::string _label,int _nvar,
		      float _x_min, float _x_max,
		      float _y_min, float _y_max,
		      std::vector<int> _skip,
		      int _shapeStart,
		      bool _addSize,
		      bool _includeMiss):
  label(_label), bounds(_x_min,_x_max,_y_min,_y_max),
  skip(_skip),
  shapeStart(_shapeStart),addSize(_addSize),
  includeMiss(_includeMiss),nvar(_nvar),outlier(false)
{
  
}




template<class T>
Exposure<T>::Exposure(std::string _label,int _nvar,
		      Bounds<float> _bounds,
		      std::vector<int> _skip,
		      int _shapeStart,
		      bool _addSize,
		      bool _includeMiss):
  label(_label), bounds(_bounds),
  skip(_skip),
  shapeStart(_shapeStart),addSize(_addSize),
  includeMiss(_includeMiss),nvar(_nvar),outlier(false)
{
  
}


template<class T>
bool Exposure<T>::readShapeletDir(std::string dir,int nchip,
				  std::string suffix,std::string exp)
{
  if (exp.empty()) exp=label;
  dbg<< "Reading exposure " << exp<<std::endl;
  for(int ichip=1;ichip<=nchip;++ichip) {
    max_ccd=ichip;
    std::vector<Detection<T> *> dets;
    // check if this chip should be skipped
    std::vector<int>::iterator iter=find(skip.begin(),skip.end(),ichip);
    if(iter!=skip.end()) continue;
    
    std::stringstream inputFile;
    inputFile << dir << "/" << exp << "_";
    
    if(ichip<10) inputFile <<0;
    inputFile << ichip << "_"+suffix;
    
    
    try {
      dbg << "opening file " << inputFile.str()<<std::endl;
      std::auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(inputFile.str(),CCfits::Read));
      
      CCfits::ExtHDU& table = pInfile->extension(1);
      
      long nTabRows=table.rows();
      dbg << "found " << nTabRows<<" objects"<<std::endl;
      long start=1;
      long end=nTabRows;
      
      std::vector<int> psf_flags;
      std::vector<double> psf_size;
      std::vector<double> xpos;
      std::vector<double> ypos;
      
      table.column("psf_flags").read(psf_flags, start, end);
      // only need the first row since they are all the same
      table.column("sigma_p").read(psf_size, start, start+1);
      table.column("x").read(xpos, start, end);
      table.column("y").read(ypos, start, end);
      
      std::vector<long> order;       // shapelet order
      table.column("psf_order").read(order, start, end);
      
      
      for (int i=0; i<nTabRows; i++) {
	if (!psf_flags[i]) {          // pass psf flags
	  
	  int row=i+1;
	  Detection<T> *det=new Detection<T>(xpos[i],ypos[i],nvar);
	  int ncoeff=(order[i]+1)*(order[i]+2)/2;      // psf values
	  std::valarray<double> coeffs;
	  table.column("shapelets").read(coeffs, row); 
	  
	  xdbg<<"adding object "<<ichip<<" "<<xpos[i]<<" "
	      <<ypos[i]<<" "<<nvar<<std::endl;

	  
	  int last_index=nvar;
	  int start_index=0;
	  if(addSize) {
	    xxdbg<<"adding size "<<psf_size[0]<<std::endl;
	    (*det)(0)=psf_size[0];
	    last_index--;
	    start_index++;
	  }
	  for(int j=0;j<last_index;++j) {
	      xxdbg<<"adding index var "<<start_index+j<<" "
		     <<"from shapelet index "<<shapeStart+j<<" value:"
		     <<coeffs[shapeStart+j]<<std::endl;
	      
	      (*det)(start_index+j)=coeffs[shapeStart+j];
	  }

	  dets.push_back(det);
	}
      }
      
    }
    
    catch (CCfits::FitsException& ) {
      
      if(!includeMiss) {
	dbg<<"Can't open chip: "<<inputFile.str()
	   <<" from exposure "<<exp<<" skipping"<<std::endl;
	return false;
      }
      else {
	dbg<<"Can't open chip: "<<inputFile.str()
	   <<" from exposure "<<exp<<" will be set to missing"<<std::endl;
      }
    }
    
    if(chips.find(ichip)==chips.end()) {
      xdbg<<"Could not find chip "<<ichip<<" creating new one"<<std::endl;
      Chip<T> *chip=new Chip<T>(ichip,bounds);
      addChip(ichip,chip);
    }
    for(int i=0;i<dets.size();++i) {
      chips[ichip]->addDet(dets[i]);
    }
    
  }

  return true;
}

template<class T>
void Exposure<T>::writeFits(std::string dir,std::string suffix)
{

  // write the exposure information
  int nrows=0;
  dbg<<"filling arrays"<<std::endl;

  std::vector<T> x,y;
  std::vector<int> ccd;
  std::vector<bool> clip;
  std::vector<valarray<T> > vals;

  ChipIter c_it=ChipBegin();
  for( ; c_it!=ChipEnd();++c_it) {

    nrows+=c_it->second->getNDet();
    typename Chip<T>::DetIter d_it=c_it->second->DetBegin();

    for( ; d_it!=c_it->second->DetEnd();++d_it) {
      x.push_back((*d_it)->getX());
      y.push_back((*d_it)->getY());
      ccd.push_back(c_it->second->getLabel());
      clip.push_back((*d_it)->isClipped());
      vals.push_back((*d_it)->getValarray());
    }
  }


  int nwvar=5;

  std::vector<string> colName(nwvar,"");
  std::vector<string> colForm(nwvar,"");
  std::vector<string> colUnit(nwvar,"");
  colName[0] = "x";
  colName[1] = "y";
  colName[2] = "ccd";
  colName[3] = "clip";
  colName[4] = "vals";

  std::stringstream v_form;
  v_form << nvar << "D";
  
  colForm[0] = "1E";
  colForm[1] = "1E";
  colForm[2] = "1J";
  colForm[3] = "1I";
  colForm[4] = v_form.str();

  
  colUnit[0] = "";
  colUnit[1] = "";
  colUnit[2] = "";
  colUnit[3] = "";
  colUnit[4] = "";

  
  string name=dir+"/"+label+suffix;
  long naxis    =   2;      
  long naxes1[2] = { 1, 1 }; 
  FITS fitfile("!"+name,CCfits::Write);
  
  dbg<<"creating file"<<ChipEnd()->first<<std::endl;  
  Table* newTable = fitfile.addTable("data",nrows,colName,colForm,colUnit);
 
  dbg<<"adding keys"<<" "<<std::endl;
  newTable->addKey("name",label,"name or exposure");
  newTable->addKey("addSize",addSize,"size included in vector");
  newTable->addKey("sStart",shapeStart,"start of shapelet vector");
  newTable->addKey("nvar",nvar,"number of variables");
  newTable->addKey("nccd",max_ccd,"number of CCDs"); 
  newTable->addKey("x_max",bounds.getXMax(),"x max"); 
  newTable->addKey("y_max",bounds.getYMax(),"y max"); 
  newTable->addKey("x_min",bounds.getXMin(),"x min"); 
  newTable->addKey("y_min",bounds.getYMin(),"y min"); 
  newTable->addKey("incMiss",includeMiss,"include missing CCDs"); 

  std::stringstream skipString;
  for(int i=0;i<skip.size();++i) {
    skipString << skip[i];
    if(i<skip.size()-1)skipString<< ",";
  }
  newTable->addKey("skip",skipString.str(),"CCDs to skip");

  dbg<<"writing arrays"<<std::endl;
  newTable->column(colName[0]).write(x,1);
  newTable->column(colName[1]).write(y,1);
  newTable->column(colName[2]).write(ccd,1);
  newTable->column(colName[3]).write(clip,1);
  newTable->column(colName[4]).writeArrays(vals,1);
  

}
  
template<class T>
bool Exposure<T>::readShapelet(std::string file,int _nvar,int _nccd)
{
  
  
  try {
    dbg << "opening file " << file<<std::endl;
    std::auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(file.c_str(),CCfits::Read));
    
    CCfits::ExtHDU& table = pInfile->extension(1);
    double x_max,y_max;
    double x_min,y_min;
    string skip_string;
    table.readKey("name",label);
    table.readKey("addSize",addSize);
    table.readKey("nvar",nvar);
    table.readKey("nccd",max_ccd);
    table.readKey("skip",skip_string);
    table.readKey("sStart",shapeStart);
    table.readKey("x_max",x_max);
    table.readKey("y_max",y_max);
    table.readKey("x_min",x_min);
    table.readKey("y_min",y_min);

    
    table.readKey("incMiss",includeMiss);

    if(_nvar>nvar) {
      std::cerr<<"Cannot read more vars"<<std::endl;
      return 0;
    }
    if(_nvar>0) nvar=_nvar;
    if(_nccd>0) max_ccd=_nccd;

    dbg<<"read keys"<<std::endl;
    skip.clear();
    std::vector<string> svals;
    Tokenize(skip_string,svals,",");
    for(int i=0;i<svals.size();++i) skip.push_back(atoi(svals[i].c_str()));
    
    long nrows=table.rows();
    dbg << "found " << nrows<<" objects"<<std::endl;
    long start=1;
    long end=nrows;
    
    std::vector<int> ccd;
    std::vector<bool> clip;
    std::vector<T> xpos;
    std::vector<T> ypos;
    std::vector<std::valarray<T> > vals;
    
    table.column("clip").read(clip, start, end);
    table.column("ccd").read(ccd, start, end);
    table.column("x").read(xpos, start, end);
    table.column("y").read(ypos, start, end);
    table.column("vals").readArrays(vals, start, end);
      
    
    for(int ichip=1;ichip<=max_ccd;++ichip) {
      
      // check if this chip should be skipped
      std::vector<int>::iterator iter=find(skip.begin(),skip.end(),ichip);
      if(iter!=skip.end()) continue;

      Chip<T> *chip=new Chip<T>(ichip,x_min,x_max,y_min,y_max);
      addChip(ichip,chip);
    }

    for(int i=0;i<nrows;++i) {
      if(ccd[i]>max_ccd) continue;
      
      Detection<T> *det=new Detection<T>(xpos[i],ypos[i],nvar);
      det->setClip(clip[i]);
      for(int j=0;j<nvar;++j) {
	(*det)(j)=vals[i][j];
      }
      chips[ccd[i]]->addDet(det);
    }
    
  }
  
  catch (CCfits::FitsException& e) {
      
    std::cout<<"Can't read file "<<file<<" "<<e.message()<<std::endl;
    return false;
  }

  ChipIter c_it=ChipBegin();
  for( ; c_it!=ChipEnd();++c_it) {
    if(c_it->second->getNDet()==0 && !includeMiss) {
      std::cout<<"Chip "<<c_it->second->getLabel()<<" is missing objects."<<std::endl;
      return false;
    }
  }
  return true;

}

template class Detection<float> ;
template class Detection<double> ;
template class Chip<float> ;
template class Chip<double> ;
template class Exposure<float> ;
template class Exposure<double> ;
