#include "PCAObjects.h"
#include <CCfits/CCfits>
#include "TMV.h"
#include <sstream>
#include "myTypeDef.h"
#include <cassert>
#include "Log.h"
#include "Image.h"
#include "NR.h"
namespace PCA {

  using std::cout;
  using std::endl;

  // reeturn legendre polynomial of specified order at single point centered on image
  static DVector definePXY(int order, float x, float xmin, float xmax)
  {
    DVector temp(order+1);
    double newx = (2.*x-xmin-xmax)/(xmax-xmin);
    temp[0] = 1.;
    if(order>0) temp[1] = newx;
    for(int i=2;i<=order;++i) {
      temp[i] = ((2.*i-1.)*newx*temp[i-1] - (i-1.)*temp[i-2])/i;
    }
    return temp;
  }
  
  // return matrix of x,y points 
  static void setPRow(int fitorder, Position<float> pos, const Bounds<float>& bounds, DVectorView prow)
  {
    //Assert(int(prow.size()) == (fitorder+1)*(fitorder+2)/2);
    DVector px =
      definePXY(fitorder,pos.x,bounds.getXMin(),bounds.getXMax());
    DVector py =
      definePXY(fitorder,pos.y,bounds.getYMin(),bounds.getYMax());
    int pq = 0;
    for(int n=0;n<=fitorder;++n) {
      for(int p=n,q=n-p;q<=n;--p,++q) {
        //Assert(pq < int(prow.size()));
        prow(pq) = px[p]*py[q];
        ++pq;
      }
    }
    //Assert(pq == int(prow.size()));
  }



  template<class T>
  int Cell<T>::getNGood() {
    int ngood=0;
    for(int i=0;i<dets.size();++i) {
      if(!dets[i]->isClipped()) ngood++;
    }
    return ngood;
  }

  template<class T>
  int Cell<T>::getNClip() {
    int nclip=0;
    for(int i=0;i<dets.size();++i) {

    }
    return nclip;
  }
  // Return the mean values of all the detections in a cell
  // if no detections found it will return a zero
  template<class T>
  std::vector<T> Cell<T>::getMeanVals()
  {

    int ngood=this->getNGood();
    std::vector<T> v(nvar,defaultVal);
    if (ngood==0) {
       FILE_LOG(logDEBUG1)<<"  this cell does not have at least one detections. "<<endl;
       this->setMissing(true);
       return v;
    }

    for(int j=0;j<nvar;++j) v[j]=0;

    for(int i=0;i<dets.size();++i) {
      if(dets[i]->isClipped()) continue;
      for(int j=0;j<nvar;++j) {
        v[j]+=dets[i]->getVal(j);
      }
    }
    
    
    for(int j=0;j<nvar;++j) {
      v[j]/=dets.size();
      
      FILE_LOG(logDEBUG1)<<"  Mean Cell var :"<<j<<" "<<v[j]<<endl;
    }

    return v;
    
  }
  
    
  
  // Return the median values of all the detections in a cell
  template<class T>
  std::vector<T> Cell<T>::getMedianVals()
  {
    std::vector<T> v(nvar,defaultVal);
    int ngood=this->getNGood();
    if(ngood<3) {
      FILE_LOG(logDEBUG1)<<"  this cell does not have at least three detections. "<<endl;
      this->setMissing(true);
      return v;
    }
    for(int j=0;j<nvar;++j) v[j]=0;

    for(int j=0;j<nvar;++j) {
      std::vector<T> tmp;      
      for(int i=0;i<dets.size();++i) {
	if(dets[i]->isClipped()) continue;
        tmp.push_back(dets[i]->getVal(j));
      }
      v[j]=median<T>(tmp);

      FILE_LOG(logDEBUG2)<<"Median Cell var "<<j<<" "<<v[j]<<endl;
    }

    return v;
    
  }

  // Return the median values of all the detections in a cell
  template<class T>
  std::vector<T> Cell<T>::getMeanClipVals(float clip)
  {

    std::vector<T> v(nvar,defaultVal);
    FILE_LOG(logDEBUG1)<<"   Meanclips with :"<<dets.size()<<endl;
    int ngood=this->getNGood();
    if(ngood<1) {
      FILE_LOG(logDEBUG1)<<"  this cell does not have at least one detections. "<<endl;
      this->setMissing(1);
      return v;
    }

    for(int j=0;j<nvar;++j) v[j]=0;

    //iterate once through the variables to label outliers in any 
    // of the variables using the the median absolute deviation
    int nclip=0;
    for(int j=0;j<nvar;++j) {

      FILE_LOG(logDEBUG1)<<"   Getting variable: "<<j<<endl;
      std::vector<T> tmp(ngood);     
      
      int cur=0;
      for(int i=0;i<dets.size();++i) {
	if(dets[i]->isClipped()) continue;
	tmp[cur]=(dets[i]->getVal(j));
	cur++;
      }
      double mad;
      double median=median_mad(tmp,mad);
      FILE_LOG(logDEBUG1)<<"   Found "<<tmp.size()
			 <<" objects with median: "<<median
			 <<" mad: "<<mad<<endl;

      cur=0;
      int ngood=0;
      for(int i=0;i<dets.size();++i) {
	if(dets[i]->isClipped()) continue;

	if( std::abs(dets[i]->getVal(j)-median) > clip*mad) {
          FILE_LOG(logDEBUG2)<<"   Clipping object "<<i
			     <<" diff from median: "
			     << std::abs(dets[i]->getVal(j)-median)
			     <<" max allowed:"<<clip*mad<<endl;
          nclip++;
        }
	else {
	  v[j]+=dets[i]->getVal(j);
	  FILE_LOG(logDEBUG2)<<"     adding "<<dets[i]->getVal(j)<<endl;
	  ngood++;
	}
	
      }
      FILE_LOG(logDEBUG1)<<"   Clipped "<<nclip<<" detections"<<endl;
      v[j]/=ngood;
    }

    
    std::stringstream str;
    str<<"   Final means :";

    for(int j=0;j<nvar;++j) {
      str<<v[j]<<" ";
    }
    FILE_LOG(logDEBUG1)<<str.str()<<endl;
    return v;

  }



  // Return a polynomial fit of all the detections in a cell
  template<class T>
  std::vector<T> Cell<T>::getFitVals(int order,float clip)
  {
    fitorder=order;
    int nfit=(fitorder+1)*(fitorder+2)/2;
    int ndet=this->getNGood();
    std::vector<T> v(nvar*nfit,defaultVal);    



    if(ndet<nfit) {
      FILE_LOG(logDEBUG1)<<"  this cell has more parameters than detections "
			 <<ndet<<" "<<nfit<<".  "<<endl;
      this->setMissing(true);
      return v;
    }
    

    FILE_LOG(logDEBUG)<<"getting fit vals for "<<ndet<<" detections "<<nfit<<endl;



    // we need to calculate first which objects we may need
    // to reject.  We don't necesarily want to clip them because
    // they may be useful in better iterations
    std::vector<bool> use_det(ndet,true);
    for(int j=0;j<nvar;++j) {

      std::vector<T> tmp;     
      
      for(int i=0;i<dets.size();++i) {
	if(dets[i]->isClipped()) continue;
	tmp.push_back(dets[i]->getVal(j));
      }
      double mad;
      double median=median_mad(tmp,mad);

      
      for(int i=0;i<dets.size();++i) {
	if(dets[i]->isClipped()) continue;

	// check that is within 5 sigma
	if(std::fabs(dets[i]->getVal(j)-median)>clip*mad) {
	   FILE_LOG(logDEBUG)<<"skipping "<<i<<" value to large "
			     <<std::fabs(dets[i]->getVal(j)-median)/mad<<" sigma"<<endl;
 	  use_det[i]=false;
	}
	
      }
    }
    
    int ndet_cur=0;
    for(int i=0;i<dets.size();++i) if(use_det[i]) ndet_cur++;

    if(ndet_cur<nfit) {
      FILE_LOG(logDEBUG1)<<"  this cell has more parameters than detections "
			 <<ndet<<" "<<nfit<<".  "<<endl;
      this->setMissing(true);
      return v;
    }

    FILE_LOG(logDEBUG)<<"using "<<ndet_cur<<" objects "<<endl;
    DMatrix b(ndet_cur,nvar);
    DMatrix A(ndet_cur,nfit);
    


    for(int j=0;j<nvar;++j) {
      int cur=0;    
      for(int i=0;i<dets.size();++i) {
	if(dets[i]->isClipped() || !use_det[i]) continue;
	b(cur,j)=dets[i]->getVal(j);
	setPRow(fitorder,dets[i]->getPos(),bounds,A.row(cur));
	cur++;
      }
    }
    FILE_LOG(logDEBUG)<<"data "<<b<<endl;

    FILE_LOG(logDEBUG)<<"poses "<<A<<endl;
    DMatrix x=b/A;//b.subMatrix(0,cur,0,nvar)/A.subMatrix(0,cur,0,ndet);
    FILE_LOG(logDEBUG)<<"x "<<x<<endl;


    int cur=0;
    for(int i=0;i<x.nrows();++i) {
      for(int j=0;j<x.ncols();++j) {
	
        v[cur]=x(i,j);
	cur++;
      }
    }


    return v;

  }
      
  template<class T>
  std::vector<T> Cell<T>::getVals(std::string type,std::vector<float> &params)
  {
    if(getTypeFromString(type)==Mean) {
      return getMeanVals();
    }
    else if(getTypeFromString(type)==MeanClip) {
      assert(params.size()>=1);
      float sigma_clip=params[0];
      return getMeanClipVals(sigma_clip);
    }
    else if(getTypeFromString(type)==Median) {
      return getMedianVals();
    }
    else if(getTypeFromString(type)==Fit) {
      assert(params.size()>=1);
      int order=static_cast<int>(params[0]);
      float clip=params[1];
      return getFitVals(order,clip);
    }
    
  }


  
  template<class T>
  int Cell<T>::getNVal(std::string type,std::vector<float> &params)
  {
  
    if(getTypeFromString(type)==Mean) {
      return nvar;
    }

    else if(getTypeFromString(type)==MeanClip) {
      return nvar;
    }
    else if(getTypeFromString(type)==Median) {
      return nvar;
    }
    else if(getTypeFromString(type)==Fit) {
      assert(params.size()>=1);
      int fitorder=static_cast<int>(params[0]);
      return nvar*(fitorder+1)*(fitorder+2)/2;
    }
    return 0;
  }


  template<class T>
  std::vector<std::vector<T> > Cell<T>::getDiff(tmv::ConstVectorView<T> &vals,std::string type,
						std::vector<float> params,
						const std::vector<double> &mean,
						const std::vector<double> &sigma,
						bool clip,
						double nclip)
  {
    std::vector<std::vector<T> > vdiff(nvar);
    if(this->isMissing()) return  vdiff;

    if(getTypeFromString(type)!=Fit) {
      for(int j=0;j<nvar;++j) {
	for(int i=0;i<dets.size();++i) {
	  if(dets[i]->isClipped()) continue;

	  if(!clip) {
	    FILE_LOG(logDEBUG1)<<"  nvar "<<j<<" det "<<i<<" "
			       <<vals[j]<<" "<<dets[i]->getVal(j)<<" "
			       <<vals[j]-dets[i]->getVal(j)<<endl;
	    
	    vdiff[j].push_back(vals[j]-dets[i]->getVal(j));
	  }
	  else {
	    double diff=std::fabs((vals[j]-dets[i]->getVal(j))-mean[j])/sigma[j];
	    if(diff>nclip) {
	     
	      FILE_LOG(logDEBUG1)<<"  setting var "<<j<<" det "<<i<<" as clipped "
				 <<diff<<" compared to "<<nclip<<endl;
	      
	      
	      dets[i]->setClip(true);
	    }
	    else  vdiff[j].push_back(vals[j]-dets[i]->getVal(j));

	  }
	}
      }
    }
    else {
      int fitorder=params[0];
      int nfit=(fitorder+1)*(fitorder+2)/2;
      int ndet=this->getNGood();
      
      DMatrix br(ndet,nvar);
      DMatrix A(ndet,nfit);
      DMatrix x(nfit,nvar);
      
      int cur=0;
      for(int n=0;n<dets.size();++n) {
	if(dets[n]->isClipped()) continue;
	setPRow(fitorder,dets[n]->getPos(),bounds,A.row(cur));
	cur++;
      }
      cur=0;
      for(int i=0;i<nfit;++i) {
	for(int j=0;j<nvar;++j) {
	  x(i,j)=vals[cur];
	  cur++;
	}
      }
      FILE_LOG(logDEBUG1)<<"x "<<x<<endl;
      FILE_LOG(logDEBUG1)<<"A "<<A<<endl;
      br=A*x;

      for(int j=0;j<nvar;++j) {
	cur=0;
	for(int i=0;i<dets.size();++i) {
	  if(dets[i]->isClipped()) continue;
	  
	  if(!clip) {
	    FILE_LOG(logDEBUG1)<<"  fit nvar "<<j<<" det "<<i<<" "
			       <<br(cur,j)<<" "<<dets[i]->getVal(j)
			       <<" "<<br(cur,j)-dets[i]->getVal(j)<<endl;
			
	    vdiff[j].push_back(br(cur,j)-dets[i]->getVal(j));
	  }
	  else {
	    double diff=std::fabs( (br(cur,j)-dets[i]->getVal(j))-mean[j])/sigma[j];
	    if(diff>nclip) {
	      dets[i]->setClip(true);
	      FILE_LOG(logDEBUG1)<<"  setting var "<<j<<" det "<<i<<" as clipped "
				 <<diff<<" compared to "<<nclip<<endl;
	    }
	    else vdiff[j].push_back(br(cur,j)-dets[i]->getVal(j));
	  }
	  cur++;
	}
              
      }
      
    }
    return vdiff;
    
  }

  template<class T>
  std::vector<std::valarray<T> > Cell<T>::getDetVals(tmv::ConstVectorView<T> &vals,
						     std::string type,std::vector<float> params)
  {
    std::vector<std::valarray<T> > vdet;

    if(getTypeFromString(type)!=Fit) {
      
	for(int i=0;i<dets.size();++i) {
	  std::valarray<T> sdet(nvar);
	  for(int j=0;j<nvar;++j) {
	    sdet[j]=vals[j];
	  }
	  vdet.push_back(sdet);
	}
    }
    else {
      int fitorder=params[0];
      int nfit=(fitorder+1)*(fitorder+2)/2;
      int ndet=this->getNDet();
      DMatrix br(ndet,nvar);
      DMatrix A(ndet,nfit);
      DMatrix x(nfit,nvar);
      
      for(int n=0;n<dets.size();++n) {
	setPRow(fitorder,dets[n]->getPos(),bounds,A.row(n));
      }
      int cur=0;
      for(int i=0;i<nfit;++i) {
	for(int j=0;j<nvar;++j) {
	  x(i,j)=vals[cur];
	  cur++;
	}
      }

      br=A*x;
      
      for(int i=0;i<dets.size();++i) {
	std::valarray<T> sdet(nvar);
	for(int j=0;j<nvar;++j) {
	  sdet[j]=br(i,j);
	}
	vdet.push_back(sdet);
      }
    }

    return vdet;
    
  }
  




  template<class T>
  int Chip<T>::getNClip()  
  {
    int nclip=0;
    for(int i=0;i<cells.size();++i) nclip+=cells[i]->getNClip();
    return nclip;
  }

  template<class T>
  int Chip<T>::getNGood()
  {
    int ngood=0;
    for(int i=0;i<cells.size();++i) ngood+=cells[i]->getNGood();
    return ngood;
  }

  template<class T>
  int Chip<T>::getNDet()
  {
    int ndet=0;
    for(int i=0;i<cells.size();++i) ndet+=cells[i]->getNDet();
    return ndet;
  }
       

  // Get the mean values of all the cells in a chip
  // The ordering of the variables are 
  // Ce1l 1 var1..varN, Cell2 var1..varN, Cell3...
  template<class T>
  std::vector<T> Chip<T>::getVals(std::string type,std::vector<float> &params)
  { 

    // assume all cells have the same number
    int ntotvar=cells[0]->getNVal(type,params);
    std::vector<T> v(ntotvar*cells.size());
    int cur_index=0;
    for(int i=0;i<cells.size();++i) {
      FILE_LOG(logDEBUG)<<"  Get vals from cell "<<i<<endl;
      std::vector<T> cv=cells[i]->getVals(type,params);
      for(int j=0;j<ntotvar;++j) {
        v[cur_index]=cv[j];
        cur_index++;
      }
    }
    return v;
  }

  template<class T>
  std::vector<bool> Chip<T>::getMissing()
  {
    std::vector<bool> v(cells.size(),false);

    for(int i=0;i<cells.size();++i) {
      if(cells[i]->isMissing()) v[i]=true;
    }
    return v;
  }

  template<class T>
  void Chip<T>::setMissing(float prob)
  {
    for(int i=0;i<cells.size();++i) {
      if(ran01()<prob) cells[i]->setMissing(true);
    }
  }

  
  template<class T>
  void Chip<T>::divide(int nvar,int _nx,int _ny) {
    nx=_nx;
    ny=_ny;
    std::vector<Bounds<float> > vb=bounds.divide(nx,ny);
    for(int i=0;i<vb.size();++i) {
      Cell<T> *cell=new Cell<T>(nvar,vb[i]);
      cells.push_back(cell);
      cbounds.push_back(vb[i]);
    }
  }
  
  template<class T>
  void Chip<T>::addDet(Detection<T> *det) {

    // The bounds class has the y as the fast moving coordinate
    int bin_x=static_cast<int>(det->getPos().x/(bounds.getXMax()/nx));
    int bin_y=static_cast<int>(det->getPos().y/(bounds.getYMax()/ny));
    int bin=bin_x*ny+bin_y;

    cells[bin]->addDet(det);
  }


  template<class T>
  Exposure<T>::Exposure (string _label,int _nchip, int _shapestart,double _ra,double _dec,float _airmass):
    label(_label),nchip(_nchip),ra(_ra),dec(_dec),airmass(_airmass),
    nx_chip(-1.),ny_chip(-1.),xmax_chip(-1.),ymax_chip(-1.),
    shapeStart(_shapestart),outlier(0) {}

  template<class T>
  int Exposure<T>::getNClip()
  {
    int nclip=0;
    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();

    for(; iter!=chips.end();++iter) nclip+=iter->second->getNClip();
    return nclip;
  }

  template<class T>
  int Exposure<T>::getNGood()
  {
    int ngood=0;
    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();

    for(; iter!=chips.end();++iter) ngood+=iter->second->getNGood();
    return ngood;
  }



  template<class T>
  int Exposure<T>::getNDet()
  {
    int ndet=0;
    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();

    for(; iter!=chips.end();++iter) ndet+=iter->second->getNDet();
    return ndet;
  }


  template<class T>
  bool Exposure<T>::readShapelet(std::string dir,int nvar,bool add_size,
				 bool include_miss,bool use_dash,string suffix,
				 std::string exp) {
    if (exp.empty()) exp=label;
    FILE_LOG(logINFO) << "Reading exposure " << exp<<endl;
    for(int ichip=1;ichip<=nchip;++ichip) {
      
      Chip<T> *chip=new Chip<T>(ichip,xmax_chip,ymax_chip);
      chip->divide(nvar,nx_chip,ny_chip);

      // check if this chip should be skipped
      std::vector<int>::iterator iter=find(skip.begin(),skip.end(),ichip);
      if(iter!=skip.end()) continue;
      
      std::stringstream inputFile;
      if(!use_dash) {
        inputFile << dir << "/" << exp << "_";
      }
      else {
        inputFile << dir << "/" << exp << "-";
      }
      if(ichip<10) inputFile <<0;
      
      if(!use_dash) {
        inputFile << ichip << "_"+suffix;
      }
      else {
        inputFile << ichip << "-"+suffix;
      }


      try {
	FILE_LOG(logDEBUG) << "opening file " << inputFile.str()<<endl;
        std::auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(inputFile.str(),CCfits::Read));
        
        CCfits::ExtHDU& table = pInfile->extension(1);
        
        long nTabRows=table.rows();
	FILE_LOG(logDEBUG) << "found " << nTabRows<<" objects"<<endl;
        long start=1;
        long end=nTabRows;
        
        std::vector<int> psf_flags;
	std::vector<double> psf_size;
        std::vector<double> xpos;
        std::vector<double> ypos;
        
        table.column("psf_flags").read(psf_flags, start, end);
        table.column("sigma_p").read(psf_size, start, start+1);
        table.column("x").read(xpos, start, end);
        table.column("y").read(ypos, start, end);
        
        std::vector<long> order;       // shapelet order
        table.column("psf_order").read(order, start, end);
        
        int icount=0;
        for (int i=0; i<nTabRows; i++) {
          if (!psf_flags[i]) {          // pass psf flags
            
            int row=i+1;
            Detection<T> *det=new Detection<T>(xpos[i],ypos[i],nvar);
            int ncoeff=(order[i]+1)*(order[i]+2)/2;      // psf values
            std::valarray<double> coeffs;
            table.column("shapelets").read(coeffs, row); 
            FILE_LOG(logDEBUG1)<<"adding object "<<ichip<<" "<<xpos[i]<<" "
                               <<ypos[i]<<" "<<coeffs[shapeStart]<<" "
                               <<coeffs[shapeStart+1]<<" "<<coeffs[shapeStart+2]<<" "<<endl;
            for(int j=0;j<nvar;++j) {
              
              det->setVal(j,coeffs[shapeStart+j]);
            }
            if(add_size) {
	      FILE_LOG(logDEBUG1)<<"adding size "<<psf_size[0]<<endl;
	      det->setVal(0,psf_size[0]);
	    }
            chip->addDet(det);
          }
        }
        
      }
      catch (CCfits::FitsException& ) {

	
        if(!include_miss) {
	  FILE_LOG(logERROR)<<"Can't open chip: "<<inputFile.str()
			    <<" from exposure "<<exp<<" skipping"<<endl;
	  return false;
	}
	else {
	  FILE_LOG(logERROR)<<"Can't open chip: "<<inputFile.str()
			    <<" from exposure "<<exp<<" will be set to missing"<<endl;
	}
      }
      addChip(ichip,chip);
      
    }

    return true;
  }
  

  template<class T>
  bool Exposure<T>::readPixels(std::string dir,int npix, int nvar,std::string sdir,
			       bool use_dash,std::string exp) {
    if (exp.empty()) exp=label;
    FILE_LOG(logINFO) << "Reading exposure " << exp<<endl;
    for(int ichip=1;ichip<=nchip;++ichip) {
      
      Chip<T> *chip=new Chip<T>(ichip,xmax_chip,ymax_chip);
      chip->divide(nvar,nx_chip,ny_chip);

      // check if this chip should be skipped
      std::vector<int>::iterator iter=find(skip.begin(),skip.end(),ichip);
      if(iter!=skip.end()) continue;
      
      std::stringstream inputFile;
      inputFile << dir << "/" << exp << "_";
      if(ichip<10) inputFile <<0;

      inputFile << ichip << ".fits.fz";      

      string image_file=inputFile.str();

      std::stringstream inputFile2;
      if(!use_dash) {
        inputFile2 << sdir << "/" << exp << "_";
      }
      else {
        inputFile2 << sdir << "/" << exp << "-";
      }
      if(ichip<10) inputFile2 <<0;
      
      if(!use_dash) {
        inputFile2 << ichip << "_psf.fits";
      }
      else {
        inputFile2 << ichip << "-psf.fits";
      }
      
      string psf_file=inputFile2.str();

      FILE_LOG(logDEBUG) << "Reading image file " << image_file<<endl;

      Image<T> im (image_file,2); // main image
      Image<T> wim(image_file,4); // weight image
      Image<T> bpm(image_file,3); // bad pixel mask

      FILE_LOG(logDEBUG) << "Reading psf file " << psf_file<<endl;
      std::vector<Position<T> > pos;
      std::vector<double> use_sky;
      try {
	std::auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(psf_file,CCfits::Read));
	
	CCfits::ExtHDU& table = pInfile->extension(1);
	
	
	long nTabRows=table.rows();
	FILE_LOG(logDEBUG) << "found " << nTabRows<<" objects"<<endl;
	long start=1;
	long end=nTabRows;
	
	std::vector<int> psf_flags;
	std::vector<int> sky;
	std::vector<double> xpos;
	std::vector<double> ypos;
	
	table.column("psf_flags").read(psf_flags, start, end);
	table.column("x").read(xpos, start, end);
	table.column("y").read(ypos, start, end);
	table.column("sky").read(sky, start, end);
	
	std::vector<long> order;       // shapelet order
	
	table.column("psf_order").read(order, start, end);
	for (int i=0; i<nTabRows; i++) {
	  if (!psf_flags[i]) {          // pass psf flags
            
	    pos.push_back(Position<T>(xpos[i],ypos[i]));
	    use_sky.push_back(sky[i]);
	  }
	}
      }
      catch (CCfits::FitsException& ) {
	FILE_LOG(logERROR)<<"Can't open chip: "<<inputFile.str()<<" from exposure "<<exp<<" skipping"<<endl;
	return false;
      }
     
      int app=npix;
      int npixu=4*(app+1)*(app+1);
      assert(npixu==nvar);
      std::vector<double> pixels(nvar),weight(nvar);
      // Loop over centers of objects and get sky-subtracted pixel lists
      for(int istar=0;istar<pos.size();++istar) {

	double xcen=pos[istar].x;
	double ycen=pos[istar].y;
	FILE_LOG(logDEBUG1)<<"Getting pixels around star : "<<istar<<" "
			   <<pos[istar]<<" "<<endl;

	int i1 = int(floor(xcen-app));
	int i2 = int(ceil(xcen+app));
	int j1 = int(floor(ycen-app));
	int j2 = int(ceil(ycen+app));

	if (i1 < 0) { i1 = 0; }
	if (i1 > xmax_chip) { i1 = xmax_chip; }
	if (j1 < 0) { j1 = 0; }
	if (j1 > ymax_chip) { i1 = ymax_chip; }

	double chipx = i1-xcen;
	double peak=0;
	int cpix=0;;
	FILE_LOG(logDEBUG1)<<"Boundary x : "<<i1<<","<<i2<<"  y: "<<j1<<","<<j2<<endl;

	for(int i=i1;i<=i2;++i) {
	  double chipy = j1-ycen;
	  for(int j=j1;j<=j2;++j) {
	    double flux = im(i,j)-use_sky[istar];
	    double inverseVariance=wim(i,j);
	    double bp=bpm(i,j);

 	    if(bp>0.0) {
	      inverseVariance=0;
	      flux=-999.0;
	    }
	    
	   //  FILE_LOG(logDEBUG2)<<"Values at pixel "<<i<<","<<j<<","<<cpix<<" "
// 			       <<flux<<" "<<inverseVariance<<" "<<bp<<endl;
	    pixels[cpix]=flux;
	    weight[cpix]=inverseVariance;
 	    cpix++;
	  }
	}

	Detection<T> *det=new Detection<T>(pos[istar].x,pos[istar].y,nvar);
	for(int j=0;j<nvar;++j) {
	  det->setVal(j,pixels[j]);
	}
	chip->addDet(det);
      }
      addChip(ichip,chip);
    }
    return true;
  }
  



  template<class T>
  std::vector<bool> Exposure<T>::getMissing()
  {
    int nchip_var=ny_chip*nx_chip;
    int ncell=chips.size()*ny_chip*nx_chip;
    std::vector<bool> v(ncell,false);
    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();
    int cur_index=0;    
    for(; iter!=chips.end();++iter) {
      std::vector<bool> cv=iter->second->getMissing();
      for(int j=0;j<cv.size();++j) {
	v[cur_index]=cv[j];
	cur_index++;
      }
    }
    
    return v;
  }




  template<class T>
  void Exposure<T>::setMissing(float prob)
  {
    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();
    
    for(; iter!=chips.end();++iter) iter->second->setMissing(prob);
  }

  template<class T>
  tmv::Vector<T> Exposure<T>::getVals(std::string type,std::vector<float> &params)
  { 
    int nchip_var=ny_chip*nx_chip;
    int nfocal=chips.size()*nchip_var;
    

    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();

    int ntotvar=iter->second->getCell(0)->getNVal(type,params);
    FILE_LOG(logDEBUG)<<"Getting data from "<<label
		      <<" total vars: "<<ntotvar*nfocal<<" chips: "
		      <<chips.size()<<endl;
    tmv::Vector<T> v(ntotvar*nfocal,0.0);
    int cur_index=0;
    int cur_chip=0;
    
    for(; iter!=chips.end();++iter,cur_chip++) {
      FILE_LOG(logDEBUG)<<" Get vals from chip: "<<iter->first<<endl;
      std::vector<T> cv=iter->second->getVals(type,params);
      for(int j=0;j<cv.size();++j) {

        // complicated to match current structure
        // this is probably wrong
        //FILE_LOG(logINFO)<<"  "<<j<<" "<<nfocal*(j%ntotvar)+cur_chip*nchip_var+j/(ntotvar)<<endl;
        //v[nfocal*(j%ntotvar)+cur_chip*nchip_var+j/(ntotvar)]=cv[j];
	
        // I want to switch to this once I am willing to change all my 
        // plotting routines
        v[cur_index]=cv[j];
        cur_index++;
      }

    }
    return v;
  }

  template<class T>
  std::vector<double> Exposure<T>::outlierReject(const tmv::Vector<T> &data_r,
						 double sigma,string type,
						 std::vector<float> params)
  {
    typename std::map<int,Chip<T>*>::const_iterator iter=chips.begin();
    int cur=0;
    std::vector<std::vector<double> > diff_all;
    FILE_LOG(logDEBUG)<<"Exposure "<<label<<" outlier "<<endl;

    int ichip=0;
    int nperchip;
    int nvartot;
    // compute median and deviation for the exposure
    for(; iter!=chips.end();++iter,++ichip) {

      if(ichip==0) {
	nvar=iter->second->getCell(0)->getNVar();
	nperchip=nx_chip*ny_chip*nvar;
	nvartot=nvar;
	if(type=="fit") {
	  int order=(params[0]+1)*(params[0]+2)/2;
	  nvartot=nvar*order;
	  nperchip*=order;
	}
      }

      
      tmv::ConstVectorView<T> data_chip=data_r.subVector(ichip*nperchip,(ichip+1)*nperchip);

      FILE_LOG(logDEBUG)<<" Finding Outliers Chip "<<iter->first<<endl
	;
      for(int icell=0;icell<iter->second->getNCell();++icell) {

	
	tmv::ConstVectorView<T> data_cell=data_chip.subVector(icell*nvartot,(icell+1)*nvartot);
		
	// do not test missing data that was added later or that may have -999
	if(iter->second->getCell(icell)->isMissing()) continue;
	
	std::vector<double> d1,d2;
	std::vector<std::vector<double> > diff=
	  iter->second->getCell(icell)->getDiff(data_cell,type,params,d1,d2,false);

	for(int ivar=0;ivar<nvar;++ivar) {
	  std::vector<double> tmp;

	  copy(diff[ivar].begin(),diff[ivar].end(),std::back_inserter(tmp));
	  double dmad;
	  double dmedian=median_mad(tmp,dmad);
	  FILE_LOG(logDEBUG)<<"  ivar "<<ivar<<" vcell1 "<<icell<<" :"
			    <<dmedian<<" "<<dmad<<endl;
	  diff_all.push_back(tmp);
	}
	
      }
    }
    
    // now remove outliers
    std::vector<double> mad(nvar);
    std::vector<double> median(nvar);
    for(int ivar=0;ivar<nvar;++ivar) {
      median[ivar]=median_mad(diff_all[ivar],mad[ivar]);
      FILE_LOG(logDEBUG)<<"   ivar "<<ivar<<" totalres1: "<<median[ivar]<<" :"<<mad[ivar]<<endl;
    }

    iter=chips.begin();
    ichip=0;
    diff_all.clear();
    for(; iter!=chips.end();++iter,ichip++) {
      tmv::ConstVectorView<T> data_chip=data_r.subVector(ichip*nperchip,(ichip+1)*nperchip);
      FILE_LOG(logDEBUG)<<" removing from Chip "<<iter->first<<endl;
      for(int icell=0;icell<iter->second->getNCell();++icell) {
	tmv::ConstVectorView<T> data_cell=data_chip.subVector(icell*nvartot,(icell+1)*nvartot);
	if(iter->second->getCell(icell)->isMissing()) continue;

	
	std::vector<std::vector<double> > diff=
	  iter->second->getCell(icell)->getDiff(data_cell,type,params,
						median,mad, true,sigma);
	
	for(int ivar=0;ivar<nvar;++ivar) {
	  std::vector<double> tmp;
	  copy(diff[ivar].begin(),diff[ivar].end(),std::back_inserter(tmp));
	  double dmad;
	  double dmedian=median_mad(diff[ivar],dmad);
	  diff_all.push_back(tmp);
	  FILE_LOG(logDEBUG)<<"  ivar "<<ivar<<"vcell2 "<<icell<<" :"<<dmedian<<" "<<dmad<<endl;
	}
      }
    }
    
    for(int ivar=0;ivar<nvar;++ivar) {
      median[ivar]=median_mad(diff_all[ivar],mad[ivar]);
      FILE_LOG(logDEBUG)<<"   ivar "<<ivar<<" totalres2: "<<median[ivar]<<" :"<<mad[ivar]<<endl;
    }
    return median;
    
  }
  

  //template class Detection<float> ;
  template class Detection<double> ;
  // template class Cell<float> ;
  template class Cell<double> ;
  //template class Chip<float> ;
  template class Chip<double> ;
  //template class Exposure<float> ;
  template class Exposure<double> ;

  
};
    
