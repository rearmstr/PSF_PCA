#include "PCAObjects.h"
#include <CCfits/CCfits>
#include "TMV.h"
#include <sstream>
#include "myTypeDef.h"
namespace PCA {

  using std::cout;
  using std::endl;

  static FVector definePXY(int order, float x, float xmin, float xmax)
  {
    FVector temp(order+1);
    double newx = (2.*x-xmin-xmax)/(xmax-xmin);
    temp[0] = 1.;
    if(order>0) temp[1] = newx;
    for(int i=2;i<=order;++i) {
      temp[i] = ((2.*i-1.)*newx*temp[i-1] - (i-1.)*temp[i-2])/i;
    }
    return temp;
  }

  static void setPRow(int fitorder, Position<float> pos, const Bounds<float>& bounds, FVectorView prow)
  {
    //Assert(int(prow.size()) == (fitorder+1)*(fitorder+2)/2);
    FVector px =
      definePXY(fitorder,pos.x,bounds.getXMin(),bounds.getXMax());
    FVector py =
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




  int Cell::getNDet() {
    int ngood=0;
    for(int i=0;i<dets.size();++i) {
      if(!dets[i]->isClipped()) ngood++;
    }
  }
  // Return the mean values of all the detections in a cell
  // if no detections found it will return a zero
  std::vector<float> Cell::getMeanVals()
  {
    std::vector<float> v(nvar,defaultVal);
    if(dets.size()==0) return v;
        
    for(int j=0;j<nvar;++j) v[j]=0;

    for(int i=0;i<dets.size();++i) {
      for(int j=0;j<nvar;++j) {
        v[j]+=dets[i]->getVal(j);
      }
    }


    for(int j=0;j<nvar;++j) {
      v[j]/=dets.size();
    }

    return v;

  }

  // Return the median values of all the detections in a cell
  std::vector<float> Cell::getMedianVals()
  {
    std::vector<float> v(nvar,defaultVal);
    if(dets.size()==0) return v;
    for(int j=0;j<nvar;++j) v[j]=0;

    for(int j=0;j<nvar;++j) {
      std::vector<float> tmp(dets.size());      
      for(int i=0;i<dets.size();++i) {
        tmp[i]=dets[i]->getVal(j);
      }
      v[j]=median<float>(tmp);
      //cout<<"    Cell var "<<j<<" "<<v[j]<<endl;
    }

    return v;

  }

  // Return the median values of all the detections in a cell
  std::vector<float> Cell::getMeanClipVals(float clip)
  {
    std::vector<float> v(nvar,defaultVal);
    if(dets.size()==0) return v;
    for(int j=0;j<nvar;++j) v[j]=0;
    //cout<<"Cell "<<endl;
    //iterate once through the variables to label outliers using the
    //the median absolute deviation
    for(int j=0;j<nvar;++j) {
      std::vector<float> tmp(dets.size());      
      for(int i=0;i<dets.size();++i) {
        tmp[i]=dets[i]->getVal(j);
      }
      double mad;
      double median=median_mad(tmp,mad);
      //cout<<nvar<<" "<<median<<" "<<mad<<endl;
      for(int i=0;i<dets.size();++i) {
        if( std::abs(tmp[i]-median) > clip*mad) {
          dets[i]->setClip(true);
          //cout<<"Clip "<<i<<" "<< std::abs(tmp[i]-median)<<" "<<clip*mad<<endl;
        }
      }
      
    }

    int ngood=0;
    // Loop through ignoring clipped guys and calculate mean
    for(int i=0;i<dets.size();++i) {
      if(dets[i]->isClipped()) continue;
      ngood++;
      for(int j=0;j<nvar;++j) {
        v[j]+=dets[i]->getVal(j);
      }
    }
    
    for(int j=0;j<nvar;++j) {
      v[j]/=ngood;
    }

    return v;

  }



  // Return the median values of all the detections in a cell
  std::vector<float> Cell::getFitVals(int order)
  {
    fitorder=order;
    int nfit=(fitorder+1)*(fitorder+2)/2;
    std::vector<float> v(nvar*nfit,defaultVal);
    if(dets.size()==0) return v;
    FMatrix b(dets.size(),nvar);
    FMatrix A(dets.size(),nfit);

    for(int j=0;j<nvar;++j) {
      std::vector<float> tmp(dets.size());      
      for(int i=0;i<dets.size();++i) {
        b(i,j)=dets[i]->getVal(j);
      }
    }
    
    for(int n=0;n<dets.size();++n) {
      setPRow(fitorder,dets[n]->getPos(),bounds,A.row(n));
    }

    FMatrix x=b/A;
    int cur=0;
    for(int i=0;i<x.nrows();++i) {
      for(int j=0;j<x.ncols();++j) {
        v[cur]=x(i,j);
        cur++;
      }
    }


    return v;

  }
      
  std::vector<float> Cell::getVals(std::string type)
  {
    if(type=="mean") {
      return getMeanVals();
    }
    else if(type=="mean_clip3") {
      return getMeanClipVals(3);
    }
    else if(type=="mean_clip4") {
      return getMeanClipVals(4);
    }
    else if(type=="median") {
      return getMedianVals();
    }
    else if(type=="plin") {
      return getFitVals(1);
    }
    else if(type=="pquad") {
      return getFitVals(2);
    }
    
  }

  int Cell::getNVal(std::string type)
  {
  
    if(type=="mean") {
      return nvar;
    }
    if(type=="mean_clip3") {
      return nvar;
    }
    if(type=="mean_clip4") {
      return nvar;
    }
    else if(type=="median") {
      return nvar;
    }
    else if(type=="plin") {
      // need to check that this is consistent
      return nvar*3;//(fitorder+1)*(fitorder+2)/2;
    }
    else if(type=="pquad") {
      // need to check that this is consistent
      return nvar*6;//(fitorder+1)*(fitorder+2)/2;
    }
  }

  // Get the mean values of all the cells in a chip
  // The ordering of the variables are 
  // Ce1l 1 var1..varN, Cell2 var1..varN, Cell3...
  std::vector<float> Chip::getVals(std::string type)
  { 

    // assume all cells have the same number
    int ntotvar=cells[0]->getNVal(type);
    std::vector<float> v(ntotvar*cells.size());
    int cur_index=0;
    for(int i=0;i<cells.size();++i) {
      //cout<<"  Get vals from cell "<<i<<endl;
      std::vector<float> cv=cells[i]->getVals(type);
      for(int j=0;j<ntotvar;++j) {
        v[cur_index]=cv[j];
        cur_index++;
      }
    }
    return v;
  }

  std::vector<bool> Chip::getMissing()
  {
    std::vector<bool> v(cells.size(),false);
    int cur_index=0;
    for(int i=0;i<cells.size();++i) {
      if(cells[i]->getNDet()==0) v[i]=true;
    }
    return v;
  }
  
  void Chip::divide(int nvar,int _nx,int _ny) {
    nx=_nx;
    ny=_ny;
    std::vector<Bounds<float> > vb=bounds.divide(nx,ny);
    for(int i=0;i<vb.size();++i) {
      Cell *cell=new Cell(nvar,vb[i]);
      cells.push_back(cell);
      cbounds.push_back(vb[i]);
    }
  }
  

  void Chip::addDet(Detection *det) {

    // The bounds class has the y as the fast moving coordinate
    int bin_x=static_cast<int>(det->getPos().x/(bounds.getXMax()/nx));
    int bin_y=static_cast<int>(det->getPos().y/(bounds.getYMax()/ny));
    int bin=bin_x*ny+bin_y;

    cells[bin]->addDet(det);
  }



  Exposure::Exposure (string _label,int _nchip, double _ra,double _dec,float _airmass):
    label(_label),nchip(_nchip),ra(_ra),dec(_dec),airmass(_airmass),
    nx_chip(-1.),ny_chip(-1.),xmax_chip(-1.),ymax_chip(-1.),shapeStart(3),outlier(0) {}

  bool Exposure::readShapelet(std::string dir,int nvar,bool use_dash,std::string exp) {
    if (exp.empty()) exp=label;
    cout << "Reading exposure " << exp<< endl;
    for(int ichip=1;ichip<=nchip;++ichip) {
      
      Chip *chip=new Chip(ichip,xmax_chip,ymax_chip);
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
        inputFile << ichip << "_psf.fits";
      }
      else {
        inputFile << ichip << "-psf.fits";
      }


      try {
        std::auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(inputFile.str(),CCfits::Read));
        
        CCfits::ExtHDU& table = pInfile->extension(1);
        
        
        long nTabRows=table.rows();
        
        long start=1;
        long end=nTabRows;
        
        std::vector<int> psf_flags;
        std::vector<double> xpos;
        std::vector<double> ypos;
        
        table.column("psf_flags").read(psf_flags, start, end);
        table.column("x").read(xpos, start, end);
        table.column("y").read(ypos, start, end);
        
        std::vector<long> order;       // shapelet order
        table.column("psf_order").read(order, start, end);
        
        int icount=0;
        for (int i=0; i<nTabRows; i++) {
          if (!psf_flags[i]) {          // pass psf flags
            
            int row=i+1;
            Detection *det=new Detection(xpos[i],ypos[i],nvar);
            int ncoeff=(order[i]+1)*(order[i]+2)/2;      // psf values
            //cout<<xpos[i]<<" "<<ypos[i]<<" "<<ncoeff<<endl;
            std::valarray<double> coeffs;
            table.column("shapelets").read(coeffs, row); 
            //cout<<" "<<ichip<<" "<<xpos[i]<<" "<<ypos[i]<<" "<<coeffs[shapeStart]<<" "<<coeffs[shapeStart+1]<<" "<<coeffs[shapeStart+2]<<" "<<endl;
            for(int j=0;j<nvar;++j) {
              
              det->setVal(j,coeffs[shapeStart+j]);
            }
            
            chip->addDet(det);
          }
        }
        
      }
      catch (CCfits::FitsException& ) {
        cout<<"Can't open chip: "<<inputFile.str()<<" from exposure "<<exp<<" skipping"<<endl;
        return false;
      }
      
      addChip(ichip,chip);
      
    }

    return true;
  }
  

  std::vector<bool> Exposure::getMissing()
  {
    int nchip_var=ny_chip*nx_chip;
    int ncell=chips.size()*ny_chip*nx_chip;
    std::vector<bool> v(ncell,false);
    std::map<int,Chip*>::const_iterator iter=chips.begin();

    int cur_index=0;    
    for(; iter!=chips.end();++iter) {
      
      std::vector<bool> cv=iter->second->getMissing();
      for(int j=0;j<cv.size();++j) v[cur_index]=cv[j];
    }
    
    return v;
  }

  tmv::Vector<float> Exposure::getVals(std::string type)
  { 
    int nchip_var=ny_chip*nx_chip;
    int nfocal=chips.size()*nchip_var;

    std::map<int,Chip*>::const_iterator iter=chips.begin();

    int ntotvar=iter->second->getCell(0)->getNVal(type);
    //cout<<"Total vars: "<<ntotvar*nfocal<<" chips: "<<chips.size()<<endl;
    tmv::Vector<float> v(ntotvar*nfocal,0.0);
    int cur_index=0;
    int cur_chip=0;
    
    for(; iter!=chips.end();++iter,cur_chip++) {
      //cout<<"Get vals from chip: "<<iter->first<<endl;
      std::vector<float> cv=iter->second->getVals(type);
      //cout<<"Got "<<cv.size()<<endl;
      for(int j=0;j<cv.size();++j) {
        
        // really complicated to match current structure
        // this is probably wrong
        //cout<<"  "<<j<<" "<<nfocal*(j%ntotvar)+cur_chip*nchip_var+j/(ntotvar)<<endl;
        v[nfocal*(j%ntotvar)+cur_chip*nchip_var+j/(ntotvar)]=cv[j];

        // I want to switch to this once I am willing to change all my 
        // plotting routines
        //v[cur_index]=cv[j];
        cur_index++;
      }

    }
    return v;
  }




};
    
