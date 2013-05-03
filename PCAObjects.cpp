#include "PCAObjects.h"
#include <CCfits/CCfits>
#include "TMV.h"
#include <sstream>
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

  static void setPRow(int fitorder, Position pos, const Bounds& bounds, DVectorView prow)
  {
    Assert(int(prow.size()) == (fitorder+1)*(fitorder+2)/2);
    FVector px =
      definePXY(fitorder,pos.getX(),bounds.getXMin(),bounds.getXMax());
    FVector py =
      definePXY(fitorder,pos.getY(),bounds.getYMin(),bounds.getYMax());
    int pq = 0;
    for(int n=0;n<=fitorder;++n) {
      for(int p=n,q=n-p;q<=n;--p,++q) {
        Assert(pq < int(prow.size()));
        prow(pq) = px[p]*py[q];
        ++pq;
      }
    }
    Assert(pq == int(prow.size()));
  }








  // Return the mean values of all the detections in a cell
  std::vector<float> Cell::getMeanVals()
  {
    std::vector<float> v(nvar,0.0);
    if(dets.size()==0) return v;

    for(int i=0;i<dets.size();++i) {
      for(int j=0;j<nvar;++j) {
        v[j]+=dets[i]->getVal(j);
      }
    }


    for(int j=0;j<nvar;++j) {
      v[j]/=dets.size();
      //cout<<"    Cell var "<<j<<" "<<v[j]<<endl;
    }
    return v;

  }

  // Return the median values of all the detections in a cell
  std::vector<float> Cell::getMedianVals()
  {
    std::vector<float> v(nvar);

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
      
  std::vector<float> Cell::getVals(std::string type)
  {
    if(type=="mean") {
      return getMeanVals();
    }
    else if(type=="median") {
      return getMedianVals();
    }
    
  }
  
  // Get the mean values of all the cells in a chip
  // The ordering of the variables are 
  // Ce1l 1 var1..varN, Cell2 var1..varN, Cell3...
  std::vector<float> Chip::getVals(std::string type)
  { 
    std::vector<float> v(nvar*cells.size());
    int cur_index=0;
    for(int i=0;i<cells.size();++i) {
      //cout<<"  Get Mean from cell "<<i<<endl;
      std::vector<float> cv=cells[i].getVals(type);
      for(int j=0;j<nvar;++j) {
        v[cur_index]=cv[j];
        cur_index++;
      }
    }
    return v;
  }


  void Chip::divide(int _nx,int _ny) {
    nx=_nx;
    ny=_ny;
    std::vector<Bounds<float> > vb=bounds.divide(nx,ny);
    for(int i=0;i<vb.size();++i) {
      Cell cell(nvar,vb[i]);
      cells.push_back(cell);
      cbounds.push_back(vb[i]);
    }
  }
  

  void Chip::addDet(Detection *det) {

    // The bounds class has the y as the fast moving coordinate
    int bin_x=static_cast<int>(det->getPos().x/(bounds.getXMax()/nx));
    int bin_y=static_cast<int>(det->getPos().y/(bounds.getYMax()/ny));
    int bin=bin_x*ny+bin_y;

    cells[bin].addDet(det);
  }



  Exposure::Exposure (string _label,int _nchip, int _nvar, double _ra,double _dec,float _airmass):
    label(_label),nchip(_nchip),nvar(_nvar),ra(_ra),dec(_dec),airmass(_airmass),
    nx_chip(-1.),ny_chip(-1.),xmax_chip(-1.),ymax_chip(-1.),shapeStart(3) {}

  bool Exposure::readShapelet(std::string dir,std::string exp) {
    if (exp.empty()) exp=label;
    cout << "Reading exposure " << exp<< endl;
    for(int ichip=1;ichip<=nchip;++ichip) {
      
      Chip *chip=new Chip(nvar,ichip,xmax_chip,ymax_chip);
      chip->divide(nx_chip,ny_chip);

      // check if this chip should be skipped
      std::vector<int>::iterator iter=find(skip.begin(),skip.end(),ichip);
      if(iter!=skip.end()) continue;
      
      std::stringstream inputFile;
      inputFile << dir << "/" << exp << "_";
      if(ichip<10) inputFile <<0;
      
      inputFile << ichip << "_psf.fits";


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
        cout<<"Can't open chip: "<<ichip<<" from exposure "<<exp<<" skipping"<<endl;
        return false;
      }
      
      addChip(ichip,chip);
      
    }

    return true;
  }
  




  tmv::Vector<float> Exposure::getVals(std::string type)
  { 
    int nchip_var=ny_chip*nx_chip;
    int nfocal=chips.size()*nchip_var;
    tmv::Vector<float> v(nvar*nfocal);
    std::map<int,Chip*>::iterator iter=chips.begin();

    int cur_index=0;
    int cur_chip=0;
    for(; iter!=chips.end();++iter,cur_chip++) {
      //cout<<"Get mean from chip: "<<iter->first<<endl;
      std::vector<float> cv=iter->second->getVals(type);
      for(int j=0;j<cv.size();++j) {
        // really complicated to match current structure
        // this is probably wrong
        v[nfocal*(j%nvar)+cur_chip*nchip_var+j/(nvar)]=cv[j];

        //v[cur_index]=cv[j];
        cur_index++;
      }
    }
    return v;
  }




};
    
