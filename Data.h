#ifndef DATA_H
#define DATA_H

#include "Detection.h"
#include "Utils.h"
#include "UseTMV.h"
#include <string>
using std::string;


namespace PCA {

  //Abstract base class
  // It must know how to return a DVector of results
  class PCData {

    virtual bool addExp(Exposure<float> *exp)=0;
    virtual void calcData(std::string type,const std::vector<float> &params)=0;
    virtual FMatrixView getData()=0;
    virtual bool isMissing(string exp,int var) const=0;
  };




  class Cell {

  public:
    
    Cell(int _nvar=0,float xmin=0, float xmax=0,float ymin=0, float ymax=0): 
      data(_nvar,defaultVal),
      nvar(_nvar),missing(false),
      bounds(xmin,xmax,ymin,ymax),
      clipped(false) {}
    
    Cell(int _nvar,Bounds<float> b): 
      data(_nvar,defaultVal),nvar(_nvar),bounds(b),
      clipped(false),missing(false) {}
    

    // needs explicit assignment operator to use push_back in vector
    // not sure I understand this?
    Cell& operator=(const Cell &a) { 
      nvar=a.nvar;
      data=a.data;
      bounds=a.bounds;
      dets=a.dets;
      clipped=a.clipped;
      missing=a.missing;
      return *this;
    }

    void addDet(Detection<float>* _det) {dets.push_back(_det);}
    int getNDet() const {return dets.size();}
    DVector getData() const {return data;}
    bool isClipped() const {return clipped;}
    void setClipped(bool miss) {clipped=clipped;}
    
    bool isMissing() const {
      if(dets.size()==0 || missing) return true;
      return false;
    }
    void setMissing(bool miss) {missing=miss;}
    
    Position<float> getCenter() const{
      return bounds.center();
    }
    
    void setBounds(Bounds<float> &b) {bounds=b;}
    void setNvar(int var) {nvar=var;data.resize(nvar);}
     
    void calcVals(std::string type,const std::vector<float> &params);

    Detection<float> * getDet(int i) {
      if (i>=0 && i<dets.size()) return dets[i];
      return 0;
    }


  private:
    int nvar;
    DVector data;
    Bounds<float> bounds;
    std::vector<Detection<float>* > dets;
    bool clipped;
    bool missing;
    
  };

      
  // Bin each chip into different bins
  class BinnedData: public PCData {
   public:
    

    class VarIndex {

    public:
      VarIndex(int _nccd,int _nvar,int _nbin):
	nccd(_nccd),nvar(_nvar),nbin(_nbin) {
	ntot=nccd*nbin*nvar;
	ntot_ccd=nbin*nvar;
      }
      
      int operator()(int ccd,int bin,int var) {
	return ccd*ntot_ccd+bin*nvar+var;
      }
      int operator()(int ccd,int var) {
	return ccd*ntot_ccd+var;
      }
      int getFirstVar(int ccd,int bin) {
	return ccd*ntot_ccd+bin*nvar;
      }
      int getCCD(int bin) {
	return ntot/ntot_ccd;
      }
      
      int nccd;
      int nvar;
      int nbin;
      int ntot;
      int ntot_ccd;
    };
      

    BinnedData(int _nx, int _ny,int  _nvar, int _nchips):
      nx(_nx),ny(_ny),nvar(_nvar),nchips(_nchips),
      indx(_nchips,nvar,nx*ny),initialized(false) {}

    // from PCData
    void calcData(std::string type,const std::vector<float> &params);
    bool addExp(Exposure<float> *exp);
    bool isMissing(string exp,int var) const {}
    FMatrixView getData() {return data.view(); }
    FVector getMean() {return mean.view(); }

    
    int getNx() const {return nx;}
    int getNy() const {return ny;}
    int getNCCD() const {return nchips;}
    int getNVar() const {return nvar;}

    Exposure<float>* getExp(int i) { return exp_list[i];}
    int getNExp() { return exp_list.size();}

    std::vector<Bounds<float> > getCellBounds() {return cbounds;}


    typedef std::vector<Cell>::iterator CellIter;

    // Iterate through all the cells on a ccd for an exposure
    CellIter BeginExpCCD(int exp,int ccd) {
      return cells[exp][ccd].begin();
    }
    CellIter EndExpCCD(int exp,int ccd) {
      return cells[exp][ccd].end();
    }


   private:
    
    
    // Each row represents one exposure 
    FMatrix data;
    FVector mean;

    // only allocate size once for matrices
    bool initialized;

    // this holds the number of exposures not flagged as outliers
    int cur_exp;

    // label each variable as missing or not
    std::vector<std::vector<bool> > missing;
    
    // this will have nexp*nccd*nx*ny cells
    std::vector< std::vector< std::vector<Cell> > > cells;
    std::vector<Exposure<float>*> exp_list;

    // store name to location in array 
    std::map<string,int> exp_map; 

    // store ccd label to index in array
    std::map<string,int> ccd_map; 

    // bounds of the cells on a single ccd
    std::vector<Bounds<float> > cbounds;
    
    int nx;
    int ny;
    int nvar;
    int nchips;
    VarIndex indx;
  };


      
  class PsfMethod {
   public:
    virtual void writeFits(string file)=0;
    virtual DVector getData(string exp,int chip, float x,float y,int k)=0;
    virtual void readFile(std::string)=0;
    virtual void solve()=0;

  };

  
  class PcaMethod: public PsfMethod {

  public:
    PcaMethod(BinnedData *_data): data(_data) {}
    
    void solve();
    
    // read output file
    void readFile(std::string) {}
    void writeFits(string file);

    // get the resulting solution using the first k PCs
    DVector getData(string exp,int chip, float x,float y,int k) {}
    
    FMatrixView getU() {return U.view();}
    FMatrixView getVt() {return Vt.view();}
    FDiagMatrixView getS() {return S.view();}
   private:
    BinnedData *data;
    FMatrix U;
    FDiagMatrix S;
    FMatrix Vt;
    FVector mean;
    
    
  };

   class EMPcaMethod: public PsfMethod {

   public:
     EMPcaMethod(int _npc,BinnedData *_data,int _min_iter,int _max_iter,
		 float _tol,double _sigma=0.01):
       min_iter(_min_iter),max_iter(_max_iter),tol(_tol),
       npc(_npc), data(_data) {}
     
     void solve();
     
     // read output file
     void readFile(std::string) {}
     void writeFits(string file);
     
     void setSigma(double _sig) {sigma=_sig;}
     double getSigma() const {return sigma;}

     // get the resulting solution using the first k PCs
     DVector getData(string exp,int chip, float x,float y,int k) {}
     
     FMatrixView getC() {return C.view();}
     FMatrixView getX() {return x.view();}
     
   private:
     int npc;
     // sigma for random generation
     float sigma;
     
     int min_iter;
     int max_iter;
     float tol;

     BinnedData *data;
     FMatrix C;
     FMatrix x;

   };



  // Fit each chip to a polynomial
//   class PolyData: public PCData {
//    public:
//     PolyData(int _order);
    
//     // from PCData
//     void calcData();
//     void addExp(Exposure *exp);
//     isMissing(string exp,int var) {return false;}
//     isMissing(string exp,int ccd,int var) {return false;}
//     DMatrixView getData() {return data.view(); }

//    private:
//     DMatrix data;
//     vector<Exposure*> exp_list;
//     map<string,int> exp_map; // store name to location in array
//     int order;
//     int fitparam;
//   };

  
//   class PcaMethod: public PsfMethod {
    
    
//   };
  
//   class EmPcaMethod: public PsfMethod {
      

//   };

//   class PsfReco {

//   public:
//     PsfReco(PsfData *_data,PsfMethod *method);
    

//   private:
//     PsfData *data;
//     PsfMethod *method;
//   };
}

#endif
