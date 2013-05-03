#ifndef MYCLASS_H
#define MYCLASS_H

#ifndef myTypeDef_H
#include "myTypeDef.h"
#endif

#include "Bounds.h"

// typedef tmv::Matrix<double> DMatrix;
// typedef tmv::Vector<double> DVector;

class c_Data {

public:
  c_Data(): Xmat(1,1,0.0), dataMask(1,1,1), zTab(1,1,0.0) { Nmiss=0; };
  
  DMatrix Xmat;              // data matrix nrows x ncols
  DMatrix dataMask;          // 0 -> yes; 1-> no. which data values missing
                             // nrows/nShapelet x ncols; initialized to 1
  std::vector<Bounds<> > bounds;
  
  std::vector<double> meanVec;       // mean data vector [nrows]
  // only define when icMean is set
  
  std::vector<int> masked;           // 0 -> no; 1 -> yes. which vector is masked [ncols]
  // replaced by Nmasked, but kept here for useTMV.cpp
  std::vector<int> Nmasked;          // number of data components masked in each exposure
  // [ncols]; can replace masked? Yes, done.
  std::vector<int> missIndex;        // 0 -> missing; 1 -> non-missing (only used by
  // Wiberg) [nrows*ncols];
  int Nmiss;                         // total number of missing data components

  DMatrix zTab;             // defocus lookup table nrows x nzTabCol
  // only define when icDefocus is set
  std::vector<int> defocusID;        // ID in the defocus lookup table best match
  std::vector<double> defocusCoeff;  // to the exposure (and the coefficient) [ncols]

  // note: meanVec, defocusID, defocusCoeff are used by ROOT only, resizing is done
  //       where necessary. The rest should be resized by all nodes.

};


class c_ControlParam {

  /* if the number of parameters in here is changed, modify N_len in the MPI datatype
     definition in the main program accordingly.
  */

public:

  int iverbose,  // print out debugging comments: 0 -> no output; 1 -> y
    icSVD,     // perform PCA using SVD (=1) or not (=0); not reset by icMissing
    icEM,      // perform PCA using EM (=1) or not (=0)
    icWiberg,  // perform PCA using Wiberg algorithm (=1) or not (=0)
  // this part haven't been updated; need to check missIndex()
    icout,     // output control: 0 -> all vectors in one file (line)
  //                 1 -> each vector a file
    icMissing, // 0 -> no missing data; 1 -> missing some data
  // SNAP stars on grid -> set by hand;
  // stars random -> set to 0 then readjust in getRandStarPSF().
  // It is necessary to BCAST icMissing
    icMean,    // 0 -> no mean subtraction; 1 -> mean subtracted
    icDefocus, // remove defocus pattern before PCA (1) or not (0)
    kEigen,    // number of eigen vectors for reconstruction
    iterMax,   // maxium number of iterations for EM & Wiberg PCA
    iSeed,     // random number seed
    icGetMat,  // select routines to read in data matrix; 0 -> getMat()
  //                                         1 -> getRandStarPSF()
    readDESorSNAP,     // read random stars of DES (=1) or SNAP (=2)
    nrows,     // data matrix is nrows x ncols; derived instead of reading in
    nChips,    // number of CCD chips;
  // for random star runs, it is over-written by chipBound.dat
  // for SNAP rect-grid tests, set this to the number of grid
  // AND mc = nc = 1
    mc,        // number of cells along x on a single chip
    nc,        // number of cells along y on a single chip
    nShapelet, // number of shapelet coefficients; iStartShapelet is
  // hardwired in getRandStarPSF()
    ncols,     // number of exposures; over-written by the number of lines
  // in DESpointingList
    nzTabCol,// defocus table number of column
    skip61;// skip ccd 61


};

class c_Result {

public:

  c_Result(): reconXmat(1,1,0.0), reconErrMat(1,1,0.0), reconErrMeanPix(1,1,0.0),
	      reconErrMeanExp(1,1,0.0), eigenCoeffMat(1,1,0.0), Umat(1,1,0.0),
	      Svec(1,0.0) { };

  double reconErr,        // reconstruction error (avg e of all exposures)
  // for non-missing data components only
    reconErrMiss,    // for missing data components only
    reconErrTot;     // for all data components

  DMatrix reconXmat;        // reconstruction of the data matrix nrows x ncols
  DMatrix reconErrMat;      // reconstruction error matrix nrows x ncols
  DMatrix reconErrMeanPix;  // recon error at each pixel 2 x nrows
  DMatrix reconErrMeanExp;  // recon error for each exposure 2 x ncols
                            // reconErrMean*(0,*) for non-missing data components only
                            // reconErrMean*(1,*) for missing data components only

  DMatrix eigenCoeffMat;    // SVD (U S V) eigen coeff matrix nrows x ncols
  DMatrix Umat;             // SVD U matrix nrows x nrows (could be large)
  DVector Svec;                 // SVD S (diagonal) [nrows]

};

class c_inFileName {

public:

  std::string inputDir,
    controlParamFile,
    chipBoundFile,
    nameBaseSNAPpsf,
    dirBaseDES,
    runID_DES,
    DESexpListFile,
    nameBaseZtab;

};

class c_outFileName {

public:

  std::string outputDir,
    starRecFile,
    gridXY2File,
    dataMatFile,
    dataMaskFile,
    eigenValSVDfile,
    eigenVecSVDfile,
    eigenCoefSVDfile,
    singularSVDfile,
    reconSVDfile,
    reconErrSVDfile,
    reconErrPixSVDfile,
    reconErrExpSVDfile,
    eigenVecEMfile,
    eigenCoefEMfile,
    reconEMfile,
    reconErrEMfile,
    reconErrPixEMfile,
    reconErrExpEMfile,
    defocusIDfile,
    defocusCoeffFile;

};

#endif
