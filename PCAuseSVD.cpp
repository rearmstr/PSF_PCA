#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TMV.h"
#include "TMV_Sym.h"

#include <math.h>   // for sqrt, log, log10 etc

#include <vector>

#include "myClass.h"
#include "myIO.h"
#include "myTypeDef.h"

using namespace std;


/* --------------------------------------------------------------------- */
void reconstruction(c_ControlParam &contParam, c_Data &myData, c_Result &RES)
{
  /* reconstruct the vertors using eigen coefficients and eigen vectors.
     Note that the columns of Umat are the eigen vectors and columns of the
     eigen coefficients matrix has the corresponding eigen coefficients */
  int jID;
  double rmsMean=0.0;

  for (size_t i=0;i<myData.Xmat.nrows();i++)
  {
    for (size_t j=0;j<myData.Xmat.ncols();j++)
    {
       RES.reconXmat(i,j)=0.0;
       // for (size_t k=0;k<myData.Xmat.nrows();k++)
       for (size_t k=0;k<contParam.kEigen;k++)
       {
          RES.reconXmat(i,j) += RES.Umat(i,k) * RES.eigenCoeffMat(k,j);
       }
    }
  }

  RES.reconErr=0.0;                          // reconstruction error
  for (size_t i=0; i<contParam.nrows; i++) { RES.reconErrMeanPix(0,i) = 0.0; }
  for (size_t j=0; j<contParam.ncols; j++) { RES.reconErrMeanExp(0,j) = 0.0; }

  for (size_t i=0; i<contParam.nrows; i++) {
    for (size_t j=0; j<contParam.ncols; j++) {
      RES.reconErrMat(i,j) = pow(myData.Xmat(i,j)-RES.reconXmat(i,j), 2);
      RES.reconErrMeanPix(0,i) += RES.reconErrMat(i,j);              // statistics
      RES.reconErrMeanExp(0,j) += RES.reconErrMat(i,j);
      RES.reconErr += RES.reconErrMat(i,j);
    }
  }
  for (size_t i=0; i<contParam.nrows; i++) {
    RES.reconErrMeanPix(0,i)=sqrt(RES.reconErrMeanPix(0,i)/contParam.ncols);
  }
  for (size_t j=0; j<contParam.ncols; j++) {
    RES.reconErrMeanExp(0,j)=sqrt(RES.reconErrMeanExp(0,j)/contParam.nrows);
    rmsMean += RES.reconErrMeanExp(0,j);
  }
  rmsMean = rmsMean/contParam.ncols;
  RES.reconErr=sqrt(RES.reconErr/(contParam.nrows*contParam.ncols) - rmsMean*rmsMean);

  cout << "#\t mean of the residual rms = " << rmsMean << endl;
  cout << "#\t rms  of the residual rms = " << RES.reconErr << endl;

  if (contParam.icMean == 1) {                     // add mean back to reconstruction
     for (size_t i=0; i<contParam.nrows; i++) {
        for (size_t j=0; j<contParam.ncols; j++) { RES.reconXmat(i,j)+=myData.meanVec.at(i); }
     }
  }

  if (contParam.icDefocus == 1) {                 // add defocus pattern back to recon
     for (size_t j=0; j<contParam.ncols; j++) {
       jID=myData.defocusID[j];
       for (size_t i=0; i<contParam.nrows; i++) {
         RES.reconXmat(i,j) += myData.defocusCoeff[j]*myData.zTab(i,jID);
       }
     }
  }
}


/* --------------------------------------------------------------------- */
int PCAuseSVD(c_ControlParam &contParam, c_Data &myData, c_outFileName &outName)
{
  fstream myfile;                        // PCA using SVD
  double ScumS=0.0;                      // cummulative signal
  DVector ScumN(contParam.nrows,0.0);    // cummulative noise
  tmv::SymMatrix<double> XcovMat(contParam.nrows);

  c_Result RES;

  RES.reconXmat.resize(contParam.nrows,contParam.ncols);
  RES.eigenCoeffMat.resize(contParam.nrows,contParam.ncols);
  RES.reconErrMat.resize(contParam.nrows,contParam.ncols);
  RES.reconErrMeanPix.resize(2,contParam.nrows);
  RES.reconErrMeanExp.resize(2,contParam.ncols);
  RES.Umat.resize(contParam.nrows,contParam.nrows);
  RES.Svec.resize(contParam.nrows);

  cout << "##### SVD ...\n";

  XcovMat = myData.Xmat * myData.Xmat.transpose();

  RES.Umat = XcovMat.svd().getU();           // eigen is the row vector
  RES.Svec = XcovMat.svd().getS().diag();

  // cout << "#\t eigen vectors (column) U = " << RES.Umat << endl;
  // cout << "#\t eigen values S = " << RES.Svec << endl;

  ScumN(contParam.nrows-1) = RES.Svec(contParam.nrows-1);   // cummulative noise
  for (int i=contParam.nrows-2; i>=0; i--) { ScumN(i) = ScumN(i+1) + RES.Svec(i); }

  myfile.open(outName.eigenValSVDfile.data(),ios::out);    // output eigen values
  for (int i=0; i<contParam.nrows; i++) {
    ScumS += RES.Svec(i);
    myfile << RES.Svec(i) << "  " << ScumS << "  " << ScumN(i) << endl;
  }

  RES.eigenCoeffMat = RES.Umat.transpose() * myData.Xmat;   // eigen coefficients
  // cout << "\t eigen coefficients = " << RES.eigenCoeffMat << endl;

  reconstruction(contParam,myData,RES);
  // cout << "\t reconstruction = "
  //      << RES.reconXmat.subMatrix(0,contParam.nrows,0,5) << endl;

  // outputToFile (RES.Umat*RES.Umat.transpose(), "results/UUT");
  if (contParam.icout == 0) { 
     outputToFile (RES.Umat,     outName.eigenVecSVDfile);
     outputToFile (RES.eigenCoeffMat, outName.eigenCoefSVDfile);
     outputToFile (RES.reconXmat,outName.reconSVDfile);
     outputToFile (RES.reconErrMat,outName.reconErrSVDfile); }
  if (contParam.icout == 1) {
     outputVectorAll(RES.Umat,     outName.eigenVecSVDfile+"_");
     outputToFile (RES.eigenCoeffMat, outName.eigenCoefSVDfile+"_");
     outputVectorAll(RES.reconXmat,outName.reconSVDfile+"_");
     outputVectorAll(RES.reconErrMat,outName.reconErrSVDfile+"_"); }

  outputToFile (RES.reconErrMeanPix,outName.reconErrPixSVDfile);
  outputToFile (RES.reconErrMeanExp,outName.reconErrExpSVDfile);

  return 0;
}
