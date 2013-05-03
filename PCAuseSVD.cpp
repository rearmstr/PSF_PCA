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
  int nvar=contParam.nrows;
  int nobs=contParam.ncols;
  int npca=contParam.kEigen;
  int jID;
  double rmsMean=0.0;
  RES.reconXmat.setZero();

  if( nvar>nobs) {
    RES.reconXmat=RES.Umat.colRange(0,npca)*
      RES.eigenCoeffMat.rowRange(0,npca);
  } 
  else {
    RES.reconXmat=RES.Umat.rowRange(0,npca)*
      RES.eigenCoeffMat.colRange(0,npca);
  }
  
  cout<<"Reconstruction done"<<endl;
  RES.reconErr=0.0;                          // reconstruction error
  RES.reconErrMeanPix.row(0).setZero();
  RES.reconErrMeanExp.row(0).setZero();


  RES.reconErrMat= myData.Xmat.transpose()-RES.reconXmat;
  rmsMean=RES.reconErrMat.sumElements()/(nobs*nvar);
  RES.reconErrMat=ElemProd(RES.reconErrMat,RES.reconErrMat);
  RES.reconErr=RES.reconErrMat.sumElements()/(nobs*nvar);
  for(int i=0;i<nobs;++i) {
    RES.reconErrMeanExp(0,i)=sqrt( RES.reconErrMat.row(i).sumElements())/double(nvar);
  }
  for(int i=0;i<nvar;++i) {
    RES.reconErrMeanPix(0,i)=sqrt(RES.reconErrMat.col(i).sumElements())/double(nobs);
  }

  RES.reconErr=std::sqrt(std::abs(RES.reconErr - rmsMean*rmsMean));

  cout << "#\t mean of the residual rms = " << rmsMean << endl;
  cout << "#\t rms  of the residual rms = " << RES.reconErr << endl;


  //RES.reconXmat.resize(nvar,nobs);
  //add mean and reverse the rows and column for mat as before
//   if (contParam.icMean == 1) {                     
//      for (size_t i=0; i<nobs; i++) {
//        for (size_t j=0; j<nvar; j++) { RES.reconXmat(j,i)=RES.reconXmat(j,i)+myData.meanVec.at(j); }
//      }
//   }
  



//   for (size_t i=0;i<myData.Xmat.nrows();i++)
//   {
//     for (size_t j=0;j<myData.Xmat.ncols();j++)
//     {
//       RES.reconXmat(i,j)=0.0;
//        // for (size_t k=0;k<myData.Xmat.nrows();k++)
//        for (size_t k=0;k<contParam.kEigen;k++)
//        {
//           RES.reconXmat(i,j) += RES.Umat(i,k) * RES.eigenCoeffMat(k,j);
//        }
//     }
//   }
  
//   RES.reconErr=0.0;                          // reconstruction error
//   for (size_t i=0; i<contParam.nrows; i++) { RES.reconErrMeanPix(0,i) = 0.0; }
//   for (size_t j=0; j<contParam.ncols; j++) { RES.reconErrMeanExp(0,j) = 0.0; }

//   for (size_t i=0; i<contParam.nrows; i++) {
//     for (size_t j=0; j<contParam.ncols; j++) {
//       RES.reconErrMat(i,j) = pow(myData.Xmat(i,j)-RES.reconXmat(i,j), 2);
//       RES.reconErrMeanPix(0,i) += RES.reconErrMat(i,j);              // statistics
//       RES.reconErrMeanExp(0,j) += RES.reconErrMat(i,j);
//       RES.reconErr += RES.reconErrMat(i,j);
//     }
//   }
//   for (size_t i=0; i<contParam.nrows; i++) {
//     RES.reconErrMeanPix(0,i)=sqrt(RES.reconErrMeanPix(0,i)/contParam.ncols);
//   }
//   for (size_t j=0; j<contParam.ncols; j++) {
//     RES.reconErrMeanExp(0,j)=sqrt(RES.reconErrMeanExp(0,j)/contParam.nrows);
//     rmsMean += RES.reconErrMeanExp(0,j);
//   }
//   rmsMean = rmsMean/contParam.ncols;
//   RES.reconErr=sqrt(RES.reconErr/(contParam.nrows*contParam.ncols) - rmsMean*rmsMean);

//   cout << "#\t mean of the residual rms = " << rmsMean << endl;
//   cout << "#\t rms  of the residual rms = " << RES.reconErr << endl;

//   if (contParam.icMean == 1) {                     // add mean back to reconstruction
//      for (size_t i=0; i<contParam.nrows; i++) {
//         for (size_t j=0; j<contParam.ncols; j++) { RES.reconXmat(i,j)+=myData.meanVec.at(i); }
//      }
//   }

//   if (contParam.icDefocus == 1) {                 // add defocus pattern back to recon
//      for (size_t j=0; j<contParam.ncols; j++) {
//        jID=myData.defocusID[j];
//        for (size_t i=0; i<contParam.nrows; i++) {
//          RES.reconXmat(i,j) += myData.defocusCoeff[j]*myData.zTab(i,jID);
//        }
//      }
//   }
}


/* --------------------------------------------------------------------- */
int PCAuseSVD(c_ControlParam &contParam, c_Data &myData, c_outFileName &outName)
{
  fstream myfile;                        // PCA using SVD
  double ScumS=0.0;                      // cummulative signal
  DVector ScumN(contParam.nrows,0.0);    // cummulative noise
  tmv::SymMatrix<double> XcovMat(contParam.nrows);

  c_Result RES;

  int nvar=contParam.nrows;
  int nobs=contParam.ncols;
  int npca=contParam.kEigen;
  RES.reconXmat.resize(nobs,nvar);
  RES.reconErrMat.resize(nobs,nvar);
  RES.reconErrMeanPix.resize(2,contParam.nrows);
  RES.reconErrMeanExp.resize(2,contParam.ncols);

  RES.Svec.resize(nvar);
  DDiagMatrix Svec(nvar);
 
  // assume average has already been subtracted off for each variable
  // Subtract of the average
  // DVector ave(nvar);
  //   for(int i=0;i<nobs; ++i) {
  //     ave+=RES.Umat.row(i);
  //   }
  
  //   for(int i=0;i<var; ++i) {
  //     RES.Umat.col(i)-=ave(i);
  //   } 

  if(nobs > nvar) {
    RES.eigenCoeffMat.resize(nvar,nvar);
    Svec.resize(nvar);
    RES.Svec.resize(nvar);
    RES.Umat.resize(nobs,nvar);
    RES.Umat=myData.Xmat.transpose();
    SV_Decompose(RES.Umat,Svec,RES.eigenCoeffMat,true);
    RES.Svec=Svec.diag();
    ScumN(nvar-1) = RES.Svec(nvar-1);   //cummulative noise
    for (int i=nvar-2; i>=0; i--) { ScumN(i) = ScumN(i+1) + RES.Svec(i); }
  }
  else {
    RES.eigenCoeffMat.resize(nobs,nvar);
    Svec.resize(nobs);
    RES.Svec.resize(nobs);
    RES.Umat.resize(nobs,nobs);

    RES.eigenCoeffMat = myData.Xmat.transpose();

    SV_Decompose(RES.eigenCoeffMat.transpose(),
		 Svec,RES.Umat.transpose());

    // reduce number of principal components
    if(contParam.kEigen > nobs) contParam.kEigen=nobs;

    RES.Svec=Svec.diag();
    ScumN(nobs-1) = RES.Svec(nobs-1);   //cummulative noise
    for (int i=nobs-2; i>=0; i--) { ScumN(i) = ScumN(i+1) + RES.Svec(i); }

  }
  
  outputToFile (RES.eigenCoeffMat.transpose(), outName.eigenVecSVDfile+"_before");
  RES.eigenCoeffMat=Svec*RES.eigenCoeffMat;


//   RES.reconXmat.resize(contParam.nrows,contParam.ncols);
//   RES.eigenCoeffMat.resize(contParam.nrows,contParam.ncols);
//   RES.reconErrMat.resize(contParam.nrows,contParam.ncols);
//   RES.reconErrMeanPix.resize(2,contParam.nrows);
//   RES.reconErrMeanExp.resize(2,contParam.ncols);
//   RES.Umat.resize(contParam.nrows,contParam.nrows);
//   RES.Svec.resize(contParam.nrows);

//   cout << "##### SVD ...\n";

//   XcovMat = myData.Xmat * myData.Xmat.transpose();

//   RES.Umat = XcovMat.svd().getU();           / eigen is the row vector
//   RES.Svec = XcovMat.svd().getS().diag();

//   /cout << "#\t eigen vectors (column) U = " << RES.Umat << endl;
//   / cout << "#\t eigen values S = " << RES.Svec << endl;

//   ScumN(contParam.nrows-1) = RES.Svec(contParam.nrows-1);   / cummulative noise
//   for (int i=contParam.nrows-2; i>=0; i--) { ScumN(i) = ScumN(i+1) + RES.Svec(i); }

  myfile.open(outName.eigenValSVDfile.data(),ios::out);    // output eigen values

  int nrow=nvar;
  if(nobs < nvar) {
    nrow=nobs;
  }
  for (int i=0; i<nrow; i++) {
    ScumS += RES.Svec(i);
    myfile << RES.Svec(i) << "  " << ScumS << "  " << ScumN(i) << endl;
  }

//   RES.eigenCoeffMat = RES.Umat.transpose() * myData.Xmat;   / eigen coefficients
//   /cout << "\t eigen coefficients = " << RES.eigenCoeffMat << endl;

  reconstruction(contParam,myData,RES);
  //cout << "\t reconstruction = "
  //     << RES.reconXmat.subMatrix(0,contParam.nrows,0,5) << endl;

  // outputToFile (RES.Umat*RES.Umat.transpose(), "results/UUT");
  if (contParam.icout == 0) { 
    
    // we need generalize these
    //outputToFile (RES.Umat.transpose(),     outName.eigenVecSVDfile);
    //outputToFile (RES.eigenCoeffMat.transpose(), outName.eigenCoefSVDfile);
    outputToFile (RES.Umat.transpose(),outName.eigenCoefSVDfile);
    outputToFile (RES.eigenCoeffMat.transpose(), outName.eigenVecSVDfile);
    outputToFile (RES.reconXmat.transpose(),outName.reconSVDfile);
    outputToFile (RES.reconErrMat.transpose(),outName.reconErrSVDfile); 
    outputToFile (RES.Svec,outName.singularSVDfile); 
    
    
  }
//   if (contParam.icout == 1) {
//      outputVectorAll(RES.Umat.transpose(),     outName.eigenVecSVDfile+"_");
//      outputToFile (RES.eigenCoeffMat.transpose(), outName.eigenCoefSVDfile+"_");
//      outputVectorAll(RES.reconXmat.transpose(),outName.reconSVDfile+"_");
//     outputVectorAll(RES.reconErrMat.transpose(),outName.reconErrSVDfile+"_"); }

  outputToFile (RES.reconErrMeanPix,outName.reconErrPixSVDfile);
  outputToFile (RES.reconErrMeanExp,outName.reconErrExpSVDfile);


  return 0;
}
