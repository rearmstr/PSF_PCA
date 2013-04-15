#include "mpi.h"
#include <iostream>

#include "TMV.h"
#include "TMV_Sym.h"

#include <math.h>   // for sqrt, log, log10 etc

#include <vector>

#include "myClass.h"
#include "NR.h"          // ran01() and sort()
#include "myIO.h"
#include "myTypeDef.h"
#include "PCAcommon.h"

using namespace std;


/* --------------------------------------------------------------------- */
int PCAuseEM(c_ControlParam &contParam, c_Data &myData, c_outFileName &outName,
             double tol, int myrank, int numprocs)
{
  // PCA using Expectation Maximization (EM)
  // what should be the stopping criterior? || Wnew-W ||/(kEigen*nrows) < tol
  // reconstruction? -- W * XeigenCoeff
  // eigen vector is not normalized?
  // eigen values?
  // In the case of missing data, kEigen+Nmasked[i] < nrows should be satisfied.
  // Otherwise, the problem is under-determined.
  // Warning: initial W should not have big chunks being zero. It might cause
  //          problem after reshuffle and carrying out division of submatrix

  int nrows=contParam.nrows;
  int ncols=contParam.ncols;
  int kEigen=contParam.kEigen;

  const int icWexpand=0;   // 1: expand W to solve least square using QR
                           // 0: smart way to do it using QR (>10 times faster)

  int i,j,iter,i1,j1,icount,indx,missCnt=0;
  double Werr,wiMag;
  cout<<"Neigen "<<kEigen<<endl;
  DMatrix W(nrows,kEigen,0.0),Wnew(nrows,kEigen),XeigenCoeff(kEigen,ncols),
          tempMat(kEigen,kEigen),tempMat1(kEigen,nrows);

  DMatrix Xvec(kEigen+nrows,1),Yvec(nrows,1);
     // DMatrix Wexpand(nrows,kEigen+nrows),Xvec(kEigen+nrows,1),Yvec(nrows,1);
          // tempWexp(kEigen+nrows,kEigen+nrows),tempWexp1(kEigen+nrows,nrows);

  int ROOT=0,st_tag,st_rank,irank,numToRecv,ibound,jID;
  MPI::Status status;
  DMatrix XvecRecv(kEigen+nrows,1);
  // DMatrix bkXmat(nrows,ncols);                   // backup copy of Xmat
             // reconErrMat took its place to save space.

  c_Result RES;        // used by Root only


  if (kEigen >= nrows) {
    if (myrank == ROOT) {
       cout << "Error: PCAuseEM cannot handle full dimension!\n";
       cout << "\t kEigen = " << kEigen << "\t nrows = " << nrows << endl;
    }
    return 1;
  }

  if (myrank == ROOT) {
    cout << "##### EM algorithm ...\n";
    // initialize eigen vectors w_i
    // cout << " ran01 one is " << ran01() << endl;
    for (i=0;i<kEigen;i++) {
      // W(i,i)=1.0;           // avoid large chunks of W being zero
      for (j=0;j<nrows;j++) { W(j,i)=ran01(); }
      W(i,i)=2.0;
      // for (j=0;j<nrows;j++) { W(j,i)=i+j; }
    }
    // cout << " ran01 two is " << ran01() << endl;
    // outputToFile(W,"results/Wmat0");
    ibound = int (ncols/numprocs) * numprocs;      // used to decide the
    // cout << "\t ibound = " << ibound << endl;      // number of processes to receive by ROOT

    RES.reconErrMat.resize(contParam.nrows,contParam.ncols);
    RES.reconErrMat = myData.Xmat;                 // temporary backup of Xmat
  }

  for (iter=0;iter<contParam.iterMax;iter++)
  {
     MPI::COMM_WORLD.Bcast(&W(0,0),nrows*kEigen,MPI::DOUBLE,ROOT);

     // I chose to calculate tempMat1 by all nodes. It can also be done by ROOT only
     // and broadcast.
     // tempMat = W.transpose() * W;
     // cout << "tempMat: " << tempMat << endl;
     // tempMat1 = tempMat.inverse() * W.transpose();
     // cout << "tempMat1: " << tempMat1 << endl;

     if (contParam.icMissing == 0 )     // no missing data only ROOT do the work
     {
          if (myrank == ROOT) {
            tempMat = W.transpose() * W;
            tempMat1 = tempMat.inverse() * W.transpose();
            XeigenCoeff = tempMat1 * myData.Xmat;
          }
     }
     else                // some data are missing jobs are sent to all nodes
     {                   // need to BCAST W; gather Xvec for Xmat(*,i) and XeigenCoeff
        for (i=myrank; i<ncols;i+=numprocs)      // loop over data vector to
        {                           // calculate XeigenCoeff & missing data
          if (myData.Nmasked[i] == 0) {      // ith vector has no missing data, then as usual
            tempMat = W.transpose() * W;
            tempMat1 = tempMat.inverse() * W.transpose();
            XeigenCoeff.subMatrix(0,kEigen,i,i+1) = 
                            tempMat1*myData.Xmat.subMatrix(0,nrows,i,i+1);
            if (myrank != ROOT) {
               Xvec.subMatrix(0,kEigen,0,1) = XeigenCoeff.subMatrix(0,kEigen,i,i+1);
               MPI::COMM_WORLD.Send(&Xvec(0,0),kEigen+nrows,MPI::DOUBLE,ROOT,i);
            }
          }
          else                      // ith vector contains missing data, then
          {                         // solve for XeigenCoeff & missing data together
             /* reshuffle rows of W so that missing data corresponding to the top portion
                of W, and non-missing data to the lower portion of W. The result is
                stored in Wnew. In later part of the code, Wnew is used to store the
                updated W after one iteration.
              */
            if (icWexpand == 0) {
             int irowMiss=0;
             int irowNonMiss=myData.Nmasked[i];
             for (i1=0; i1<nrows; i1++) {
	       // BA: this index is due to the fact that dataMask only has one value per cell
	       // since it is not necesary to to for all shapelet components.
                indx=i1%(contParam.nrows/contParam.nShapelet);      // for dataMask
		// BA this should be changed to bool, not 0.5
                if (myData.dataMask(indx,i) < 0.5) {     // missing; put corresponding row
                                                        // of W at the top part of Wnew
                  Wnew.subMatrix(irowMiss,irowMiss+1,0,kEigen)=W.subMatrix(i1,i1+1,0,kEigen);
                  irowMiss++;
                }
                else {     // non-missing put the W row to the lower portion of Wnew
                  Wnew.subMatrix(irowNonMiss,irowNonMiss+1,0,kEigen)=W.subMatrix(i1,i1+1,0,kEigen);
                  Yvec(irowNonMiss,0)=myData.Xmat(i1,i);
                  irowNonMiss++;
                }
             }
             Xvec.subMatrix(0,kEigen,0,1) = Yvec.subMatrix(myData.Nmasked[i],nrows,0,1) /
                                            Wnew.subMatrix(myData.Nmasked[i],nrows,0,kEigen);
             Xvec.subMatrix(kEigen,kEigen+myData.Nmasked[i],0,1) = 
                                            Wnew.subMatrix(0,myData.Nmasked[i],0,kEigen) *
                                            Xvec.subMatrix(0,kEigen,0,1);
            }

            if (icWexpand == 1) {
             DMatrix Wexpand(nrows,kEigen+nrows);
             Wexpand.subMatrix(0,nrows,0,kEigen)=W;
             icount=myData.Nmasked[i];
             for (i1=0; i1<nrows; i1++)
              {
                indx=i1%(contParam.nrows/contParam.nShapelet);      // for dataMask
                Yvec(i1,0)=myData.Xmat(i1,i)*myData.dataMask(indx,i);   // set data vector to 0 at missing
                for (j1=kEigen; j1<kEigen+myData.Nmasked[i]; j1++)
                  {
                    Wexpand(i1,j1)=0.0;
                   }
                if (myData.dataMask(indx,i) == 0)
                  {
                    Wexpand(i1,kEigen+myData.Nmasked[i]-icount)=-1.0;
                    icount=icount-1;
                   }
               }

             // Now find Xvec by minimizing chi square,
             //      chi^2 = | Wexpand*Xvec - Yvec |^2
             // possible speed up QR, Cholesky or TMV func (x=b/A)
             // SVD is the most stable one, but slower (NR 15.4.16-20)
             // QR might be the fastest
             // my first method (solving the normal equation) is the slowest
             //       (Wexpand^T * Wexpand)*Xvec = Wexpand^T * Yvec
             // if ( 1 == 0 ) {     // normal equation (slowest method)

                // cout << "tempWexp takes up too much memory, commented out\n";

                // tempWexp.subMatrix(0,kEigen+Nmasked[i],0,kEigen+Nmasked[i]) = 
                //       Wexpand.subMatrix(0,nrows,0,kEigen+Nmasked[i]).transpose()
                //       * Wexpand.subMatrix(0,nrows,0,kEigen+Nmasked[i]);
                // tempWexp1.subMatrix(0,kEigen+Nmasked[i],0,nrows) = 
                //       tempWexp.subMatrix(0,kEigen+Nmasked[i],0,kEigen+Nmasked[i]).inverse() 
                //       * Wexpand.subMatrix(0,nrows,0,kEigen+Nmasked[i]).transpose();
                // Xvec.subMatrix(0,kEigen+Nmasked[i],0,1) = 
                //      tempWexp1.subMatrix(0,kEigen+Nmasked[i],0,nrows) * Yvec;
             // }
             // else {             // use TMV (default to QR)
               // Wexpand.subMatrix(0,nrows,0,kEigen+myData.Nmasked[i]).divideUsing(tmv::DivType dt);
                                // dt=tmv::LU, QR, QRP, SV
               Xvec.subMatrix(0,kEigen+myData.Nmasked[i],0,1) = 
                    Yvec / Wexpand.subMatrix(0,nrows,0,kEigen+myData.Nmasked[i]);
	    }


             if (myrank == ROOT) {
               // 1st take care of the result from ROOT itself
               XeigenCoeff.subMatrix(0,kEigen,i,i+1) = Xvec.subMatrix(0,kEigen,0,1);

               icount=0;             // fill in missing elements in data vector
               for (i1=0; i1<nrows; i1++) {
                  indx=i1%(contParam.nrows/contParam.nShapelet);      // for dataMask
                  if (myData.dataMask(indx,i) < 0.5) {
                      myData.Xmat(i1,i)=Xvec(kEigen+icount,0);
                      icount=icount+1;
                    }
                }
             }
             else {
               MPI::COMM_WORLD.Send(&Xvec(0,0),kEigen+nrows,MPI::DOUBLE,ROOT,i);
             }

          }                     // closing the else where masked[i]=1

          if (myrank == ROOT) {            // ROOT recieving from other ranks

              if (i < ibound) { numToRecv = numprocs - 1; }
              else { numToRecv = ncols-ibound - 1; }

              for ( irank=1; irank<numToRecv+1; irank++ ) {     // receive from other ranks
                 st_tag=i+irank;                               // including these w/o missing data
                 st_rank=irank;
                 MPI::COMM_WORLD.Recv(&XvecRecv(0,0),kEigen+nrows,MPI::DOUBLE,
                                      st_rank,st_tag,status);
 
                 XeigenCoeff.subMatrix(0,kEigen,st_tag,st_tag+1)=XvecRecv.subMatrix(0,kEigen,0,1);

                 icount=0;             // fill in missing elements in data vector
                 for (i1=0; i1<nrows; i1++) {
                    indx=i1%(contParam.nrows/contParam.nShapelet);      // for dataMask
                    if (myData.dataMask(indx,st_tag) < 0.5) {
                        myData.Xmat(i1,st_tag)=XvecRecv(kEigen+icount,0);
                        icount=icount+1;
                      }
                 }                   // end of data recovery for loop
              }                     // end of receiving from other ranks
          }                        // end of ROOT recieving

        }                        // end of for loop over data vector
     }                         // closing the else of some data are missing

     MPI::COMM_WORLD.Barrier();

     if (myrank == ROOT) {
       tempMat = XeigenCoeff * XeigenCoeff.transpose();
       Wnew = myData.Xmat * XeigenCoeff.transpose() * tempMat.inverse();

       Werr=0;              // check for convergence |Wnew-W|/(kEigen*nrows) < tol
       for (i=0; i<kEigen; i++)
          {
             wiMag = 0.0;
             for (j=0; j<nrows; j++)
                {
                   wiMag += pow(Wnew(j,i)-W(j,i), 2);
                 }
             Werr += wiMag;
           }
       Werr=sqrt(Werr)/(kEigen*nrows);

       W = Wnew;
     }

     MPI::COMM_WORLD.Bcast(&Werr,1,MPI::DOUBLE,ROOT);
     MPI::COMM_WORLD.Barrier();

     if (Werr <= tol ) { break; }

     if (myrank == ROOT) {
       if ((iter+1) % 10 == 0) {
          cout << "\t \t   iter = " << iter+1 << "\t W error = " << Werr << endl;
          // outputToFile(Wnew,"results/W");
       }
     }

   }                              // end of the iteration loop

   myData.Xmat.resize(1,1);

  if (myrank == ROOT) {

    cout << "\t Warning: data matrix is destroyed\n";
    cout << "#\t number of iteration spent = " << iter << endl;
    // cout << "#\t eigen vectors (column) = " << Wnew << endl;
    // cout << "\t reconstruction from EM = " << (W*XeigenCoeff).subMatrix(0,nrows,0,5) << endl;

    double rms[ncols],rmsMiss[ncols];     // rms residual of each exposures.
                         // duplicate of reconErrMeanExp, for quicksort.
    double rmsMean,rmsMissMean,rmsMed,rmsMissMed,rmsMeanTot;
    int numMaskedExp;    // total number of exposures that have missing data

    vector<int> icnt(nrows,0), jcnt(ncols,0);

    RES.reconXmat.resize(contParam.nrows,contParam.ncols);
    RES.reconErrMeanPix.resize(2,contParam.nrows);
    initializeTMVmat(RES.reconErrMeanPix, 0.0);
    RES.reconErrMeanExp.resize(2,contParam.ncols);
    initializeTMVmat(RES.reconErrMeanExp, 0.0);

    RES.reconXmat = W*XeigenCoeff;
    RES.reconErr=0.0;                          // reconstruction error
    RES.reconErrMiss=0.0;                          // reconstruction error of missing comp

    for (i=0; i<nrows; i++) {
      indx=i%(contParam.nrows/contParam.nShapelet);      // for dataMask
      for (j=0; j<ncols; j++) {
        RES.reconErrMat(i,j) = pow(RES.reconErrMat(i,j)-RES.reconXmat(i,j), 2);
        if (myData.dataMask(indx,j) < 0.5) {              // statistics missed data components
          RES.reconErrMeanPix(1,i) += RES.reconErrMat(i,j);
          RES.reconErrMeanExp(1,j) += RES.reconErrMat(i,j);
          RES.reconErrMiss += RES.reconErrMat(i,j);
          missCnt++;
          icnt[i]++;
          jcnt[j]++;
        }
        else {                                 // statistics non-missing data components
          RES.reconErrMeanPix(0,i) += RES.reconErrMat(i,j);
          RES.reconErrMeanExp(0,j) += RES.reconErrMat(i,j);
          RES.reconErr += RES.reconErrMat(i,j);
        }
      }
    }
    RES.reconErrTot = RES.reconErr + RES.reconErrMiss;
    // RES.reconErrTot = sqrt(RES.reconErrTot/(nrows*ncols));

    RES.reconErr=sqrt(RES.reconErr/(nrows*ncols-missCnt));        // rms of rms of single exposures
    if (missCnt == 0) { RES.reconErrMiss=0; }
    else { RES.reconErrMiss=sqrt(RES.reconErrMiss/missCnt); }

    for (i=0; i<nrows; i++) {                             // per pixel recon err
      RES.reconErrMeanPix(0,i) = sqrt(RES.reconErrMeanPix(0,i)/(ncols-icnt[i]));
      if (icnt[i] == 0) { RES.reconErrMeanPix(1,i)=0; }
      else { RES.reconErrMeanPix(1,i) = sqrt(RES.reconErrMeanPix(1,i)/icnt[i]); }
    }

    // cout << "#\t rms of rms: " << RES.reconErr << " " << RES.reconErrMiss
    //      << " " << RES.reconErrTot << endl;
    // cout << "#\t missCnt = " << missCnt << endl;

    numMaskedExp=0;                                       // mean rms'
    rmsMean=0.0;
    rmsMissMean=0.0;
    rmsMeanTot=0.0;
    RES.reconErr=0.0;               // used for rms of rms of single expousres
    RES.reconErrMiss=0.0;           // similar role as reconErrTot but should calc differently

    for (j=0; j<ncols; j++) {                             // rms of single exposures
      rmsMeanTot += sqrt((RES.reconErrMeanExp(0,j)+RES.reconErrMeanExp(1,j))/nrows);
      rms[j] = sqrt(RES.reconErrMeanExp(0,j)/(nrows-jcnt[j]));
      RES.reconErrMeanExp(0,j) = rms[j];
      if (jcnt[j] == 0) { RES.reconErrMeanExp(1,j)=0; rmsMiss[j]=0; }
      else {
         RES.reconErrMeanExp(1,j) = sqrt(RES.reconErrMeanExp(1,j)/jcnt[j]);
         rmsMiss[j] = RES.reconErrMeanExp(1,j);
      }
      if (myData.Nmasked[j] > 0) numMaskedExp++;
      rmsMean += rms[j];
      rmsMissMean += rmsMiss[j];
      RES.reconErr += pow(rms[j],2);
      RES.reconErrMiss += pow(rmsMiss[j],2);
    }

    rmsMean = rmsMean/ncols;
    RES.reconErr = RES.reconErr/ncols;
    if (numMaskedExp > 0) {
       rmsMissMean = rmsMissMean/numMaskedExp;     // zero otherwise
       RES.reconErrMiss = RES.reconErrMiss/numMaskedExp;
    }
    rmsMeanTot = rmsMeanTot/ncols;

    sort(ncols,rms);                                      // median of the rms'
    sort(ncols,rmsMiss);
    rmsMed=rms[int(ncols/2)];
    rmsMissMed=rmsMiss[int(ncols/2)];

    // cout << "#\t rms of rms: " << RES.reconErr << " " << RES.reconErrMiss
    //      << " " << RES.reconErrTot << endl;
    // cout << "#\t numMaskedExp = " << numMaskedExp << endl;

    RES.reconErr=sqrt(RES.reconErr - rmsMean*rmsMean);
    RES.reconErrMiss=sqrt(RES.reconErrMiss - rmsMissMean*rmsMissMean);
    RES.reconErrTot=sqrt(RES.reconErrTot/(nrows*ncols) - rmsMeanTot*rmsMeanTot);

    cout << "#\t mean of the residual rms   = " << rmsMean << " " << rmsMissMean << endl;
    cout << "#\t median of the residual rms = " << rmsMed << " " << rmsMissMed << endl;
    cout << "#\t rms of the residual rms    = " << RES.reconErr 
         << " " << RES.reconErrMiss << " " << RES.reconErrTot << endl;

    if (contParam.icMean == 1) {                     // add mean back to reconstruction
       for (i=0; i<nrows; i++) {
          for (j=0; j<ncols; j++) { RES.reconXmat(i,j)+=myData.meanVec.at(i); }
       }
    }

    if (contParam.icDefocus == 1) {                 // add defocus pattern back to recon
       for (j=0; j<ncols; j++) {
         jID=myData.defocusID[j];
         for (i=0; i<nrows; i++) {
           RES.reconXmat(i,j) += myData.defocusCoeff[j]*myData.zTab(i,jID);
         }
       }
    }

    // tempMat = XeigenCoeff * XeigenCoeff.transpose();
    // tempMat = W.transpose() * W;          // normalise eigen vectors
    // for (i=0; i<kEigen; i++) {
    //   tempMat(i,i) = 1.0/sqrt(tempMat(i,i));
    //   for (j=0; j<nrows; j++) { W(j,i) = W(j,i) * tempMat(i,i); }
    // }

    if (contParam.icout == 0) {
       outputToFile (W,     outName.eigenVecEMfile);
       outputToFile (XeigenCoeff, outName.eigenCoefEMfile);
       outputToFile (RES.reconXmat,outName.reconEMfile);
       outputToFile (RES.reconErrMat,outName.reconErrEMfile); }
    if (contParam.icout == 1) {
       outputVectorAll(W,     outName.eigenVecEMfile+"_");
       outputToFile (XeigenCoeff, outName.eigenCoefEMfile+"_");
       outputVectorAll(RES.reconXmat,outName.reconEMfile+"_");
       outputVectorAll(RES.reconErrMat,outName.reconErrEMfile+"_"); }

    outputToFile (RES.reconErrMeanPix,outName.reconErrPixEMfile);
    outputToFile (RES.reconErrMeanExp,outName.reconErrExpEMfile);
  }    // closing ROOT

  return 0;
}
