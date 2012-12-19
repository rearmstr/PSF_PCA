#include "mpi.h"
#include <iostream>
#include "TMV.h"
#include "TMV_Sym.h"

#include <math.h>   // for sqrt, log, log10 etc

#include <vector>

#include "myClass.h"
#include "myTypeDef.h"

using namespace std;


/* --------------------------------------------------------------------- */
void calRemoveZresidual(int Iexp, int Jtab, double& coeff, double& rms,
                        c_ControlParam &contParam, c_Data &myData)
{                        // given exposure (Iexp) & defocus pattern (Jtab),
                         // find the best fit coeff & calc residual (rms)
                         // note: Iexp & Jtab start counting from 0
   int i,j;
   double sumEobs2=0.,sumEzTab2=0.,sumEobsEzTab=0.;

   if (1 == 1) {         // minimize coeff
      for (i = 0; i < contParam.nrows; ++i) {
        sumEobs2 += myData.dataMask(i,Iexp)*pow(myData.Xmat(i,Iexp),2);
        sumEobsEzTab += myData.dataMask(i,Iexp) * myData.Xmat(i,Iexp) * 
                        myData.zTab(i,Jtab);
        sumEzTab2 += myData.dataMask(i,Iexp)*pow(myData.zTab(i,Jtab),2);
      }

      coeff = sumEobsEzTab/sumEzTab2;
      if (coeff < 0) { rms=1.E24; }
      else { rms = sqrt(sumEobs2 - pow(sumEobsEzTab,2)/sumEzTab2); }
      }
   else {                // coeff is set to 1
      for (i = 0; i < contParam.nrows; ++i) {
        sumEobsEzTab += myData.dataMask(i,Iexp) *
                        pow(myData.Xmat(i,Iexp)-myData.zTab(i,Jtab),2);
      }
      coeff = 1.0;
      rms = sqrt(sumEobsEzTab);
   }

   if (Iexp == -25) {
      cout << Iexp << "\t" << Jtab << "\t" << coeff << "\t" << rms << endl;
   }

}


/* --------------------------------------------------------------------- */
void removeDefocus(int myrank, int numprocs, c_ControlParam &contParam, c_Data &myData)
{
   int i,j,ROOT=0;
   int st_tag,st_rank,irank,numToRecv,ibound;
   MPI::Status status;

   const int Nsamp=contParam.nzTabCol;      // number of points to sample the defocus table
   int npt[Nsamp], dn, imin, nLeft, iminRecv;
   double c[Nsamp],rms[Nsamp],rmsMin,cMin,cVal,rmsVal,cMinRecv;

   if (myrank == ROOT) {
      cout << "##### removing defocus pattern ...\n";
      ibound = int (contParam.ncols/numprocs) * numprocs;      // used to decide the
                                            // number of processes to receive by ROOT
      myData.defocusID.resize(contParam.ncols);
      myData.defocusCoeff.resize(contParam.ncols);
   }

   for (i=myrank; i<contParam.ncols; i+=numprocs)     // loop over data vector to
   {                                        // subtract defocus pattern
      npt[0]=0;
      npt[Nsamp-1]=contParam.nzTabCol-1;

      continueIter:

      nLeft=npt[Nsamp-1]-npt[0]+1;
      dn=int(nLeft/Nsamp);
      for (j=1;j<Nsamp-1;++j) {
         npt[j]=npt[0]+j*dn;
        // if (myrank == ROOT) { cout << j << "\t" << npt[j] << "\t" << i <<
        //                                    "\t" << dn << "\t" << npt[0] << endl; }
      }

      rmsMin=1.0E24;
      imin=npt[0];
      for (j=0;j<Nsamp;++j) {
         calRemoveZresidual(i,npt[j],c[j],rms[j],contParam,myData);
         if (rms[j] < rmsMin) { rmsMin=rms[j]; imin=j; }
         if (i == -25) {
            cout << i << "\t" << j << "\t" << imin << "\t" << rms[j] << "\t" << rmsMin << endl;
         }
      }

      // if (i == -1) {
      //   cout << i << "\t" << imin << "\t" << rmsMin << endl;
      // }
      if (i == -25) {
         if (myrank == ROOT) { cout << endl; }
      }

      if (imin > 0) { npt[0]=npt[imin-1]; }
      if (imin < Nsamp-1) { npt[Nsamp-1]=npt[imin+1]; }
      nLeft=npt[Nsamp-1]-npt[0]+1;

      if (nLeft > 2*Nsamp) { goto continueIter; }
      else {
         rmsMin=1.0E24;
         imin=npt[0];
         cMin=0.0;
         for (j=0;j<nLeft;++j) {
            calRemoveZresidual(i,npt[0]+j,cVal,rmsVal,contParam,myData);
            if (rmsVal < rmsMin) { rmsMin=rmsVal; imin=npt[0]+j; cMin=cVal; }
            if (i == -25) {
               cout << i << "\t" << j << "\t" << imin << "\t" << rmsVal << "\t" << rmsMin
                    << "\t" << cVal << "\t" << cMin << endl;
            }
         }
      }

      // cout << i << "\t" << imin << "\t" << cMin << endl;

      // MPI_SEND cMin & imin for calc residual Xmat (pass small amount of data) ...
      if (myrank == ROOT) {  
         myData.defocusID[i]=imin;     // 1st take care of the result from ROOT itself
         myData.defocusCoeff[i]=cMin;
               // 2nd deal with the results from other nodes

         if (i < ibound) { numToRecv = numprocs - 1; }
         else { numToRecv = contParam.ncols-ibound-1; }

         for ( irank=1; irank<numToRecv+1; irank++ ) {   // receive from all the other ranks
            st_tag=i+irank;
            st_rank=irank;
            MPI::COMM_WORLD.Recv(&cMinRecv,1,MPI::DOUBLE,st_rank,st_tag,status);
            MPI::COMM_WORLD.Recv(&iminRecv,1,MPI::INT,st_rank,st_tag+contParam.ncols,status);
            myData.defocusCoeff[st_tag]=cMinRecv;
            myData.defocusID[st_tag]=iminRecv;
         }
      }
      else {
         MPI::COMM_WORLD.Send(&cMin,1,MPI::DOUBLE,ROOT,i);
         MPI::COMM_WORLD.Send(&imin,1,MPI::INT,ROOT,i+contParam.ncols);
      }
   }
   
}
