/* this version includes assigning stars at random locations to the grid cells.
   The previous version is saved in MBPbnl as v6_useTMV.cpp
   The idea is to setup a grid on each chip in the focal plane, assign stars to 
   the grid cells, take the average PSF of the stars in a cell as its PSF.
   The data is read in and processed exposure by exposure.
   Every star has exposure number, chip number, x & y in the chip, and PSF.
   Cell assignment gives each star a cell id (i,j) in the chip the it belongs.
      nrows = mc*nc * nChips
      ncols = nExp
   A separate subroutine is written in parallel to getMat() to accomplish this.
*/
#include "mpi.h"
#include <iostream>
#include <string>
#include "TMV.h"
#include "TMV_Sym.h"

#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()
#include <math.h>   // for sqrt, log, log10 etc

// #include <sys/types.h>   // for mkdir
// #include <sys/stat.h>    // for mkdir

#include "myClass.h"
#include "myTypeDef.h"
#include "myIO.h"
#include "NR.h"
#include "initialize.h"
#include "PCAcommon.h"
#include "PCAuseSVD.h"
#include "PCAuseEM.h"
#include "PCAuseWiberg.h"
#include "rmDefocus.h"

using namespace std;

/* Matrix dimmension assumptions:
     Each column of Xmat is a data vector. It is for these vectors that
     the eigens are sought after. The eigens can be derived from Xmat
     or the covariance matrix. I concluded that for our application,
     it is faster to work in covariance matrix space. In this case,
     the resulting U and V matrix are the same (X=USVtranspose).
     Each column vector of U correspondes to an eigen vector. 

   Note that TMV's svd().getV() returns Vtranspose, or one can think
     of TMV using convention of X = U S V.

   PCAuseSVD is implemented
   PCAuseEM  is implemented including missing data (icMissing=1)
   PCAuseWiberg Y -> U VT, each row of Y is a data vector
*/


// int main( int argc, char *argv[] )      // the same as ** w/o []
int main( int argc, char **argv )
{
  int jID,i,j;
  double tol;        // || eigenVecs ||/(kEigen*nrows) < tol
  double SNAPmaskPct;
  // string s,nameBase,fileName;

  c_Data myData;
  c_ControlParam contParam;
  c_inFileName inName;
  c_outFileName outName;

  int numprocs, myrank, namelen, ROOT=0;            // for MPI
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  clock_t startT=clock(),readT,svdT,emT,wT;

  MPI::Init(argc, argv);
  numprocs = MPI::COMM_WORLD.Get_size();
  myrank   = MPI::COMM_WORLD.Get_rank();
  MPI::Get_processor_name(processor_name, namelen);

  MPI::Datatype INT_N;       // define MPI datatype N integers for pass control parameters

  int N_len=20;                       // hardwired; should match c_ControlParam

  int INT_N_len[N_len];
  MPI::Datatype INT_N_type[N_len];
  MPI::Aint INT_N_displace[N_len];
  MPI::Aint intSize;

  intSize=(MPI::Aint) MPI::INT.Get_size();
  INT_N_displace[0]=(MPI::Aint) 0;
  for (i=1; i<N_len; i++) { INT_N_displace[i]=i*intSize; }
  for (i=0; i<N_len; i++) { INT_N_len[i]=1; INT_N_type[i]=MPI::INT; }

  INT_N=MPI::Datatype::Create_struct(N_len,INT_N_len,INT_N_displace,INT_N_type);
  INT_N.Commit();            // end define MPI datatpe

  // int status=mkdir("/data/mzm/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if ( myrank == ROOT ) {               // input parameters and data

    readFileNames(inName,outName,"fileName.par");
    readInputParams(contParam,SNAPmaskPct,tol,inName.controlParamFile.data());

    srand(contParam.iSeed);
    // srand(time(0));    // set initial seed value to system clock
    // testRandNum();

    if (contParam.icGetMat == 0) { getMat(contParam,SNAPmaskPct,myData,inName); }
    if (contParam.icGetMat == 1) { getRandStarPSF(contParam,myData,inName,outName,SNAPmaskPct); }

    if (contParam.icMissing == 1) {
       outputToFile(myData.dataMask,outName.dataMaskFile);
    }
    outputToFile(myData.Xmat,outName.dataMatFile);

  }

  /* BCAST control parameters and data */

  MPI::COMM_WORLD.Bcast(&contParam,1,INT_N,ROOT);

  if (myrank != ROOT) { resizeDataMat(contParam, myData); }

  if (contParam.icDefocus == 1) {
      MPI::COMM_WORLD.Bcast(&myData.Xmat(0,0),
             contParam.nrows*contParam.ncols,MPI::DOUBLE,ROOT);
      MPI::COMM_WORLD.Bcast(&myData.zTab(0,0),
             contParam.nrows*contParam.nzTabCol,MPI::DOUBLE,ROOT);
  }
  int dataLen=contParam.nrows/contParam.nShapelet*contParam.ncols;
  MPI::COMM_WORLD.Bcast(&myData.dataMask(0,0),dataLen,MPI::DOUBLE,ROOT);
  MPI::COMM_WORLD.Bcast(&myData.Nmasked.front(),contParam.ncols,MPI::INT,ROOT);
  // MPI::COMM_WORLD.Bcast(&myData.masked.front(),contParam.ncols,MPI::INT,ROOT);

  /* remove defocus pattern before applying PCA */

  if (contParam.icDefocus == 1) {
    removeDefocus(myrank,numprocs,contParam,myData);    // calc coeff and id
    if ( myrank == ROOT ) {            // remove from Xmat
      for (j=0; j<contParam.ncols; j++) {
        jID=myData.defocusID[j];
        for (i=0; i<contParam.nrows; i++) {
          myData.Xmat(i,j) -= myData.defocusCoeff[j]*myData.zTab(i,jID);
        }
      }
    }
    if ( myrank == ROOT ) {
       outputIntArr(contParam.ncols,myData.defocusID,outName.defocusIDfile);
       outputArr(contParam.ncols,myData.defocusCoeff,outName.defocusCoeffFile);
    }
  }

  /* subtract mean before PCA; perform PCA using SVD */

  if ( myrank == ROOT ) {
    if (contParam.icMean == 1) { subtractMean(contParam,myData); }
                                // mean subtracted from data
    readT=clock();
    cout << "\t read time is " << ((readT-startT)/(double)CLOCKS_PER_SEC) << endl;

    if (contParam.icSVD == 1) { PCAuseSVD(contParam,myData,outName); }

    svdT=clock();
    if (contParam.icSVD == 1) {
      cout << "#\t SVD time is " << ((svdT-readT)/(double)CLOCKS_PER_SEC) << endl;
    }
  }

  MPI::COMM_WORLD.Barrier();

  /* perform EM PCA */

  // rebroadcast Xmat after mean and defocus subtraction
  // MPI::COMM_WORLD.Bcast(&Xmat(0,0),nrows*ncols,MPI::DOUBLE,ROOT);
  MPI::COMM_WORLD.Bcast(&myData.Xmat(0,0),
             contParam.nrows*contParam.ncols,MPI::DOUBLE,ROOT);
  MPI::COMM_WORLD.Bcast(&tol,1,MPI::DOUBLE,ROOT);

  if (contParam.icEM == 1) {
     PCAuseEM(contParam,myData,outName,tol,myrank,numprocs);
     if ( myrank == ROOT ) {
        emT=clock();
        cout << "#\t EM time is " << ((emT-svdT)/(double)CLOCKS_PER_SEC) << endl;
     }
  }

  /* perform PCA using Wiberg's algorithm */

  if ( myrank == ROOT ) {
     if (contParam.icWiberg == 1) {
        PCAuseWiberg(contParam,myData,tol);
        wT=clock();
        cout << "#\t Wiberg time is " << ((wT-emT)/(double)CLOCKS_PER_SEC) << endl;
     }
  }

  MPI::Finalize();

  return 0;
}
