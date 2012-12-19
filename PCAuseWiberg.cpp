#include <iostream>
#include "TMV.h"
#include "TMV_Sym.h"

#include <math.h>   // for sqrt, log, log10 etc

#include "myClass.h"
#include "myTypeDef.h"
#include "myIO.h"

using namespace std;


/* --------------------------------------------------------------------- */
int PCAuseWiberg(c_ControlParam &contParam, c_Data &myData, double tol)
{
  // PCA using Wiberg's algorithm
  //    Y -> U VT, each row of Y is a data vector, differs from EM and SVD
  // steps:
  //   1) initialize tilda_v, i.e. V and vec_m (columns of V are eigen vectors)
  //   2) calculate F from tilda_v
  //   3) hat_u = inverse(F^T F)F^T (y-m)
  //   4) calculate tilda_G from hat_u
  //   5) d(tilda_v) from (Q_F tilda_G)^+ Q_F (y-m)
  //      here Q_F = I - F inverse(F^T F)F^T and
  //           (Q_F tilda_G)^+ is T tilda_D^(-1) S^T w/
  //           SDT^T as the SVD decomposition of Q_F tilda_G and
  //           tilda_D^(-1) = diag[1/d1, ..., 1/dq, 0, ..., 0] there are
  //           (n-kEigen)(kEigen+1) nonzero values.
  //   6) tilda_v + d(tilda_v) --> tilda_v and goto step 2) if not yet converge
  //
  // Not tested when there is missing data
  //
  int i,j,k,ncount,iMiss,indx,irow0,icol0,iter;
  int nrows=contParam.nrows;
  int ncols=contParam.ncols;
  int kEigen=contParam.kEigen;
  int m=ncols, n=nrows,                  // Y is mxn matrix
      mn=nrows*ncols-myData.Nmiss,
      mr=ncols*kEigen,
      nr=nrows*kEigen, nr1=nr+nrows;
  DMatrix F(mn,mr,0.0),G(mn,nr1,0.0),indMat_mn(mn,mn,0.0),
          vVec(nr1,1),mVec(mn,1),uVec(mr,1),yVec(mn,1),
          dvVec(nr1,1),QF(mn,mn),QFG(mn,nr1),FTFinvFT(mr,mn),tempMat(mr,mr),
          S(mn,nr1),D(nr1,nr1),T(nr1,nr1),QFGdagger(nr1,mn),reconY(m,n);
  DMatrix V0mat(nrows,kEigen,0.0),m0Vec(nrows,1,0.0),Umat(ncols,kEigen),
          Vmat(nrows,kEigen,0.0),dV(nrows,kEigen);
  DVector Dvec(nr1);

  float Verr,viMag,reconErr;

  cout << "##### Wiberg's algorithm ...\n";
  if ( contParam.icMissing == 1 ) {
    cout << "\t not tested for missing data case\n";
    return 0;
  }

  for (i=0; i<mn; i++) { indMat_mn(i,i)=1.0; }

  cout << "\t 1) initialize eigen vectors, and mean vector\n";

  for (i=0;i<kEigen;i++) { V0mat(i,i)=1.0; }       // eigen vectors
  V0mat(1,0)=0.5;
  V0mat(2,0)=0.2;
  V0mat(0,1)=0.1;
  V0mat(1,1)=0.3;
  V0mat(2,1)=-1.25;

  iMiss=0;
  for (i=0;i<nrows;i++) {
    ncount=0;
    for (j=0; j<ncols; j++) {
      indx=j*nrows+i;                     // mean vec using existing data
      if (myData.missIndex[indx] == 0) { iMiss++; }
      else {
        yVec(indx-iMiss,0)=myData.Xmat(i,j);               // yVec
        ncount+=1;
        m0Vec(i,0)+=myData.Xmat(i,j);
      }
    }
    m0Vec(i,0)=m0Vec(i,0)/ncount;
  }
  cout << "\t iMiss = " << iMiss << endl;

  for (iter=0;iter<contParam.iterMax;iter++)       // iteration starts here
  {
  
    ncount=0;
    for (j=0; j<m; j++) {                   // mVec
      for (i=0; i<n; i++) {
        indx=j*n+i;
        if (myData.missIndex[indx] == 0) { ncount++; }
        else { mVec(indx-ncount,0)=m0Vec(i,0); }
      }
    }

    for (i=0; i<n; i++) {                  // tilda_v
      for (j=0; j<kEigen; j++) {
        indx=i*kEigen+j + i;                   // extra i is counting mean
        vVec(indx,0)=V0mat(i,j);
      }
      vVec(indx+1,0)=m0Vec(i,0);
    }

    cout << "\t 2) construct F then u\n";
    cout << "\t \t F ... \n";

    ncount=0;
    for (i=0; i<m; i++) {                  // construct F from V0mat
      irow0=i*n;                           // F doesn't need to be initialized
      icol0=i*kEigen;                      // every iteration. 0 elements never
      for (j=0; j<n; j++) {                // been touched.
        indx=irow0+j;
        if (myData.missIndex[indx] == 0) { ncount++; }
        else {
        for (k=0; k<kEigen; k++) { F(indx-ncount,icol0+k)=V0mat(j,k); }
        }
      }
    }

    cout << "\t \t u ... \n";

    tempMat=F.transpose() * F;            // uVec
    FTFinvFT=tempMat.inverse() * F.transpose();
    uVec=FTFinvFT * (yVec - mVec);
    outputToFile (uVec, "results/uVec");

    cout << "\t 3) calculate G and dv\n";
    // outputToFile (F, "results/F");

    ncount=0;
    for (i=0; i<m; i++) {                  // construct G from uVec
      for (j=0; j<n; j++) {                // G doesn't need to be initialized
        irow0=i*n + j;                     // every iteration. 0 elements never
        icol0=j*(kEigen+1);                // been touched
        if (myData.missIndex[irow0] == 0) { ncount++; }
        else {
          for (k=0; k<kEigen; k++) {
            indx=i*kEigen+k;
            G(irow0-ncount,icol0+k)=uVec(indx,0);
          }
          G(irow0-ncount,icol0+kEigen)=1.0;
        }
      }
    }
    outputToFile (G, "results/G");

    cout << "\t\t got here 1\n";

    QF = - F * FTFinvFT;
    for (i=0; i<mn; i++) { QF(i,i)+=1.0; }

    cout << "\t\t got here 2\n";

    QFG = QF * G;
    S = QFG.svd().getU();
    D = QFG.svd().getS();
    // Dvec = QFG.svd().getS().diag();
    T = QFG.svd().getV();
    outputToFile (D, "results/D");

    cout << "\t\t got here 3\n";

    for (i=0; i<(n-kEigen)*(kEigen+1); i++) { D(i,i) = 1.0/D(i,i); }
    QFGdagger = T.transpose() * D * S.transpose();

    dvVec=QFGdagger * QF * (yVec - mVec);

    cout << "\t 4) update v and test convergence\n";
    // for (i=0; i<nr1; i++) { dvVec(i,0) *= 0.5; }
    vVec = vVec + dvVec;

    cout << "\t \t construct updated eigen vectors and mean vectors\n";
    for (i=0; i<n; i++) {                  // V0mat and m0Vec
      for (j=0; j<kEigen; j++) {
        indx=i*(kEigen+1)+j;               // reverse of tilda_v construction
        Vmat(i,j) = vVec(indx,0);
      }
      m0Vec(i,0) = vVec(indx+1,0);
    }

    // test convergence and iterate if failed
     dV = Vmat - V0mat;               // check for convergence |dV|^2 < tol
     Verr=0;
     for (i=0; i<kEigen; i++) {
       viMag = 0.0;
       for (j=0; j<nrows; j++) { viMag += pow(dV(j,i), 2); }
       Verr += viMag;
       // cout << i << "\t" << j << "\t" << Verr << endl;
     }

     if (Verr <= tol ) { break; }
     V0mat=Vmat;

  }                                   // end of iteration

  cout << "#\t number of iteration spent = " << iter << endl;
  cout << "#\t error on the eigen vector is " << Verr << endl;
  outputToFile (Vmat, "results/eigenVecWiberg");
  outputToFile (V0mat, "results/V0mat");
  outputToFile (m0Vec, "results/m0Vec");

  for (i=0; i<m; i++) {                  // Umat for reconstruction
    for (j=0; j<kEigen; j++) {
      indx=i*kEigen+j;
      Umat(i,j) = uVec(indx,0);
    }
  }

  reconY = Umat*V0mat.transpose();
  // outputToFile (reconY.transpose(), "results/reconY");  // compiling error

  reconY = reconY - myData.Xmat.transpose();        // reconstruction error
  reconErr=0.0;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      reconErr += pow(reconY(i,j), 2);
    }
  }
  cout << "#\t reconstruction error = " << reconErr << endl;

  return 0;
}
