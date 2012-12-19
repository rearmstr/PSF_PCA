#include <iostream>
#include "TMV.h"
#include "TMV_Sym.h"

#include <vector>

#include "myClass.h"
#include "myTypeDef.h"
#include "NR.h"

using namespace std;


/* --------------------------------------------------------------------- */
void initializeTMVmat(DMatrix &mat, double initialVal)
{
    /* initialize a TMV matrix (double) to initial value */
    for (size_t i=0; i<mat.nrows(); i++) {
       for (size_t j=0; j<mat.ncols(); j++) { mat(i,j) = initialVal; }
    }

}

/* --------------------------------------------------------------------- */
void resizeDataMat(c_ControlParam &contParam, c_Data &myData)
{    
     /* resize data matrices accord to input
        meanVec, defocusID, defocusCoeff are used by ROOT only, resizing
        is done at mean subtraction and defocus removal respectively.
        The rest should be resized by all nodes.
        missIndex() is used by Wiberg only and commented out for saving space.

        Since getMat() or getRandStarPSF() might modify nrows, ncols and at
        the same time assigning Xmat etc, ROOT resizes the data matrices in
        these routines, while other nodes resize after the broadcasting of
        control parameters.
      */

    myData.Xmat.resize(contParam.nrows,contParam.ncols);
    initializeTMVmat(myData.Xmat,0.0);
    myData.dataMask.resize(contParam.nrows/contParam.nShapelet,contParam.ncols);
    initializeTMVmat(myData.dataMask,1);
    myData.Nmasked.resize(contParam.ncols, 0);
    // myData.missIndex.resize(contParam.nrows*contParam.ncols);

    if ( contParam.icDefocus == 1 ) {
       myData.zTab.resize(contParam.nrows,contParam.nzTabCol);
       initializeTMVmat(myData.zTab,0.0);
    }
}


/* --------------------------------------------------------------------- */
void subtractMean(c_ControlParam &contParam, c_Data &myData)
{
     int i,j,icount,indx;

     myData.meanVec.resize(contParam.nrows);

     for (i=0; i<contParam.nrows; i++) {           // calculate mean vector
        myData.meanVec.at(i)=0.0;
        indx=i%(contParam.nrows/contParam.nShapelet);      // for dataMask
        icount=0;
        for (j=0; j<contParam.ncols; j++) {
           myData.meanVec.at(i)+=myData.Xmat(i,j)*myData.dataMask(indx,j);
           icount += int (myData.dataMask(indx,j));
        }
        if (icount == 0) {
           cout << endl;
           cout << "Error: missing pixel " << i << endl;
           cout << endl;
           exit(1);
        }
        myData.meanVec.at(i) /= icount;
     }
     for (i=0; i<contParam.nrows; i++) {           // subtract mean from Xmat
        for (j=0; j<contParam.ncols; j++) { myData.Xmat(i,j)-=myData.meanVec.at(i); }
     }

}


/* --------------------------------------------------------------------- */
void checkEmptyCell(c_ControlParam &contParam, c_Data &myData)
{
   /* if a cell is empty, fill it with random [0,1] numbers, Nmask[]-=1, and
      set dataMask to 1 (not masked).
   */
     int i,j,icount,indx;

     int nCellTot = contParam.nrows/contParam.nShapelet;
     for (indx=0; indx<nCellTot; indx++) {
        icount=0;
        for (j=0; j<contParam.ncols; j++) {
           icount += int (myData.dataMask(indx,j));
        }
        if (icount == 0) {
           cout << "\t Warning: missing pixel at " << indx << endl;
           for (int icol=0; icol<contParam.ncols; icol++) {
               myData.dataMask(indx,icol) = 1;
               for (int iShape=0; iShape<contParam.nShapelet; iShape++) {
                  myData.Nmasked.at(icol) -= 1;
                  i=indx+iShape*nCellTot;
                  myData.Xmat(i,icol) = ran01();
               }
           }
        }
     }

}
