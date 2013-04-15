#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TMV.h"
#include "TMV_Sym.h"
#include <CCfits/CCfits>

#include <math.h>   // for sqrt, log, log10 etc

#include <vector>

#include <ctype.h>

#include "myClass.h"
#include "myTypeDef.h"
#include "NR.h"
#include "myIO.h"
#include "PCAcommon.h"

using namespace std;


/* --------------------------------------------------------------------- */
void getMat(c_ControlParam &contParam, double& SNAPmaskPct, c_Data &myData,
            c_inFileName &inName)
{
  int i,j,icRead=2,icount,indx,i0;  // 1: read from file; 0: generate here.
  string s,nameBase,fileName;       // 2: read vector by vector; 3: X=VUT
  double tempVal;

  int nrows=contParam.nrows;
  int ncols=contParam.ncols;
  int kEigen=contParam.kEigen;

  resizeDataMat(contParam, myData);

  // for (int i=0; i<nrows*ncols; i++) { missIndex[i]=1; }   // initialize

  if ( icRead == 1 )                       // read Xmat from a file
    {
       fileName="data/MM.dat";
       // fileName="data/randXY.dat";
       // fileName="data/4x10Mat.dat";
       inputFromFile(fileName,contParam,myData);
       // cout << "data: " << myData.Xmat << endl;
    }

  if ( icRead == 0 )                      // generate random data set
    {
       double pi,a,b,c,r,theta,phi;

       pi=4.0*atan(1.0);
       a=100;
       b=10;
       c=10;

       for (int i=0; i<ncols; i++)
          {
             // Xmat(0,i)=randg(2.0);                     // randXY
             // Xmat(1,i)=2.0*Xmat(0,i) + 1.0*randg(0.1);
             //   Xmat(0,i)=randg(2.0);                     // y=2x, z=sin(x+y)
             //   Xmat(1,i)=2.0*Xmat(0,i) + 1.0*randg(0.1);
             //   Xmat(2,i)=1.0*(Xmat(0,i) + 2.0*Xmat(1,i)) + 1.0*randg(0.1);
             theta=pi*(ran01()-0.5);                      // ellipsoid
             phi=2.0*pi*ran01();
             r=ran01();
             myData.Xmat(0,i)=a*r*sin(theta)*cos(phi);
             myData.Xmat(1,i)=b*r*sin(theta)*sin(phi);
             myData.Xmat(2,i)=c*r*cos(theta);
           }
    }

  if ( icRead == 3 )                      // Xmat=V * U_transpose + noise
    {                                     // V is nrows x kEigen
       // int kEigen=3;
       DMatrix Umat(ncols,kEigen),Vmat(nrows,kEigen);
       for (int k=0; k<kEigen; k++)
          {
             for (int i=0; i<nrows; i++) { Vmat(i,k)=randg(1.0); }
             for (int j=0; j<ncols; j++) { Umat(j,k)=randg(1.0); }
           }
       myData.Xmat = Vmat * Umat.transpose();

       for (int i=0; i<nrows; i++) {      // add noise
         for (int j=0; j<ncols; j++) {
           myData.Xmat(i,j) += randg(0.0025);
         }
       }
       outputToFile(myData.Xmat,"results/Xmat");
    }

  if ( icRead == 2 )                       // read Xmat a vector at a time
    {
       // nameBase="data/QSO/500x500/spec";
       // nameBase="data/QSO/164x2000/spec";
       // nameBase="data/SNAP/brightStarsNoNoise/e12_";
       for (i=0; i<ncols; i++)
         {
            s=convertInt(i+1);
            fileName=inName.nameBaseSNAPpsf+s;
            readVector(i,fileName,myData.Xmat,nrows);
         }
    }

  if ( contParam.icDefocus == 1 )                    // read defocus pattern lookup table
    {
       // nameBase="data/SNAP/defocusTab_1k/e12_";
       for (i=0; i<contParam.nzTabCol; i++)
         {
            s=convertInt(i+1);
            fileName=inName.nameBaseZtab+s;
            readVectorDefocus(i,fileName,myData.zTab,nrows);
         }
    }

  if ( 0 == 1 )                       // QSO spectra mask
  {
    // for (j=0; j<100; j++) {           // cut 100 spectra by 10% at large lambda
    // for (i=1800; i<nrows-1; i++) {
    for (j=0; j<ncols; j++) {       // cut all spectra by 10% at random lambda
      i0 = int (ran01()*1780);
      for (i=i0; i<i0+200; i++) {
         // indx=j*nrows+i;
         myData.dataMask(i,j)=0;                // 0 means missing
         // myData.missIndex[indx]=0;
         myData.Nmiss++;
      }
      // myData.masked[j]=1;
      myData.Nmasked[j]=200;
      cout << "\t \t" << j << "\t" << i0 << endl;
   }
  // myData.missIndex[0]=0;
  // myData.dataMask(0,0)=0;
  // myData.masked[0]=1;
  // myData.Nmasked[0]=1;
  }

  if ( contParam.icMissing == 1 )                       // SNAP PSF mask
  {
    for (j=0; j<ncols; j++) {       // cut all PSF by 10% at random location
      myData.Nmasked[j]=0;
      for (i=0; i<nrows/2; i++) {
        tempVal=ran01();
        // if ( tempVal < 0.01 )               // mask if 0-1 random number < 0.1
        if ( tempVal < 0.01*SNAPmaskPct )      // mask if 0-1 random number < SNAPmaskPct/100
        {
          // indx=j*nrows+i;                   // e1
          myData.dataMask(i,j)=0;                // 0 means missing
          // myData.missIndex[indx]=0;
          myData.Nmiss++;
          // myData.masked[j]=1;
          myData.Nmasked[j]+=1;

          // indx=j*nrows+i+nrows/2;                   // e2
          // myData.missIndex[indx]=0;
          myData.Nmiss++;
          myData.Nmasked[j]+=1;
        }
      }
   }
  }
}

/* --------------------------------------------------------------------- */
void read_des_exp(vector< vector<double> > &starPSF, vector<int> &starChip,
                  vector<double> & starX, vector<double> & starY,
                  vector<string> exposureList, int iExp, int nChips,
                  c_inFileName &inName,c_ControlParam &param)
{   /* read a chip at a time, then use insert to concatenate the vectors
            icReadFromHDFS: =1 data files in hadoop file system; =0 otherwise
            DESorBCS: =1 DES; =2 BCS
     */

       int iChip,icReadFromHDFS=0,DESorBCS=2;

       string dirBase,sRunID,sExp,sChip,sFilter,dir,pointing,psf_fileName,fileBase,
              cmdStr,tmpFile;
       const char* inputFile;

       // sRunID="wlse0010t";
       sRunID=inName.runID_DES;
       sFilter="i";

       // dirBase="/astro/tutti1/DES/wlbnl/";
       dirBase=inName.dirBaseDES;

       sExp=convertInt(iExp+1);
       pointing = exposureList.at(iExp);

       if (DESorBCS == 1) {          // read DES
          dirBase+=sRunID+"/";                  // or use dirBase.append(sRunID);
          dir=dirBase+pointing+"/";
          if (! FileExists(dir)) {          // disabled after adding the option of HDFS
              cout << "Error: exposure " << pointing << " does not exist!\n";
              exit(1);
           }
          fileBase=dir + sRunID + "-" + pointing + "-";
       }
       else {                        // read BCS
          fileBase=dirBase + pointing + "_";
       }

       // cout << "\t exposure " << iExp << ": " << dir << endl;
       // cout << fileBase << endl;

       int hduNum = 1;                     // parameters for PSF fits file
       int nCoeff, i;

       std::string xCol="x";
       std::string yCol="y";
       std::string orderCol="psf_order";
       std::string psfCol="shapelets";

       std::string psfFlagCol="psf_flags";

       for (iChip=0; iChip<nChips; iChip++) {
          // read tempStarChip, tempStarPSF for a chip
	 
	 
	 if(param.skip61 && iChip==60) 
	   sChip=convertInt(iChip+2);
	 else 
	   sChip=convertInt(iChip+1);
	 
          if (iChip < 9) { sChip="0"+sChip; }
          if (DESorBCS == 1) {                 // DES
              psf_fileName = fileBase + sChip + "-psf.fits";
          }
          else {                               // BCS
              psf_fileName = fileBase + sChip + "_psf.fits";
          }

          if (icReadFromHDFS != 1) {    // psf fits files are accessible by NFS or on local disk
             inputFile = psf_fileName.data();        // or use fileName.c_str()
          }
          else {                        // psf fits files are in HDFS (cp to local disk then read)
             tmpFile="/tmp/mzm_des_psf.fits";
             inputFile = tmpFile.data();
             cmdStr="hdfsCopy.scr " + psf_fileName + " " + tmpFile;
             system(cmdStr.c_str());
          }

          // cout << "\t iChip = " << iChip << "  " << inputFile << "\n";

          // std::vector< std::vector<double> > psf;     // replace w/ Xmat later

	  cout << "Reading from " << inputFile<< endl;
          std::auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(inputFile,CCfits::Read));

          CCfits::ExtHDU& table = pInfile->extension(hduNum);


          long nTabRows=table.rows();

          long start=1;
          long end=nTabRows;

          std::vector<int> psf_flags;
          table.column(psfFlagCol).read(psf_flags, start, end);

          std::vector<double> x,xtemp;
          std::vector<double> y,ytemp;
          table.column(xCol).read(xtemp, start, end);
          table.column(yCol).read(ytemp, start, end);

          std::vector<long> order;       // shapelet order
          table.column(orderCol).read(order, start, end);

          vector<int> tempStarChip;
          vector< vector<double> > tempStarPSF;
          tempStarPSF.resize(nTabRows);

          int icount=0;
          for (i=0; i<nTabRows; i++) {
             if (psf_flags.at(i) == 0) {          // pass psf flag test

                int row=i+1;

                tempStarChip.push_back(iChip);           // ichip
                x.push_back(xtemp.at(i));                // x & y
                y.push_back(ytemp.at(i));

                nCoeff=(order[i]+1)*(order[i]+2)/2;      // psf values
                std::valarray<double> coeffs;
                table.column(psfCol).read(coeffs, row);
                for (int j=0; j<nCoeff; ++j) tempStarPSF.at(icount).push_back(coeffs[j]);

                icount++;
             }
          }

          // add to the exposure
          starChip.insert(starChip.end(), tempStarChip.begin(), tempStarChip.end());
          starPSF.insert(starPSF.end(), tempStarPSF.begin(), tempStarPSF.begin()+icount);
          starX.insert(starX.end(), x.begin(), x.end());
          starY.insert(starY.end(), y.begin(), y.end());

          if (icReadFromHDFS == 1) {       // remove temp psf file on local disk
             cmdStr="rm " + tmpFile;
             system(cmdStr.c_str());
          }
       }         // end of iChip

}


/* --------------------------------------------------------------------- */
void getRandStarPSF(c_ControlParam &contParam, c_Data &myData,
                    c_inFileName &inName, c_outFileName &outName,
                    double &validateSetPct)
{  /* In parallel to getMat(), this routine read in stars located randomly
      on the focal plane and assign their PSF to the data matrix.
      The idea is to setup a grid on each chip in the focal plane, assign stars
      to the grid cells, take the average PSF of the stars in a cell as its PSF.
      The data is read in and processed exposure by exposure.
      Every star has exposure number, chip number, x & y in the chip, and PSF.
      Cell assignment gives each star a cell id (i,j) in the chip it belongs.
         nrows = mc*nc * nChips * nShapelet
         ncols = nExp
      starRecFile is the output of star information:
         exposureID, chipID, cell(i,j), x, y, e1, e2
         These are used to calculate the residual of the reconstruction. 
      if a cell is empty, fill it with random [0,1] numbers, Nmask[]-=1, and
      set dataMask to 1 (not masked).
      11-25-2011 add ValidateSetPct: the percentage of data that is used as
                 validation sample. 0 or negative means no validation set.
   */

        /* +++++++++++++++++++++++++ CCD chips +++++++++++++++++++++++++ */
        int mc=contParam.mc, nc=contParam.nc, npt2=4;
        int iChip,im,in,istar,i,irow,icount,iExp,iShape,nCellTot,indx;

        // const char* chipBoundFile="chipBound.dat";

        vector<double> rowVec1d(npt2,0.0);
        vector< vector<double> > rowVec2d(nc,rowVec1d);
        vector< vector< vector<double> > > rowVec3d(mc,rowVec2d);
        vector< vector<double> > chipBound;                 // (nChips,4)
        vector< vector< vector< vector<double> > > > cell;  // all cells on focal plane
                                                            // (nChips,mc,nc,4)

        double dx, dy, xcmin, ycmin;

        /* +++++++++++++++++++++ Stars in an exposure +++++++++++++++++++ */
        int Nstar;                           // number of stars in one exposure
        int iStartShapelet;                  // starting index of the shapelet coeff (0,1...)
        iStartShapelet=3;    // DES
        //if (contParam.readDESorSNAP == 2 ) iStartShapelet=5;    // SNAP
        vector< vector<int> > grid2D(mc, vector<int>(nc,0));   // temp 2D vector (mc,nc)
        vector< vector< vector<int> > > starCount; // (nChips,mc,nc) star count in a cell

        vector<string> exposureList;        // exposure (pointing) name list

        string nameBase, s, fileName;
        const char* starPSFfile;

        fstream outFile;
        // char starRecFile[] = outName.starRecFile.data();
        // char gridXY2File[] = outName.gridXY2File.data();

        cout << "##### reading data files ...\n";
        cout << "\t collecting chip info ...\n";

        chipBound = readIn2dData(inName.chipBoundFile.data());

        contParam.nChips=chipBound.size();
        nCellTot=contParam.nChips*contParam.mc*contParam.nc;
        contParam.nrows=nCellTot*contParam.nShapelet;           // reset nrows

        cout << "\t nChips = " << contParam.nChips << endl;
        cout << "\t nShapelet = " << contParam.nShapelet << endl;
        cout << "\t data vector dimension (nrows) is reset to " << contParam.nrows << endl;

        cell.resize( contParam.nChips, rowVec3d );

        if (1 == 0) {         // output chip bounds
          for (long row=0; row<(long)chipBound.size(); ++row) {
              for (long col=0; col<(long)chipBound.at(row).size(); ++col) {
                 cout << chipBound.at(row).at(col) << "\t";
              }
              cout << "\n";
          }
        }

        // from chip bounds and mc and nc to generate cell coordinates

        /* coordinates of the lower left and upper right corner of each cell
           are kept. This is to include the general case where the cells
           are not distributed as a regular grid.
           However, the way I setup the cells and the placement of a star
           in a cell is NOT general. Regular grid cells are assumed.
        */

        cout << "\t generating cells in chips ...\n";
	
        for (i=0; i<contParam.nChips; i++) {
            // cout << "i = " << i << endl;
            xcmin=chipBound.at(i).at(0);
            ycmin=chipBound.at(i).at(1);
            dx=(chipBound.at(i).at(2)-xcmin)/mc;
            dy=(chipBound.at(i).at(3)-ycmin)/nc;
            for (im=0; im<mc; im++) {
                for (in=0; in<nc; in++) {
                    cell.at(i).at(im).at(in).at(0)=xcmin+im*dx;
                    cell.at(i).at(im).at(in).at(1)=ycmin+in*dy;
                    cell.at(i).at(im).at(in).at(2)=xcmin+(im+1)*dx;
                    cell.at(i).at(im).at(in).at(3)=ycmin+(in+1)*dy;
		    myData.bounds.push_back(Bounds<>(xcmin+im*dx,xcmin+(im+1)*dx,
		    				     ycmin+in*dy,ycmin+(in+1)*dy));
                }
            }
        }

        if (1 == 1) {
          cout << "\t output the coordinates of the cells to " << outName.gridXY2File << endl;
          outFile.open(outName.gridXY2File.data(), ios::out);

          for (i=0; i<contParam.nChips; i++) {         // make sure the ordering match with
              for (in=0; in<nc; in++) {      // that in filling Xmat
                  for (im=0; im<mc; im++) {
                      for (long col=0; col<(long)rowVec1d.size(); ++col) {
                          outFile << cell.at(i).at(im).at(in).at(col) << " ";
                      }
                      outFile << "\n";
                  }
              }
          }
          outFile.close();
        }

        cout << "\t place stars on the grid (regular grid is assumed) ... \n";
        cout << "\t\t one exposure at a time ... \n";

        /* For SNAP, all stars in an exposure is read in all at once.
         * For DES, stars are read in a chip at a time, so the read need to
         * loop over chips. After all chips have been read, the exposure is
         * processed further.
         */

        if (contParam.readDESorSNAP == 1 ) {               // DES pointing list
           ifstream infileExpList;
           string pointingVal;

           infileExpList.open(inName.DESexpListFile.data(), std::ios_base::in);
           if (infileExpList.is_open()) {
              while (! infileExpList.eof())
              {
                 while (getline(infileExpList, pointingVal, '\n')) {
                     exposureList.push_back (pointingVal);
                 }
              }
              infileExpList.close();
           }
           else {
              cout << "Error opening file " << inName.DESexpListFile << endl;
              exit(1);
           }

           contParam.ncols = exposureList.size();           // reset ncols
           cout << "\t number of exposures (ncols) is reset to " << contParam.ncols << endl;
        }

        resizeDataMat(contParam, myData);

        outFile.open(outName.starRecFile.data(), ios::out);

        if (!outFile) {
           cout << "Can't open output file " << outName.starRecFile << endl;
           exit(1);
        }

        for (iExp=0; iExp<contParam.ncols; iExp++) {      // loop over exposures

           vector< vector<double> > starPSF;    // (Nstar,PSF+2) x,y & PSF
           vector<int> starChip;                // (Nstar) which chip a star is on
           vector<double> starX, starY;         // (Nstar) star coordinates
           vector< vector<int> > starCell;      // (Nstar,2) which cell star is in
           vector<int> validateFlag;            // (Nstar) star in validation set (0) or not (1)

           if (contParam.readDESorSNAP == 2) {                 // SNAP
              // nameBase="data/SNAP/BS_randXY/PSFmomFP_";
              // nameBase="/astro/tutti1/mzm/SNAPdata/BS_randXY/PSFmomFP_";

              s=convertInt(iExp+1);
              fileName=inName.nameBaseSNAPpsf+s;
              starPSFfile = fileName.data();        // or use fileName.c_str()

              starPSF = readIn2dData(starPSFfile);
              Nstar=starPSF.size();

              starX.resize(Nstar);
              starY.resize(Nstar);
              starChip.resize(Nstar);

              for (istar=0; istar<Nstar; istar++) {   // loop over all stars in an exposure
                  starX.at(istar)=starPSF.at(istar).at(0);     // to assign star coordinate
                  starY.at(istar)=starPSF.at(istar).at(1);
              }

              for (istar=0; istar<Nstar; istar++) {       // assign star chips
                 starChip.at(istar)=-999;       // if a star doesn't land on
                                                // any chips, starChip=-999
                 for (iChip=0; iChip<contParam.nChips; iChip++) {
                   if (starPSF[istar][0] > chipBound[iChip][0] &&
                       starPSF[istar][0] < chipBound[iChip][2] &&
                       starPSF[istar][1] > chipBound[iChip][1] &&
                       starPSF[istar][1] < chipBound[iChip][3] ) {
                       starChip.at(istar)=iChip;
                   }
                 }
              }
           }
           else {                                  // DES
              read_des_exp(starPSF,starChip,starX,starY,exposureList,iExp,contParam.nChips,
                           inName,contParam);
              Nstar=starPSF.size();
              // cout << "\t Exp = " << iExp << "\t" << starPSF.size() << "\t" << starChip.size()
              //      << "\t" << starX.size() << "\t" << starY.size() << endl;
           }

           starCell.resize(Nstar, vector<int>(2, 0));
           starCount.resize( contParam.nChips, grid2D );      // reset to zero for each exposure
           validateFlag.resize(Nstar, 0);

           for (istar=0; istar<Nstar; istar++) {   // loop over all stars in an exposure
                                                   // to assign stars to cell
             if (starChip.at(istar) != -999 ) {
               for (im=0; im<mc; im++) {           // find the cell id in x direction
                   if (cell.at(starChip.at(istar)).at(im).at(0).at(2) >= starX.at(istar)) {
                      starCell.at(istar).at(0)=im;
                      break;
                   }
                   // add guard that star out of bound; --> starChip(istar)=-999 as signal
                   // only for SNAP. assume stars all land on chips for DES
               }

               for (in=0; in<nc; in++) {           // find the cell id in y direction
                   if (cell.at(starChip.at(istar)).at(0).at(in).at(3) >= starY.at(istar)) {
                      starCell.at(istar).at(1)=in;
                      break;
                   }
               }

               //if (ran01() > 0.01*validateSetPct) {
	       if (1) {
                  validateFlag.at(istar)=1;                     // 1 means not in the validation set
                  starCount.at(starChip.at(istar)).at(im).at(in) += 1;

                  // add PSF values to the cell ... need to generalize

                  irow=starChip.at(istar)*mc*nc + in*mc + im;    // position in data vector
                                // this ordering should match the cell coordinate output
                  for (iShape=0; iShape<contParam.nShapelet; iShape++) {
                     myData.Xmat(irow+iShape*nCellTot,iExp) += starPSF.at(istar).at(iStartShapelet+iShape);
                     // myData.Xmat(irow+nCellTot,iExp) += starPSF.at(istar).at(6); // e2 for SNAP for eg
                  }
		  cout << "star " << istar << " is in cell " <<starChip.at(istar)<<" "<< im <<" "<< in << " "
                       << starPSF.at(istar).at(3) <<" "<<starPSF.at(istar).at(4) <<endl;
                  // cout << "star " << istar << " is in cell (" << starCell.at(istar).at(0)
                  //      << ", " << starCell.at(istar).at(0) << ")\n";
               }
             }

             else { im=0; in=0; }

             outFile << iExp+1 << " " << starChip.at(istar)+1 << " " << im+1 << " " << in+1
                     << " " << starX.at(istar) << " " << starY.at(istar)
                     << " " << starPSF.at(istar).at(iStartShapelet)          // e1
                     << " " << starPSF.at(istar).at(iStartShapelet+1)        // e2
                     << " " << validateFlag.at(istar) << endl;

           }        // end of istar loop (assign star to cell)

           // last step is to find the average PSF in a cell

           myData.Nmasked[iExp]=0;
           for (iChip=0; iChip<contParam.nChips; iChip++) {
               for (im=0; im<mc; im++) {
                   for (in=0; in<nc; in++) {
                       irow=iChip*mc*nc + in*mc + im;     // position in data vector
                       icount=starCount.at(iChip).at(im).at(in);
                       if (icount == 0) {            // missing data
                          myData.dataMask(irow,iExp)=0;
                          for (iShape=0; iShape<contParam.nShapelet; iShape++) {

                             // indx=iExp*contParam.nrows+irow+iShape*nCellTot;
                             // myData.missIndex[indx]=0;
                             myData.Nmiss++;
                             myData.Nmasked[iExp]+=1;
                          }

                          contParam.icMissing=1;                          // set icMissing
                       }
                       else {                        // non-missing data
                          myData.dataMask(irow,iExp)=1;
                          for (iShape=0; iShape<contParam.nShapelet; iShape++) {
                             myData.Xmat(irow+iShape*nCellTot,iExp) /= icount;   // PSF = PSF/icount;
                          }
                       }
                   }    // end of in
               }        // end of im
           }            // end of iChip

           // starPSF.clear();         // doesn't need this b/c vector is out of scope at
           // starChip.clear();        // the end of the closing } of the loop
           // starCell.clear();
           // starX.clear();
           // starY.clear();

           starCount.clear();         // empty star count

        }          // end of loop over exposures


	outFile.close();

        checkEmptyCell(contParam,myData);

}


/* --------------------------------------------------------------------- */
void readInputParams(c_ControlParam &contParam, double & SNAPmaskPct, double & tol,
                     const char *fileName)
{
     std::ifstream ifs;
     vector<int> tempVec;
     int icount,Ibuf;

     vector<double> tempDoubleVec;
     double Dbuf;

     int numOfInt=20;

     ifs.open(fileName);

     if (ifs.is_open()) {
        icount=0;
        while ( icount < numOfInt ) {     // read in integers
            std::string line;
            getline(ifs, line);
            std::stringstream ss(line, std::ios_base::in);

            if (line[0] == '#' || line.empty()) { continue; }

            ss >> Ibuf;
            tempVec.push_back(Ibuf);
            icount++;
        }

        while ( ! ifs.eof() ) {         // read in misc double parameters
            std::string line;
            getline(ifs, line);
            std::stringstream ss(line, std::ios_base::in);

            if (line[0] == '#' || line.empty()) { continue; }

            ss >> Dbuf;
            tempDoubleVec.push_back(Dbuf);
        }

        ifs.close();
     }
     else {
        cout << "Unable to open file " << fileName << endl;
        exit(1);
     }

     contParam.iverbose  = tempVec.at(0);
     contParam.icSVD     = tempVec.at(1);
     contParam.icEM      = tempVec.at(2);
     contParam.icWiberg  = tempVec.at(3);
     contParam.icout     = tempVec.at(4);
     contParam.icMissing = tempVec.at(5);
     contParam.icMean    = tempVec.at(6);
     contParam.icDefocus = tempVec.at(7);
     contParam.kEigen    = tempVec.at(8);
     contParam.iterMax   = tempVec.at(9);
     contParam.iSeed     = tempVec.at(10);
     contParam.icGetMat  = tempVec.at(11);
     contParam.readDESorSNAP = tempVec.at(12);
     contParam.nChips    = tempVec.at(13);
     contParam.mc        = tempVec.at(14);
     contParam.nc        = tempVec.at(15);
     contParam.nShapelet = tempVec.at(16);
     contParam.nrows     = contParam.nChips*contParam.mc*contParam.nc
                           *contParam.nShapelet;    // modified with nChips
     contParam.ncols     = tempVec.at(17);
     contParam.nzTabCol  = tempVec.at(18);
     contParam.skip61    = tempVec.at(19);   // skip CCD 61
     
     SNAPmaskPct         = tempDoubleVec.at(0);   // non-random SNAP stars; % of masked
     tol                 = tempDoubleVec.at(1);   // EM iteration convergence tolerance


}

/* --------------------------------------------------------------------- */
void readFileNames(c_inFileName &inName, c_outFileName &outName,
                     const char *fileName)
{
     std::ifstream ifs;
     vector<std::string> tempVec;
     int icount;

     ifs.open(fileName);

     if (ifs.is_open()) {
        icount=0;
        while ( ! ifs.eof() ) {
            std::string line;
            getline(ifs, line);

            if (line[0] == '#' || line.empty()) { continue; }

            // for (int i=0; i<line.length(); i++)         // remove white spaces
            //     if (line[i] == ' ') line.erase(i,1);
            // OR use the following to remove white spaces
            line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

            tempVec.push_back(line);
        }
        ifs.close();
     }
     else {
        cout << "Unable to open file " << fileName << endl;
        exit(1);
     }

     if(tempVec.size()!=28) {
       cout << "Wrong number of parameters "<<" "<<tempVec.size()<<" in file " << fileName << endl;
       cout << "Should be 28" << endl;
        exit(1);
     }

     inName.inputDir         = tempVec.at(0);                  // input file names
     inName.controlParamFile = inName.inputDir + tempVec.at(1);
     inName.chipBoundFile    = inName.inputDir + tempVec.at(2);
     inName.nameBaseSNAPpsf  = tempVec.at(3);
     inName.dirBaseDES       = tempVec.at(4);
     inName.runID_DES        = tempVec.at(5);
     inName.DESexpListFile   = inName.inputDir + tempVec.at(6);
     inName.nameBaseZtab     = tempVec.at(7);

     outName.outputDir    = tempVec.at(8);                    // output file names
     outName.starRecFile  = outName.outputDir + tempVec.at(9);         // SVD & EM
     outName.gridXY2File  = outName.outputDir + tempVec.at(10);
     outName.dataMatFile  = outName.outputDir + tempVec.at(11);
     outName.dataMaskFile = outName.outputDir + tempVec.at(12);

     outName.eigenValSVDfile    = outName.outputDir + tempVec.at(13);  // SVD only
     outName.eigenVecSVDfile    = outName.outputDir + tempVec.at(14);
     outName.eigenCoefSVDfile   = outName.outputDir + tempVec.at(15);
     outName.reconSVDfile       = outName.outputDir + tempVec.at(16);
     outName.reconErrSVDfile    = outName.outputDir + tempVec.at(17);
     outName.reconErrPixSVDfile = outName.outputDir + tempVec.at(18);
     outName.reconErrExpSVDfile = outName.outputDir + tempVec.at(19);

     outName.eigenVecEMfile     = outName.outputDir + tempVec.at(20);  // EM only
     outName.eigenCoefEMfile    = outName.outputDir + tempVec.at(21);
     outName.reconEMfile        = outName.outputDir + tempVec.at(22);
     outName.reconErrEMfile     = outName.outputDir + tempVec.at(23);
     outName.reconErrPixEMfile  = outName.outputDir + tempVec.at(24);
     outName.reconErrExpEMfile  = outName.outputDir + tempVec.at(25);

     outName.defocusIDfile     = outName.outputDir + tempVec.at(26);   // defocus
     outName.defocusCoeffFile  = outName.outputDir + tempVec.at(27);

}
