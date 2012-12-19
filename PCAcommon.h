#ifndef PCAcommon_H
#define PCAcommon_H

void subtractMean(c_ControlParam &,c_Data &);          // subtract mean vector from data
void initializeTMVmat(DMatrix &mat, double initialVal);
void resizeDataMat(c_ControlParam &, c_Data &);        // resize data matrices accord to input
void checkEmptyCell(c_ControlParam &, c_Data &);       // fill empty cells with [0,1] random number

#endif
