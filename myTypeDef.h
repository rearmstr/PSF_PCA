#ifndef myTypeDef_H
#define myTypeDef_H

#include "TMV.h"

typedef tmv::Matrix<double> DMatrix;
typedef tmv::MatrixView<double> DMatrixView;
typedef tmv::Matrix<float> FMatrix;
typedef tmv::MatrixView<float> FMatrixView;
typedef tmv::DiagMatrix<double> DDiagMatrix;
typedef tmv::DiagMatrix<float> FDiagMatrix;
typedef tmv::Matrix<double,tmv::RowMajor> DMatrixRowM;
typedef tmv::Vector<double> DVector;
typedef tmv::Vector<float> FVector;
typedef tmv::VectorView<double> DVectorView;
typedef tmv::VectorView<float> FVectorView;

#endif
