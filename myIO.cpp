#include "myIO.h"
using namespace std;

ExtHDU* writeToFits(CCfits::FITS* file,const DMatrix &mat,string name)
{
  std::vector<long> ax(2,0);
  ax[0]=mat.nrows();
  ax[1]=mat.ncols();
  int n=mat.nrows()*mat.ncols();

  std::valarray<double> data(n);
  std::copy(mat.cptr(),mat.cptr()+n,&data[0]);
  ExtHDU* ext = file->addImage(name,DOUBLE_IMG,ax);
  long fpixel(1);
  ext->write(fpixel,n,data);
  return ext;
}
