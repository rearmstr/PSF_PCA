
#include <sstream>
#include <valarray>
#include <CCfits/CCfits>
#include "dbg.h"
#include <fitsio.h>
#include <sys/stat.h>
#include "Image.h"
#include <stdexcept>


struct FileNotFoundException : public std::runtime_error
{
    FileNotFoundException(const std::string& filename) throw() :
        std::runtime_error("Error: file "+filename+" not found")
    {}
};

struct ReadException : public std::runtime_error
{
    ReadException(const std::string& msg) throw() :
        std::runtime_error(msg)
    {}
};

struct WriteException : public std::runtime_error
{
    WriteException(const std::string& msg) throw() :
        std::runtime_error(msg)
    {}
};

bool static DoesFileExist(const std::string& file_name)
{
    struct stat file_info;
    int stat_value;

    // Attempt to get the file attributes
    stat_value = stat(file_name.c_str(),&file_info);
    if (stat_value == 0) {
        // We were able to get the file attributes
        // so the file obviously exists.
        return true;
    } else {
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        return false;
    }
}



template <typename T> 
Image<T>::Image(std::string filename, int hdu) :
    _filename(filename), _hdu(hdu)
{
    readFits();
    _loaded=true;
}

template <typename T> 
void Image<T>::load(std::string filename, int hdu) 
{
    _filename = filename;
    _hdu = hdu;
    readFits();
    _loaded=true;
}

template <typename T> 
inline int getBitPix() { return 0; }

template <> 
inline int getBitPix<double>() { return DOUBLE_IMG; }

template <> 
inline int getBitPix<float>() { return FLOAT_IMG; }

template <typename T> 
inline int getDataType() { return 0; }

template <> 
inline int getDataType<double>() { return TDOUBLE; }

template <> 
inline int getDataType<float>() { return TFLOAT; }

template <typename T> 
void Image<T>::readFits()
{
    std::stringstream err_msg;

    // Need a check here
    //xdbg<<"Start Image::readFits"<<std::endl;
    //xdbg<<"filename = "<<_filename<<std::endl;
    //xdbg<<"hdu = "<<_hdu<<std::endl;

    if (!DoesFileExist(_filename)) {
        throw FileNotFoundException(_filename);
    }

#if 0
    // New CCFits implementation
    // This isn't working for me.  
    // I think it is because all three extensions have the same name:
    // COMPRESSED_IMAGE.
    // I've sent email to the CCfits help desk to see if they can help
    // me figure out the right way to do this.
    // But for now, stick with cfitsio commands.
    try {
        CCfits::FITS fits(_filename, CCfits::Read);
        //xdbg<<"Made fits object\n";
       //  if (XDEBUG) {
// 	  typedef std::multimap<std::string, CCfits::ExtHDU*> extmap_type;
// 	  const extmap_type& extmap = const_cast<const CCfits::FITS&>(fits).extension();
// 	  for(extmap_type::const_iterator i=extmap.begin();i!=extmap.end();++i) {
// 	    //xdbg<<"Extension "<<i->first<<" = "<<i->second<<std::endl;
// 	  }
//         }

        std::valarray<double> data;
        if (_hdu == 1) {
            CCfits::PHDU& image=fits.pHDU();
            //xdbg<<"Got primary hdu object"<<std::endl;
            image.readAllKeys();
            //xdbg<<"read all keys:"<<std::endl;
            //xdbg<<image<<std::endl;
            //xdbg<<"axes = "<<image.axes()<<std::endl;
            if (image.axes() != 2)
                throw std::runtime_error("Number of axes != 2");
            image.read(data);
            //xdbg<<"read data into valarray\n";
            _xmax = image.axis(0);
            _ymax = image.axis(1);
        } else {
            fits.read(_hdu-1);
            //xdbg<<"read extension"<<std::endl;
            CCfits::ExtHDU& image=fits.extension(_hdu-1);
            //xdbg<<"Got extension hdu object"<<std::endl;
            image.readAllKeys();
            //xdbg<<"read all keys:"<<std::endl;
            //xdbg<<image<<std::endl;
            //xdbg<<"axes = "<<image.axes()<<std::endl;
            if (image.axes() != 2)
                throw std::runtime_error("Number of axes != 2");
            image.read(data);
            //xdbg<<"read data into valarray\n";
            _xmax = image.axis(0);
            _ymax = image.axis(1);
        }
        _xmin = 0;
        _ymin = 0;
        //xdbg<<"size = "<<_xmax<<" , "<<_ymax<<std::endl;
        _source.reset(new tmv::Matrix<T>(_xmax,_ymax));
        //xdbg<<"done make matrix of image"<<std::endl;
        //xdbg<<"data.size = "<<data.size()<<" =? "<<_xmax*_ymax<<std::endl;
        Assert(int(data.size()) == _xmax*_ymax);
        std::copy(&data[0],&data[data.size()-1],_source->ptr());
        //xdbg<<"done copy data to _source\n";
        _m.reset(new tmv::MatrixView<T>(_source->view()));
        //xdbg<<"Done make matrixview"<<std::endl;
    } catch (CCfits::FitsException& e) {
      //xdbg<<"Caught FitsException: "<<e.message()<<std::endl;
      throw ReadException(
			  "Error reading from " + _filename + 
			  //" hdu " + ConvertibleString(_hdu) +
			  " -- caught error\n" + e.message());
    } catch (std::exception& e) {
      //xdbg<<"Caught std::exception: "<<e.what()<<std::endl;
      throw ReadException(
			  "Error reading from " + _filename + 
			  //" hdu " + ConvertibleString(_hdu) +
			  " -- caught error\n" + e.what());
    } catch (...) {
        //xdbg<<"Caught other exception: "<<std::endl;
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " -- caught unknown error\n");
    }
#else
    // Old cfitsio implementation.
    // Remove once above is working correctly.
    fitsfile *fPtr;
    int fitsErr=0;

    fits_open_file(&fPtr,_filename.c_str(),READONLY,&fitsErr);
    //xdbg<<"Done open"<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error opening fits file " + _filename);
    }

    fits_movabs_hdu(fPtr,_hdu,0,&fitsErr);
    //xdbg<<"Moved to hdu "<<_hdu<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
			    "Error reading from " + _filename);// + 
	//" moving to hdu " + ConvertibleString(_hdu));
    }

    int bitPix, nAxes;
    long sizes[2];
    fits_get_img_param(fPtr, int(2), &bitPix, &nAxes, sizes, &fitsErr);
    //xdbg<<"done getimgparam"<<std::endl;
    //xdbg<<"naxes = "<<nAxes<<std::endl;
    //xdbg<<"bitpix = "<<bitPix<<std::endl;
    //xdbg<<"FLOAT_IMG = "<<FLOAT_IMG<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " reading image parameters");
    }
    if (nAxes != 2) {
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " Number of axes != 2");
    }
    //xdbg<<"sizes = "<<sizes[0]<<"  "<<sizes[1]<<std::endl;

    _xmin = 0;
    _xmax = sizes[0];
    _ymin = 0;
    _ymax = sizes[1];
    _source.reset(new tmv::Matrix<T>(_xmax,_ymax));
    //xdbg<<"done make matrix of image"<<std::endl;

    long fPixel[2] = {1,1};
    int anynul;
    //xdbg<<"Before read_pix\n";
    Assert(getDataType<T>());
    fits_read_pix(fPtr,getDataType<T>(),fPixel,long(_xmax*_ymax),0,
                  _source->ptr(),&anynul,&fitsErr);
    //xdbg<<"done readpix\n";
    //xdbg<<"anynul = "<<anynul<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " reading pixel data");
    }

    _m.reset(new tmv::MatrixView<T>(_source->view()));
    //xdbg<<"Done make matrixview"<<std::endl;

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error closing fits file " + _filename);
    }
#endif

    //xdbg<<"Done Image ReadFits"<<std::endl;
}

// Load Subimage
// These next three are basically the same as the above constructors,
// and readFits function, but with a range of pixels to read in.
// I could probably combine these pretty easily, since most of the
// code is identical, but for now there is some significant code 
// repetition here.
// template <typename T>
// Image<T>::Image(
//     const ConfigFile& params, std::auto_ptr<Image<T> >& weight_image,
//     int x1, int x2, int y1, int y2)
// {
//     _filename = MakeName(params,"image",true,true);
//     _hdu = GetHdu(params,"image",_filename,1);
//     readFits(x1,x2,y1,y2);
//     //xdbg<<"Opened image "<<_filename<<std::endl;

//     // Load weight image (if necessary)
//     if (params["noise_method"] == "WEIGHT_IMAGE") {
//         std::string weightName = MakeName(params,"weight",true,true);
//         int weight_hdu = GetHdu(params,"weight",weightName,1);
//         weight_image.reset(new Image<T>(weightName,weight_hdu,x1,x2,y1,y2));
//        //dbg<<"Opened weight image.\n";

//         // Make sure any bad pixels are marked with 0 variance.
//         if (params.keyExists("badpix_file") || params.keyExists("badpix_ext")) {
//             std::string badpixName = MakeName(params,"badpix",true,true);
//             int badpix_hdu = GetHdu(params,"badpix",badpixName,1);
//            //dbg<<"badpix name = "<<badpixName<<std::endl;
//            //dbg<<"hdu = "<<badpix_hdu<<std::endl;
//             Image<double> badpixIm(badpixName,badpix_hdu,x1,x2,y1,y2);
//            //dbg<<"Opened badpix image.\n";

//             for(int i=0;i<=weight_image->getMaxI();++i) {
//                 for(int j=0;j<=weight_image->getMaxJ();++j) {
//                     if (badpixIm(i,j) > 0.0) (*weight_image)(i,j) = 0.0;
//                 }
//             }
//         }
//     }
// }

// template <typename T>
// Image<T>::Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weight_image,
//                 const Bounds& bounds)
// {
//     int x1 = int(floor(bounds.getXMin()));
//     int x2 = int(ceil(bounds.getXMax()));
//     int y1 = int(floor(bounds.getYMin()));
//     int y2 = int(ceil(bounds.getYMax()));
//     _filename = MakeName(params,"image",true,true);
//     _hdu = GetHdu(params,"image",_filename,1);
//     readFits(x1,x2,y1,y2);
//     //xdbg<<"Opened image "<<_filename<<std::endl;

//     // Load weight image (if necessary)
//     if (params["noise_method"] == "WEIGHT_IMAGE") {
//         std::string weightName = MakeName(params,"weight",true,true);
//         int weight_hdu = GetHdu(params,"weight",weightName,1);
//         weight_image.reset(new Image<T>(weightName,weight_hdu,bounds));
//        //dbg<<"Opened weight image.\n";

//         // Make sure any bad pixels are marked with 0 variance.
//         if (params.keyExists("badpix_file") || params.keyExists("badpix_ext")) {
//             std::string badpixName = MakeName(params,"badpix",true,true);
//             int badpix_hdu = GetHdu(params,"badpix",badpixName,1);
//            //dbg<<"badpix name = "<<badpixName<<std::endl;
//            //dbg<<"hdu = "<<badpix_hdu<<std::endl;
//             Image<double> badpixIm(badpixName,badpix_hdu,bounds);
//            //dbg<<"Opened badpix image.\n";

//             for(int i=0;i<=weight_image->getMaxI();++i) {
//                 for(int j=0;j<=weight_image->getMaxJ();++j) {
//                     if (badpixIm(i,j) > 0.0) (*weight_image)(i,j) = 0.0;
//                 }
//             }
//         }
//     }
// }

// template <typename T>
// Image<T>::Image(const ConfigFile& params, int x1, int x2, int y1, int y2)
// {
//     _filename = MakeName(params,"image",true,true);
//     _hdu = GetHdu(params,"image",_filename,1);
//     readFits(x1,x2,y1,y2);
//     //xdbg<<"Opened image "<<_filename<<std::endl;
// }

// template <typename T>
// Image<T>::Image(const ConfigFile& params, const Bounds& bounds)
// {
//     int x1 = int(floor(bounds.getXMin()));
//     int x2 = int(ceil(bounds.getXMax()));
//     int y1 = int(floor(bounds.getYMin()));
//     int y2 = int(ceil(bounds.getYMax()));
//     _filename = MakeName(params,"image",true,true);
//     _hdu = GetHdu(params,"image",_filename,1);
//     readFits(x1,x2,y1,y2);
//     //xdbg<<"Opened image "<<_filename<<std::endl;
// }

template <typename T> 
Image<T>::Image(
    std::string filename, int hdu,
    int x1, int x2, int y1, int y2) :
    _filename(filename), _hdu(hdu)
{
    readFits(x1,x2,y1,y2);
}

template <typename T> 
Image<T>::Image(std::string filename, int hdu, const Bounds<T>& bounds) :
    _filename(filename), _hdu(hdu)
{
    int x1 = int(floor(bounds.getXMin()));
    int x2 = int(ceil(bounds.getXMax()));
    int y1 = int(floor(bounds.getYMin()));
    int y2 = int(ceil(bounds.getYMax()));
    readFits(x1,x2,y1,y2);
}

template <typename T> 
void Image<T>::readFits(int x1, int x2, int y1, int y2) 
{
    //xdbg<<"Start read fitsimage"<<std::endl;
    //xdbg<<"filename = "<<_filename<<std::endl;
    //xdbg<<"hdu = "<<_hdu<<std::endl;
    // TODO: Use CCFits
    fitsfile *fPtr;
    int fitsErr=0;

    fits_open_file(&fPtr,_filename.c_str(),READONLY,&fitsErr);
    //xdbg<<"Done open"<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error opening fits file " + _filename);
    }


    fits_movabs_hdu(fPtr,_hdu,0,&fitsErr);
    //xdbg<<"Moved to hdu "<<_hdu<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
			    "Error reading from " + _filename);// + 
            //" moving to hdu " + ConvertibleString(_hdu));
    }

    int bitPix, nAxes;
    long sizes[2];
    fits_get_img_param(fPtr, int(2), &bitPix, &nAxes, sizes, &fitsErr);
    //xdbg<<"done getimgparam"<<std::endl;
    //xdbg<<"naxes = "<<nAxes<<std::endl;
    //xdbg<<"bitpix = "<<bitPix<<std::endl;
    //xdbg<<"FLOAT_IMG = "<<FLOAT_IMG<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " getting image parameters");
    }
    if (nAxes != 2) {
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " number of axes != 2");
    }
    //xdbg<<"sizes = "<<sizes[0]<<"  "<<sizes[1]<<std::endl;

    _xmin = 0;
    _xmax = sizes[0];
    _ymin = 0;
    _ymax = sizes[1];

    if (_xmin < x1) _xmin = x1;
    if (_xmax > x2) _xmax = x2; if (_xmax < _xmin+1) _xmax = _xmin+1;
    if (_ymin < y1) _ymin = y1;
    if (_ymax > y2) _ymax = y2; if (_ymax < _ymin+1) _ymax = _ymin+1;
    _source.reset(new tmv::Matrix<T>(_xmax-_xmin,_ymax-_ymin));
    //xdbg<<"done make matrix of image"<<std::endl;

    long fPixel[2] = {_xmin+1,_ymin+1};
    long lPixel[2] = {_xmax,_ymax};
    long inc[2] = {1,1};
    int anynul;
    //xdbg<<"Before read_subset\n";
    Assert(getDataType<T>());
    fits_read_subset(fPtr,getDataType<T>(),fPixel,lPixel,inc,
                     0,_source->ptr(),&anynul,&fitsErr);
    //xdbg<<"done read_subset\n";
    //xdbg<<"anynul = "<<anynul<<std::endl;
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error reading from " + _filename + 
            //" hdu " + ConvertibleString(_hdu) +
            " reading pixel data");
    }

    _m.reset(new tmv::MatrixView<T>(_source->view()));
    //xdbg<<"Done make matrixview"<<std::endl;

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw ReadException(
            "Error closing fits file " + _filename);
    }
    //xdb<<"Leaving Image ReadFits"<<std::endl;
}

template <typename T> 
void Image<T>::flush(std::string filename, int hdu) const
{
    _filename = filename;
    _hdu = hdu;
    flush();
}

template <typename T> 
void Image<T>::flush() const
{
    Assert(_filename != "");

    fitsfile *fPtr;
    int fitsErr=0;
    fits_open_file(&fPtr,_filename.c_str(),READWRITE,&fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
            "Error opening fits file " + _filename);
    }

    if (_hdu != 1) {
        int hduType;
        fits_movabs_hdu(fPtr,_hdu,&hduType,&fitsErr);
        if (fitsErr != 0) {
            fits_report_error(stderr,fitsErr);
            throw WriteException(
				 "Error reading from " + _filename);// + 
                //" moving to hdu " + ConvertibleString(_hdu));
        }
    }

    long fPixel[2] = {1,1};
    Assert(getDataType<T>());
    fits_write_pix(fPtr,getDataType<T>(),fPixel,long(_xmax*_ymax),
                   _m->ptr(),&fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
			     "Error reading from " + _filename);//  + 
//             " hdu " + ConvertibleString(_hdu) +
//             " writing pixel data");
    }

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
            "Error closing fits file " + _filename);
    }
}

template <typename T> 
void Image<T>::write(std::string filename) const
{
    _filename = filename;
    _hdu = 1;

    fitsfile *fPtr;
    int fitsErr=0;
    fits_create_file(&fPtr,("!"+_filename).c_str(),&fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
            "Error creating fits file " + _filename);
    }

    Assert(getBitPix<T>());
    int bitPix = getBitPix<T>();
    int nAxes = 2;
    long sizes[2] = { _m->colsize(), _m->rowsize() };
    fits_create_img(fPtr, bitPix, nAxes, sizes, &fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
            "Error writing " + _filename + 
            " creating image");
    }

    long fPixel[2] = {1,1};
    Assert(getDataType<T>());
    fits_write_pix(fPtr,getDataType<T>(),fPixel,long(_xmax*_ymax),
		   _m->ptr(),&fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
            "Error reading from " + _filename + 
            " writing pixel data");
    }

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
            "Error closing fits file " + _filename);
    }
}

template <typename T> 
std::vector<Image<T>*> Image<T>::divide(int nX, int nY) const
{
    std::vector<int> x(nX+1);
    std::vector<int> y(nY+1);
    x[0] = _xmin;  x[nX] = _xmax;
    y[0] = _ymin;  y[nY] = _ymax;
    int xstep = (_xmax-_xmin)/nX;
    int ystep = (_ymax-_ymin)/nY;
    for(int i=1;i<nX;++i) x[i] = x[i-1]+xstep;
    for(int j=1;j<nY;++j) y[j] = y[j-1]+ystep;
    std::vector<Image*> blockImages;
    blockImages.reserve(nX*nY);
    for(int i=0;i<nX;++i) for(int j=0;j<nY;++j) {
      blockImages.push_back(
			    new Image(_m->subMatrix(x[i],x[i+1],y[j],y[j+1]),
				      x[i],x[i+1],y[j],y[j+1]));
    }
    return blockImages;
}

template <typename T> 
T Image<T>::interpolate(double x, double y) const
{
    Assert(x>=double(_xmin) && x<double(_xmax));
    Assert(y>=double(_ymin) && y<double(_ymax));
    int i = int(floor(x-0.5));
    double dx = x - (i+0.5);
    int j = int(floor(y-0.5));
    double dy = y - (j+0.5);
    Assert(dx >= 0. && dx < 1.);
    Assert(dy >= 0. && dy < 1.);

    /* 
       2    3
       x      the point (x,y) is within square of points 0,1,2,3 as shown

       0    1 

       Since the points are really the values at the center of the pixel, it
       is possible for x to fall closer to the edge of the chip than any points.
       In this case, the function does an extrapolation given the nearest
       square of pixels.
       */
    if (i==int(_xmax-1)) {--i; dx += 1.;}
    if (j==int(_ymax-1)) {--j; dy += 1.;}
    if (i==-1) {++i; dx -= 1.;}
    if (j==-1) {++j; dy -= 1.;}
    Assert(i>=0 && j>=0 && i+1<int(_m->colsize()) && 
           j+1<=int(_m->rowsize()));

    T f0 = (*_m)(i,j);
    T f1 = (*_m)(i+1,j);
    T f2 = (*_m)(i+1,j);
    T f3 = (*_m)(i+1,j+1);
    T dfdx = f1-f0;
    T dfdy = f2-f0;
    T d2fdxdy = f3+f0-f1-f2;
    return f0 + dfdx*dx + dfdy*dy + d2fdxdy*dx*dy;
}

template <typename T> 
T Image<T>::quadInterpolate(double x, double y) const
{
    Assert(x>=_xmin && x< _xmax);
    Assert(y>=_ymin && y< _ymax);
    int i = int (floor(x));
    double dx = x - (i+0.5);
    int j = int (floor(y));
    double dy = y - (j+0.5);
    Assert(i<int(_m->colsize()));
    Assert(j<int(_m->rowsize()));
    Assert (std::abs(dx) <= 0.5);
    Assert (std::abs(dy) <= 0.5);

    /*
       7   4   8

       x          (x,y) is closer to point 0 than any other.
       1   0   2 


       5   3   6

       If any points are off the edge, we set them to the value they would
       have if the second derivative were 0 there.
       */
    T f0 = (*_m)(i,j);
    T f1 = (i > 0) ? (*_m)(i-1,j) : 0.;
    T f2 = (i < int(_m->colsize())-1) ? (*_m)(i+1,j) : 0.;
    T f3 = (j > 0) ? (*_m)(i,j-1) : 0.;
    T f4 = (j < int(_m->rowsize())-1) ? (*_m)(i,j+1) : 0.;
    T f5 = (i > 0 && j > 0) ? (*_m)(i-1,j-1) : 0.;
    T f6 = (i < int(_m->colsize())-1 && j > 0) ? (*_m)(i+1,j-1) : 0.;
    T f7 = (i > 0 && j < int(_m->rowsize())-1) ? (*_m)(i-1,j+1) : 0.;
    T f8 = (i < int(_m->colsize())-1 && j < int(_m->rowsize())-1) ?
        (*_m)(i+1,j+1) : 0.;
    if (i == 0) {
        f1 = 2*f0 - f2;
        f5 = 2*f3 - f6;
        f7 = 2*f4 - f8;
    }
    if (i == int(_m->colsize())-1) {
        f2 = 2*f0 - f1;
        f6 = 2*f3 - f5;
        f8 = 2*f4 - f7;
    }
    if (j == 0) {
        f3 = 2*f0 - f4;
        f5 = 2*f1 - f7;
        f6 = 2*f2 - f8;
    }
    if (j == int(_m->rowsize())-1) {
        f4 = 2*f0 - f3;
        f7 = 2*f1 - f5;
        f8 = 2*f2 - f6;
    }
    T dfdx = (f2-f1)/2.;
    T dfdy = (f4-f3)/2.;
    T d2fdx2 = (f1+f2-2.*f0);
    T d2fdy2 = (f3+f4-2.*f0);
    T d2fdxdy = (f5+f8-f7-f6)/4.;
    T temp = f0 + dfdx*dx + dfdy*dy + 0.5*d2fdx2*dx*dx + 0.5*d2fdy2*dy*dy +
        d2fdxdy*dx*dy;
    return temp;
}

template <typename T> 
T Image<T>::median() const
{
    std::vector<T> pixels;
    const int n1 = _m->colsize();
    const int n2 = _m->rowsize();
    pixels.reserve(n1*n2);
    for(int i=0;i<n1;++i) for(int j=0;j<n2;++j) {
      pixels.push_back((*_m)(i,j));
    }
    sort(pixels.begin(),pixels.end());
    return pixels[pixels.size()/2];
}

template class Image<double>;
//template class Image<float>;
