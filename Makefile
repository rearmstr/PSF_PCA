# TMV is compiled using g++. c++ is an older version of g++
# GCC = g++
# CXXFLAGS = -pg -g
GCC = mpicxx
CXXFLAGS = -g -O3

# ============ MBPbnl =======================================================
# INCLUDE = -I/Users/mzm/opt/include -I/usr/local/include/CCfits
# LIBS = -L/Users/mzm/opt/lib -ltmv -ltmv_symband -lblas -lpthread -fopenmp -lCCfits -lcfitsio
#
# ============ BACH =========================================================
# INCLUDE = -I/global/data/products/Linux64/tmv/v0.64/include
# LIBS = -ltmv -ltmv_symband -lblas -lpthread -fopenmp
#
# ============ LSST (echo $TMV_DIR) =========================================
# INCLUDE = -I/lsst/u/esheldon/local/products/Linux64/tmv/v0.63/include
# INCLUDE = -I/astro/u/esheldon/exports/tmv-work/include -I/astro/u/esheldon/local/products/Linux64/ccfits/2.3/include
# LIBS = -L/astro/u/esheldon/exports/tmv-work/lib -ltmv -ltmv_symband -lblas -lpthread -lCCfits -lcfitsio

# for folio
#-I/usr/global/cfitsio/include
#INCLUDE = -I${TMV_DIR}/include -I/global/homes/r/rarmst/soft/include 
INCLUDE = -I${TMV_DIR}/include -I/usr/global/cfitsio/include -I/usr/global/CCfits/include
LIBS = -L${TMV_DIR}/lib -ltmv -ltmv_symband -lmkl_intel_lp64 -lmkl_core   \
       -lmkl_sequential -lpthread -Wl,-rpath=/usr/global/tmv0.71/lib \
       -openmp -L/usr/global/cfitsio/lib -L/usr/global/CCfits/lib -lCCfits -lcfitsio

fileName=runPCA

all: $(fileName)
clean:
	rm *.o
	rm $(fileName)

OBJtmv = $(fileName).o NR.o initialize.o myIO.o PCAcommon.o PCAuseSVD.o PCAuseEM.o \
	PCAuseWiberg.o rmDefocus.o ConfigFile.o PCAObjects.o

$(fileName).o: $(fileName).cpp myClass.h myIO.h Log.h
	$(GCC) $(CXXFLAGS) -c -o $@ $(fileName).cpp $(INCLUDE) 
$(fileName): $(OBJtmv)
	$(GCC) $(CXXFLAGS) -o $(fileName) $(OBJtmv) $(LIBS) -o $(fileName)

NR.o: NR.cpp NR.h
	$(GCC) $(CXXFLAGS) -c -o $@ NR.cpp

initialize.o: initialize.cpp initialize.h
	$(GCC) $(CXXFLAGS) -c -o $@ initialize.cpp $(INCLUDE)

myIO.o: myIO.cpp myIO.h
	$(GCC) $(CXXFLAGS) -c -o $@ myIO.cpp $(INCLUDE)

PCAcommon.o: PCAcommon.cpp PCAcommon.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAcommon.cpp $(INCLUDE)

PCAuseSVD.o: PCAuseSVD.cpp PCAuseSVD.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAuseSVD.cpp $(INCLUDE)

PCAuseEM.o: PCAuseEM.cpp PCAuseEM.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAuseEM.cpp $(INCLUDE)

PCAuseWiberg.o: PCAuseWiberg.cpp PCAuseWiberg.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAuseWiberg.cpp $(INCLUDE)

rmDefocus.o: rmDefocus.cpp rmDefocus.h
	$(GCC) $(CXXFLAGS) -c -o $@ rmDefocus.cpp $(INCLUDE)
PCAObjects.o: PCAObjects.cpp PCAObjects.h Log.h 
	$(GCC) $(CXXFLAGS) -c -o $@ PCAObjects.cpp $(INCLUDE)
