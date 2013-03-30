# TMV is compiled using g++. c++ is an older version of g++
# GCC = g++
# CXXFLAGS = -pg -g
GCC = mpicxx
CXXFLAGS = -O3 -g

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
INCLUDE = -I /data2/home/rarmst/soft/tmv0.71_icpc/include/ -I/usr/global/CCfits/include -I/usr/global/cfitsio/include
LIBS = -L/data2/home/rarmst/soft/tmv0.71_icpc/lib -ltmv -ltmv_symband -lmkl_intel_lp64 -lmkl_core   \
       -lmkl_sequential -lpthread -Wl,-rpath=/usr/global/tmv0.71/lib \
       -openmp -L/usr/global/CCfits/lib -lCCfits -L/usr/global/CCfits/lib -lcfitsio

fileName=PSF_PCA

all: $(fileName)
clean:
	rm *.o
	rm $(fileName)

OBJtmv = $(fileName).o NR.o initialize.o myIO.o PCAcommon.o PCAuseSVD.o PCAuseEM.o \
PCAuseWiberg.o rmDefocus.o ConfigFile.o

$(fileName).o: $(fileName).cpp myClass.h
	$(GCC) $(CXXFLAGS) -c -o $@ $(fileName).cpp $(INCLUDE)
$(fileName): $(OBJtmv)
	$(GCC) $(CXXFLAGS) -o $(fileName) $(OBJtmv) $(LIBS)

#test_svd.o: test_svd.cc myClass.h
#	$(GCC) $(CXXFLAGS) -c -o $@ test_svd.cc $(INCLUDE)
#test_svd: test_svd.o NR.o initialize.o myIO.o PCAcommon.o PCAuseSVD.o PCAuseEM.o \
#	PCAuseWiberg.o rmDefocus.o
#	$(GCC) $(CXXFLAGS) -o test_svd $^ $(LIBS)

NR.o: NR.cpp NR.h
	$(GCC) $(CXXFLAGS) -c -o $@ NR.cpp

initialize.o: initialize.cpp initialize.h
	$(GCC) $(CXXFLAGS) -c -o $@ initialize.cpp $(INCLUDE)

myIO.o: myIO.cpp myIO.h
	$(GCC) $(CXXFLAGS) -c -o $@ myIO.cpp $(INCLUDE)

PCAcommon.o: PCAcommon.cpp PCAcommon.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAcommon.cpp $(INCLUDE)

ConfigFile.o: ConfigFile.cpp ConfigFile.h
	$(GCC) $(CXXFLAGS) -c -o $@ ConfigFile.cpp $(INCLUDE)

PCAuseSVD.o: PCAuseSVD.cpp PCAuseSVD.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAuseSVD.cpp $(INCLUDE)

PCAuseEM.o: PCAuseEM.cpp PCAuseEM.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAuseEM.cpp $(INCLUDE)

PCAuseWiberg.o: PCAuseWiberg.cpp PCAuseWiberg.h
	$(GCC) $(CXXFLAGS) -c -o $@ PCAuseWiberg.cpp $(INCLUDE)

rmDefocus.o: rmDefocus.cpp rmDefocus.h
	$(GCC) $(CXXFLAGS) -c -o $@ rmDefocus.cpp $(INCLUDE)
