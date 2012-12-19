FORTRAN = ifort
FFLAGS  =  -O3 -vec-report0

all:
	make prepPSF
	make prepPSFrecon
	make prepPSFdiff
	make prepDESpsfDiff
	make prepPSFdiffMask
	make prepEigen
	make prepEigenDiff
	make prepMisAlignDiff
	make binData
	make calDES2pt
	make calDES_Rowe
	make calSNAP_Rowe
	make calReconErrSNAP

clean:
	rm -f *.o
	rm prepPSF prepPSFrecon prepPSFdiff prepEigen prepEigenDiff \
	   prepPSFdiffMask binData calDES2pt calDES_Rowe \
           calReconErrSNAP calSNAP_Rowe

common.o: commonRoutines.f
	$(FORTRAN) $(FFLAGS) -c -o $@ commonRoutines.f

prepPSF: prepPSF.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepPSF.o common.o
prepPSF.o: prepPSF.f
	$(FORTRAN) $(FFLAGS) -c prepPSF.f

prepPSFrecon: prepPSFrecon.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepPSFrecon.o common.o
prepPSFrecon.o: prepPSFrecon.f
	$(FORTRAN) $(FFLAGS) -c prepPSFrecon.f

prepPSFdiff: prepPSFdiff.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepPSFdiff.o common.o
prepPSFdiff.o: prepPSFdiff.f
	$(FORTRAN) $(FFLAGS) -c prepPSFdiff.f

prepDESpsfDiff: prepDESpsfDiff.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepDESpsfDiff.o common.o
prepDESpsfDiff.o: prepDESpsfDiff.f
	$(FORTRAN) $(FFLAGS) -c prepDESpsfDiff.f

prepPSFdiffMask: prepPSFdiffMask.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepPSFdiffMask.o common.o
prepPSFdiffMask.o: prepPSFdiffMask.f
	$(FORTRAN) $(FFLAGS) -c prepPSFdiffMask.f

prepEigen: prepEigen.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepEigen.o common.o
prepEigen.o: prepEigen.f
	$(FORTRAN) $(FFLAGS) -c prepEigen.f

prepEigenDiff: prepEigenDiff.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepEigenDiff.o common.o
prepEigenDiff.o: prepEigenDiff.f
	$(FORTRAN) $(FFLAGS) -c prepEigenDiff.f

prepMisAlignDiff: prepMisAlignDiff.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ prepMisAlignDiff.o common.o
prepMisAlignDiff.o: prepMisAlignDiff.f
	$(FORTRAN) $(FFLAGS) -c prepMisAlignDiff.f

binData: binData.o
	$(FORTRAN) $(FFLAGS) -o $@ binData.o
binData.o: binData.f
	$(FORTRAN) $(FFLAGS) -c binData.f

calDES2pt: calDES2pt.o
	$(FORTRAN) $(FFLAGS) -o $@ calDES2pt.o
calDES2pt.o: calDES2pt.f
	$(FORTRAN) $(FFLAGS) -c calDES2pt.f

calDES_Rowe: calDES_Rowe.o
	$(FORTRAN) $(FFLAGS) -o $@ calDES_Rowe.o
calDES_Rowe.o: calDES_Rowe.f
	$(FORTRAN) $(FFLAGS) -c calDES_Rowe.f

calSNAP_Rowe: calSNAP_Rowe.o
	$(FORTRAN) $(FFLAGS) -o $@ calSNAP_Rowe.o
calSNAP_Rowe.o: calSNAP_Rowe.f
	$(FORTRAN) $(FFLAGS) -c calSNAP_Rowe.f

calReconErrSNAP: calReconErrSNAP.o common.o
	$(FORTRAN) $(FFLAGS) -o $@ calReconErrSNAP.o common.o
calReconErrSNAP.o: calReconErrSNAP.f
	$(FORTRAN) $(FFLAGS) -c calReconErrSNAP.f
