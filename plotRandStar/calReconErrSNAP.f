c ----------- prepare data for PSF whisker plot ------------------------
c input x, y, e1, e2 of the stars from given exposure, and the
c       reconstructed PSF on a grid, calculate the difference.
c       Star locations are assumed to be random; We use the PSF at
c       the grid cell where the star is located as the fitted PSF.
c       Star cell information is read in from starRec.
c output x-dx/2, y-dy/2, dx, dy
c gnuplot plots the output by
c   plot 'file.dat' using 1:2:3:4 with vectors head filled lt 2
c-----------------------------------------------------------------------
	program main

        integer nstarMax,nstar,i,IDexp,ichar
        parameter (nstarMax=100000,Nexp=2000)
        real*8 x(nstarMax),y(nstarMax),x2,y2,
     c         e1data(nstarMax),e2data(nstarMax),
     c         de1(nstarMax),de2(nstarMax),RMS(Nexp),dum,RMSavg,RMSstd
        character*50 inputStarRec, inputE12recon, aChar,
     c               inputData, outputFile, dataDir,reconDir

        integer nChips,mc,nc,nShapeLet,im,in,iExp,iChip,npt,ipt,
     c          iExpBUF,iChipBUF,imBUF,inBUF,nLine,nDataDir,nReconDir
        parameter (nChips=44,mc=35,nc=35,nShapeLet=2,npt=nChips*mc*nc)

        real*8 e1recon(npt),e2recon(npt)
        real*8 dataMask(nstarMax),avg     ! dataMask is not used

        ! dataDir="../data/SNAP/BS_randXY/"
        ! reconDir="../results/"
        dataDir="/astro/tutti1/mzm/SNAPdata/BS_randXY/"
        reconDir="/data/mzm/SNAP_35x35_20Eigen/"
        call fileNameLen(dataDir,nDataDir)
        call fileNameLen(reconDir,nReconDir)

        inputStarRec=reconDir(1:nReconDir)//"starRec"
        inputE12recon=reconDir(1:nReconDir)//"reconEM"

        outputFile="SNAPreconErr.dat"

        RMSavg=0.d0
        RMSstd=0.d0

        open(20,file=inputStarRec)
        open(22,file=inputE12recon)
        open(23,file=outputFile)

        iExpBUF=-1
        nLine=0
        do IDexp = 1, Nexp

          read(22,*)(e1recon(i),i=1,npt), (e2recon(j),j=1,npt)

          call numToChar(IDexp,aChar,ichar)
          inputData=dataDir(1:nDataDir)//"PSFmomFP_"//aChar(1:ichar)
          open(21,file=inputData)

          nstar=0
          RMS(IDexp)=0.d0

          if (iExpBUF .eq. IDexp) then
            if (iChipBUF .ne. -998) then
              nstar=nstar+1
              read(21,*)x(nstar),y(nstar),dum,dum,dum,
     c              e1data(nstar),e2data(nstar)
              ipt=(iChipBUF-1)*mc*nc + (inBUF-1)*mc + imBUF
              de1(nstar)=e1recon(ipt)-e1data(nstar)
              de2(nstar)=e2recon(ipt)-e2data(nstar)
              RMS(IDexp)=RMS(IDexp)+de1(nstar)**2+de2(nstar)**2
            else
              read(21,*)
            endif
          endif

111       continue
            read(20,*,end=100)iExp,iChip,im,in,x2,y2
            nLine=nLine+1
            if (iExp .eq. IDexp) then
              if (iChip .ne. -998) then
                nstar=nstar+1
                read(21,*)x(nstar),y(nstar),dum,dum,dum,
     c                e1data(nstar),e2data(nstar)
                ipt=(iChip-1)*mc*nc + (in-1)*mc + im
                de1(nstar)=e1recon(ipt)-e1data(nstar)
                de2(nstar)=e2recon(ipt)-e2data(nstar)
                RMS(IDexp)=RMS(IDexp)+de1(nstar)**2+de2(nstar)**2
              else
                read(21,*)
              endif
              goto 111
            else            ! store the over-read line in buffer
              iExpBUF=iExp
              iChipBUF=iChip
              imBUF=im
              inBUF=in
            endif

100       continue
          close(21)

          RMS(IDexp)=DSQRT(RMS(IDexp)/2.d0/nstar)
          write(23,'(E15.6,2I8)')RMS(IDexp),nstar,nLine
          RMSavg=RMSavg+RMS(IDexp)
          RMSstd=RMSstd+RMS(IDexp)**2

        enddo

        close(20)
        close(22)
        close(23)

        RMSavg=RMSavg/Nexp
        RMSstd=DSQRT(RMSstd/(Nexp-1) - RMSavg**2*Nexp/(Nexp-1))
        write(*,'(A,2E15.6)')"  mean and std: ",RMSavg,RMSstd

        stop
        end
