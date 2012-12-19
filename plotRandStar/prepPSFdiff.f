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
        parameter (nstarMax=100000,IDexp=16)
        real*8 x(nstarMax),y(nstarMax),x2,y2,
     c         e1data(nstarMax),e2data(nstarMax),
     c         de1(nstarMax),de2(nstarMax),scaleData,scaleRes,dum
        character*50 inputStarRec, inputE12recon, outputData, aChar,
     c               inputData, outputDiff, dataDir,reconDir

        integer nChips,mc,nc,nShapeLet,im,in,iExp,iChip,npt,ipt,
     c          nDataDir,nReconDir
        parameter (nChips=44,mc=16,nc=16,nShapeLet=2,npt=nChips*mc*nc)

        real*8 e1recon(npt),e2recon(npt)
        real*8 dataMask(nstarMax),avg     ! dataMask is not used


        ! dataDir="../data/SNAP/BS_randXY/"
        ! reconDir="../results/"
        dataDir="/astro/tutti1/mzm/SNAPdata/BS_randXY/"
        reconDir="/data/mzm/SNAP_16x16_20Eigen/"
        call fileNameLen(dataDir,nDataDir)
        call fileNameLen(reconDir,nReconDir)

        inputStarRec=reconDir(1:nReconDir)//"starRec"
        inputE12recon=reconDir(1:nReconDir)//"reconEM"

        call numToChar(IDexp,aChar,ichar)
        inputData=dataDir(1:nDataDir)//"PSFmomFP_"//aChar(1:ichar)

        outputData="whiskerData.dat"
        outputDiff="whiskerDiff.dat"

        ! scaleData=0.3
        ! scaleRes=3.0
        scaleData=0.5
        scaleRes=2.0

        open(21,file=inputE12recon)
        do i = 1, IDexp-1
          read(21,*)
        enddo
        read(21,*)(e1recon(i),i=1,npt), (e2recon(j),j=1,npt)
        close(21)

        open(20,file=inputStarRec)
        open(21,file=inputData)
        nstar=0
111     continue
          read(20,*,end=100)iExp,iChip,im,in,x2,y2
          ! iExp=iExp+1          ! change start point from 0 --> 1
          ! iChip=iChip+1
          ! im=im+1
          ! in=in+1
          if (iExp .eq. IDexp) then
            if (iChip .ne. -998) then
              nstar=nstar+1
              read(21,*)x(nstar),y(nstar),dum,dum,dum,
     c                e1data(nstar),e2data(nstar)
              ipt=(iChip-1)*mc*nc + (in-1)*mc + im
              de1(nstar)=e1recon(ipt)-e1data(nstar)
              de2(nstar)=e2recon(ipt)-e2data(nstar)
              ! write(*,'(4E15.6)')x(nstar),x2,y(nstar),y2
            else
              read(21,*)
            endif
          endif
        goto 111

100     continue
        close(20)
        close(21)
        write(*,'(A,I5)')"      nstar = ",nstar

        call genWhisker(x,y,e1data,e2data,nstar,dataMask,0,scaleData,
     c                        outputData,avg)
        write(*,'(A,f12.6)')"      average e ",avg

        call genWhisker(x,y,de1,de2,nstar,dataMask,0,scaleRes,
     c                        outputDiff,avg)
        write(*,'(A,f12.6)')"      average de ",avg

        stop
        end
