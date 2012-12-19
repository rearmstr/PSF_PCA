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
        parameter (nstarMax=100000,IDexp=3)
        real*8 x(nstarMax),y(nstarMax),x2,y2,
     c         e1data(nstarMax),e2data(nstarMax),
     c         de1(nstarMax),de2(nstarMax),scale,dum
        character*50 inputStarRec, inputE12recon, outputData, aChar,
     c               inputData, outputDiff

        integer nChips,mc,nc,nShapeLet,im,in,iExp,iChip,npt,ipt
        parameter (nChips=62,mc=20,nc=40,nShapeLet=2)
        parameter (npt=nChips*mc*nc)

        real*8 e1recon(npt),e2recon(npt)
        real*8 dataMask(nstarMax),avg     ! dataMask is not used

        call numToChar(IDexp,aChar,ichar)

        inputStarRec="starXYe12.dat"
        inputE12recon="/data/mzm/DES_20x40_681_tol-7_11t/reconEM"

        outputData="whiskerData.dat"
        outputDiff="whiskerDiff.dat"

        scale=5.0

        open(21,file=inputE12recon)
        do i = 1, IDexp-1
          read(21,*)
        enddo
        read(21,*)(e1recon(i),i=1,npt), (e2recon(j),j=1,npt)
        close(21)

        open(20,file=inputStarRec)
        nstar=0
111     continue
          nstar=nstar+1
          read(20,*,end=100)iExp,iChip,im,in,x2,y2,
     c            e1data(nstar),e2data(nstar),x(nstar),y(nstar)
          e1data(nstar)=e1data(nstar)/sqrt(2.0)
          e2data(nstar)=-e2data(nstar)/sqrt(2.0)
          ! iExp=iExp+1          ! change start point from 0 --> 1
          ! iChip=iChip+1
          ! im=im+1
          ! in=in+1
          ipt=(iChip-1)*mc*nc + (in-1)*mc + im
          de1(nstar)=e1recon(ipt)/sqrt(2.0)-e1data(nstar)
          de2(nstar)=-e2recon(ipt)/sqrt(2.0)-e2data(nstar)
          ! write(*,'(4E15.6)')x(nstar),x2,y(nstar),y2
        goto 111

100     continue
        close(20)
        nstar=nstar-1
        write(*,'(A,I5)')"      nstar = ",nstar

        call genWhisker(x,y,e1data,e2data,nstar,dataMask,0,scale,
     c                        outputData,avg)
        write(*,'(A,f12.6)')"      average e  ",avg

        call genWhisker(x,y,de1,de2,nstar,dataMask,0,scale,
     c                        outputDiff,avg)
        write(*,'(A,f12.6)')"      average de ",avg

        stop
        end
