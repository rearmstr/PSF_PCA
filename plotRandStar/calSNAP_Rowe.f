c ----------- calculate the 2pt function of residual -------------------
c input x, y, e1, e2 from PCA reconstruction and data
c       this requires exposure id
c       calculate the differences in ellipticities
c       assign each pairs to bins in r 
c       calculate Rowe 2009 correlation functions D1 & D2
c number of r bins and dr is hardwired
c output the 2pt function (r,2pt)
c-----------------------------------------------------------------------
	program main

        integer nstarMax,nstar,i,j,IDexp,NumExp,ichar,nOutputFileLen
        parameter (nstarMax=4200,NumExp=2000)

        real*8 x(nstarMax),y(nstarMax),
     c         e1data(nstarMax),e2data(nstarMax),e1dataVal,e2dataVal,
     c         de1(nstarMax),de2(nstarMax),x2,y2

        integer nChips,mc,nc,nShapeLet,im,in,iExp,iChip,ngrid,igrid,nr
        parameter (nChips=44,mc=10,nc=10,nShapeLet=2,nr=100)
        parameter (ngrid=nChips*mc*nc)

        integer ibin,nDataDir,nReconDir,valFlag(nstarMax),iFlag

        real*8 dr,rmin,rmax,r(nr),rval,twoPtFun(nr),error(nr),e1e2,
     c         D1(nr),D2(nr),xi_p(nr),D1err(nr),D2err(nr),xiErr(nr),
     c         npt(nr)
        character*100 inputStarRec, inputE12recon, outputFile,
     c                inputData,aChar,dataDir,reconDir

        real*8 e1recon(ngrid),e2recon(ngrid)
        real*8 dataMask(nstarMax),avg     ! dataMask is not used

        integer icLog,icLogP,nrP       ! =0 linear; =1 log
        parameter (icLog=1)

        dataDir=
     c    "/astro/tutti1/mzm/SNAPdata/new_allMag_4200_randXY_one8thN/"
        reconDir="/data/mzm/SNAP_10x10_2k_1-8thN_50pctValSet/"
        call fileNameLen(dataDir,nDataDir)
        call fileNameLen(reconDir,nReconDir)

        inputStarRec=reconDir(1:nReconDir)//"starRec"
        inputE12recon=reconDir(1:nReconDir)//"reconEM"

        ! outputFile="D1_D2_SNAP_valSet_2k.dat"
        outputFile="D1_D2_SNAP_2k_10x10_1-8thN_valSet.dat"
        call fileNameLen(outputFile,nOutputFileLen)

        icValidate=0          ! 0: validation set
        icLogP=icLog
        nrP=nr

        rmax=0.6              ! SNAP focal plane is 0.6m in diameter
        if (icLog .eq. 0) then
           rmin=0.0
           dr=(rmax-rmin)/nr
        else
           rmin=1.0E-5        ! close to the PSF size?
           dr=dlog10(rmax/rmin)/nr
        endif

        do i = 1, nr
          if (icLog .eq. 0) then
             r(i)=rmin+(i-1/2)*dr
          else
             r(i)=rmin*10**((i-1/2)*dr)
          endif
          npt(i)=0.0
          xi_p(i)=0.0
          D1(i)=0.0
          D2(i)=0.0
          xiErr(i)=0.0
          D1Err(i)=0.0
          D2Err(i)=0.0
        enddo

        open(20,file=inputStarRec)
        open(21,file=inputE12recon)

        do IDexp = 1, NumExp

           if (mod(IDexp,100) .eq. 0) write(*,*)"    exposure ",IDexp
           read(21,*)(e1recon(i),i=1,ngrid), (e2recon(j),j=1,ngrid)

           nstar=0
           do istar = 1, nstarMax
              read(20,*)iExp,iChip,im,in,x2,y2,
     c                  e1dataVal,e2dataVal,iFlag
              if (iChip .ne. -998) then
                nstar=nstar+1
                x(nstar)=x2
                y(nstar)=y2
                e1data(nstar)=e1dataVal
                e2data(nstar)=e2dataVal
                valFlag(nstar)=iFlag
                igrid=(iChip-1)*mc*nc + (in-1)*mc + im
                if (igrid .lt. 1 .OR. igrid .gt. ngrid) then
                   write(*,'(6I8)')iExp,iChip,im,in,igrid,ngrid
                endif
                de1(nstar)=e1recon(igrid)-e1data(nstar)
                de2(nstar)=e2recon(igrid)-e2data(nstar)
              endif
           enddo

           call calRoweStat(x,y,e1data,e2data,de1,de2,nstar,
     c                      valFlag,icValidate,
     c                      rmin,rmax,dr,r,nrP,outputFile,icLogP,
     c                      npt,D1,D2,xi_p,D1err,D2err,xiErr)

        enddo         ! end of looping over exposures

        close(20)
        close(21)

        open(22,file=outputFile(1:nOutputFileLen))
        do i = 1, nr
          if (npt(i) .GT. 0.0) then
            xi_p(i)=xi_p(i)/npt(i)
            D1(i)=D1(i)/npt(i)
            D2(i)=D2(i)/npt(i)
            xiErr(i)=sqrt(xiErr(i)/npt(i)-xi_p(i)**2)/sqrt(npt(i))
            D1err(i)=sqrt(D1err(i)/npt(i)-D1(i)**2) / sqrt(npt(i))
            D2err(i)=sqrt(D2err(i)/npt(i)-D2(i)**2) / sqrt(npt(i))
            write(22,'(8E15.6)')r(i),xi_p(i),xiErr(i),D1(i),
     c                           D1err(i),D2(i),D2err(i),npt(i)
          endif
        enddo
        close(22)

        stop
        end


c-----------------------------------------------------------------------
        subroutine calRoweStat(x,y,e1,e2,de1,de2,nstar,valFlag,ic,
     c                         rmin,rmax,dr,r,nr,outputFile,icLog,
     c                         npt,D1,D2,xi_p,D1err,D2err,xiErr)
c
c take x,y, e1,e2, dataMask, binning in r, calculate D1, D2 in Barnaby
c Rowe paper 0904.3056 and output to outputFile.
c ic controls the use of validation set (ic=1 no; ic=0 validation set)
c-----------------------------------------------------------------------
        integer nstar,ic,i,j,ncount,nr,ibin,icLog
        integer valFlag(nstar)
        real*8 e1(nstar),e2(nstar),de1(nstar),de2(nstar),
     c         dx,dy,x(nstar),y(nstar),npt(nr),
     c         rmin,rmax,dr,r(nr),rval,eAmp,phi,tempVal,
     c         eTanV_1,eXv_1,deTanV_1,deXv_1,
     c         eTanV_2,eXv_2,deTanV_2,deXv_2,
     c         D1(nr),D2(nr),xi_p(nr),D1err(nr),D2err(nr),xiErr(nr)

        character*100 outputFile


        do i = 1, nstar-1
          if (valFlag(i) .eq. ic) then
            do j = i+1, nstar
              if (valFlag(j) .eq. ic) then

                rval=dsqrt((x(i)-x(j))**2 + (y(i)-y(j))**2)
                if (icLog .eq. 0) ibin=(rval-rmin)/dr + 1
                if (icLog .eq. 1) ibin=dlog10(rval/rmin)/dr + 1

                if (ibin .gt. nr) then
                   write(*,*)"   ibin = nr at ",i,rval
                   ibin=nr
                endif
                if (ibin .lt. 1) then
                   write(*,'(2I6,3E15.6,I8)')i,j,rval,rmin,dr,ibin
                   ibin=1
                endif

                npt(ibin)=npt(ibin)+1.0
                phi=dacos((y(i)-y(j))/rval)
                phi=phi*2.0

                eTanV_1=-e1(i)*dcos(phi)-e2(i)*dsin(phi)
                eXv_1 =  e1(i)*dsin(phi)-e2(i)*dcos(phi)

                deTanV_1=-de1(i)*dcos(phi)-de2(i)*dsin(phi)
                deXv_1 =  de1(i)*dsin(phi)-de2(i)*dcos(phi)

                eTanV_2=-e1(j)*dcos(phi)-e2(j)*dsin(phi)
                eXv_2 =  e1(j)*dsin(phi)-e2(j)*dcos(phi)

                deTanV_2=-de1(j)*dcos(phi)-de2(j)*dsin(phi)
                deXv_2 =  de1(j)*dsin(phi)-de2(j)*dcos(phi)

                tempVal = eTanV_1*eTanV_2 + eXv_1*eXv_2
                xi_p(ibin)=xi_p(ibin) + tempVal
                xiErr(ibin)=xiErr(ibin) + tempVal**2

                tempVal = deTanV_1*deTanV_2 + deXv_1*deXv_2
                D1(ibin)=D1(ibin) + tempVal
                D1Err(ibin)=D1Err(ibin) + tempVal**2

                tempVal = eTanV_1*deTanV_2 + eXv_1*deXv_2
     c                  + deTanV_1*eTanV_2 + deXv_1*eXv_2
                D2(ibin)=D2(ibin) + tempVal
                D2Err(ibin)=D2Err(ibin) + tempVal**2
              endif
            enddo
          endif
        enddo

        return
        end



c-----------------------------------------------------------------------
        subroutine numToChar(aNum,aChar,iChar)
c
c Convert a number to a character and
c determine the length of the character.
c-----------------------------------------------------------------------
        integer aNum,i,iChar
        character*32 aChar

        open(unit=25,file="fileTable")
        write(25,*)aNum
        close(unit=25)

        open(unit=25,file="fileTable")
        read(25,*)aChar
        close(unit=25)

        do i = 1, 32
          if(aChar(i:i) .EQ. " ") goto 90
        enddo
90      iChar=i-1

        return
        end


c-----------------------------------------------------------------------
        subroutine fileNameLen(fileName,nameLen)
c count the number of characters in fileName.
c-----------------------------------------------------------------------
        integer nameLen,i
        character*100 fileName

        do i = 1, 100
          if(fileName(i:i) .EQ. " ") goto 50
        enddo
        write(*,'(a)') 'fileName has more than 100 characters!'
        stop
50      nameLen = i - 1

        return
        end
c-----------------------------------------------------------------------
