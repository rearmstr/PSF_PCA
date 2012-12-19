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

        integer nstarMax,nstar,i,j,IDexp,NumExp,ichar
        parameter (nstarMax=100000,NumExp=1)
        ! parameter (nstarMax=100000,IDexp=200)

        real*8 x(nstarMax),y(nstarMax),
     c         e1data(nstarMax),e2data(nstarMax),
     c         de1(nstarMax),de2(nstarMax),x2,y2

        integer nChips,mc,nc,nShapeLet,im,in,iExp,iChip,ngrid,igrid,nr
        parameter (nChips=62,mc=20,nc=40,nShapeLet=2,nr=100)
        parameter (ngrid=nChips*mc*nc)

        integer npt(nr),ibin

        real*8 dr,rmin,rmax,r(nr),rval,twoPtFun(nr),error(nr),e1e2,
     c         D1(nr),D2(nr),xi_p(nr),D1err(nr),D2err(nr),xiErr(nr)
        character*50 inputStarRec, inputE12recon, outputFile,
     c               aChar, inputMask, inputXY

        real*8 e1recon(ngrid),e2recon(ngrid)
        real*8 dataMask(nstarMax),avg     ! dataMask is not used

        inputStarRec="/data/mzm/DES_20x40_400/starRec"
        inputXY="/data/mzm/DES_20x40_400/ra_dec.dat"
        ! inputStarRec="starXYe12.dat"
        inputE12recon="/data/mzm/DES_20x40_400/reconEM"
        ! inputMask="../results/MPI/50pct/dataMask"

        outputFile="D1_D2.dat"

        rmin=0.0
        rmax=2.5

        dr=(rmax-rmin)/nr
        do i = 1, nr
          r(i)=rmin+(i-1/2)*dr
          npt(i)=0
          xi_p(i)=0.0
          D1(i)=0.0
          D2(i)=0.0
          xiErr(i)=0.0
          D1Err(i)=0.0
          D2Err(i)=0.0
        enddo

        open(19,file="Nstar")
        open(20,file=inputStarRec)
        open(21,file=inputE12recon)
        open(22,file=inputXY)

        do IDexp = 1, NumExp
           read(21,*)(e1recon(i),i=1,ngrid), (e2recon(j),j=1,ngrid)

           do i = 1, ngrid
              e1recon(i)=e1recon(i)/sqrt(2.0)
              e2recon(i)=-e2recon(i)/sqrt(2.0)
           enddo

           read(19,*)nstar
           do istar = 1, nstar
              read(20,*)iExp,iChip,im,in,x2,y2,
     c                  e1data(istar),e2data(istar)
              read(22,*)x(istar),y(istar)
              e1data(istar)=e1data(istar)/sqrt(2.0)
              e2data(istar)=-e2data(istar)/sqrt(2.0)
              igrid=(iChip-1)*mc*nc + (in-1)*mc + im
              de1(istar)=e1recon(igrid)-e1data(istar)
              de2(istar)=e2recon(igrid)-e2data(istar)
           enddo

           call calRoweStat(x,y,e1data,e2data,de1,de2,nstar,dataMask,0,
     c                      rmin,rmax,dr,r,nr,outputFile,
     c                      npt,D1,D2,xi_p,D1err,D2err,xiErr)

        enddo         ! end of looping over exposures

        close(19)
        close(20)
        close(21)
        close(22)

        open(22,file=outputFile)
        do i = 1, nr
          if (npt(i) .GT. 0) then
            xi_p(i)=xi_p(i)/npt(i)
            D1(i)=D1(i)/npt(i)
            D2(i)=D2(i)/npt(i)
            xiErr(i)=sqrt(xiErr(i)/npt(i)-xi_p(i)**2)/sqrt(real(npt(i)))
            D1err(i)=sqrt(D1err(i)/npt(i)-D1(i)**2) / sqrt(real(npt(i)))
            D2err(i)=sqrt(D2err(i)/npt(i)-D2(i)**2) / sqrt(real(npt(i)))
          endif
          write(22,'(10E15.6)')r(i),xi_p(i),xiErr(i),D1(i),D1err(i),
     c                         D2(i),D2err(i)
        enddo
        close(22)

        stop
        end


c-----------------------------------------------------------------------
        subroutine calRoweStat(x,y,e1,e2,de1,de2,nstar,dataMask,ic,
     c                         rmin,rmax,dr,r,nr,outputFile,
     c                         npt,D1,D2,xi_p,D1err,D2err,xiErr)
c
c take x,y, e1,e2, dataMask, binning in r, calculate D1, D2 in Barnaby
c Rowe paper 0904.3056 and output to outputFile.
c ic controls the use of dataMask (ic=0 no mask; ic=1 use mask;
c    ic=-1 reverse mask)
c-----------------------------------------------------------------------
        integer nstar,ic,i,j,ncount,nr,npt(nr),ibin
        real*8 e1(nstar),e2(nstar),de1(nstar),de2(nstar),
     c         dx,dy,x(nstar),y(nstar),dataMask(nstar),
     c         rmin,rmax,dr,r(nr),rval,eAmp,phi,tempVal,
     c         eTanV_1,eXv_1,deTanV_1,deXv_1,
     c         eTanV_2,eXv_2,deTanV_2,deXv_2,
     c         D1(nr),D2(nr),xi_p(nr),D1err(nr),D2err(nr),xiErr(nr)

        character*32 outputFile


        do i = 1, nstar-1
          do j = i+1, nstar
            rval=dsqrt((x(i)-x(j))**2 + (y(i)-y(j))**2)
            ibin=(rval-rmin)/dr + 1
            if (ibin .gt. nr) then
               write(*,*)"   ibin = nr at ",i,rval
               ibin=nr
            endif
            npt(ibin)=npt(ibin)+1
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
     c              + deTanV_1*eTanV_2 + deXv_1*eXv_2
            D2(ibin)=D2(ibin) + tempVal
            D2Err(ibin)=D2Err(ibin) + tempVal**2
          enddo
        enddo

        return
        end


c-----------------------------------------------------------------------
        subroutine cal2pt(x,y,e1,e2,nstar,dataMask,ic,
     c                    rmin,rmax,dr,r,nr,outputFile)
c
c take x,y, e1,e2, dataMask, binning in r, calculate 2pt func
c and output to outputFile.
c ic controls the use of dataMask (ic=0 no mask; ic=1 use mask;
c    ic=-1 reverse mask)
c-----------------------------------------------------------------------
        integer nstar,ic,i,j,ncount,nr,npt(nr),ibin
        real*8 e1(nstar),e2(nstar),dataMask(nstar),
     c         dx,dy,x(nstar),y(nstar),
     c         rmin,rmax,dr,r(nr),rval,eAmp,twoPtFun(nr),error(nr)
        character*32 outputFile

        do i = 1, nr
          npt(i)=0
          twoPtFun(i)=0.0
        enddo

        do i = 1, nstar-1
          do j = i+1, nstar
            rval=dsqrt((x(i)-x(j))**2 + (y(i)-y(j))**2)
            ibin=(rval-rmin)/dr + 1
            if (ibin .gt. nr) ibin=nr
            npt(ibin)=npt(ibin)+1
            eAmp=dsqrt((e1(i)**2+e2(i)**2)*(e1(j)**2+e2(j)**2))
            twoPtFun(ibin)=twoPtFun(ibin)+
     c                     (e1(i)*e1(j)+e2(i)*e2(j)) ! /eAmp
!     c                     e1(i)*e1(j)
          enddo
        enddo

        open(20,file=outputFile)
        do i = 1, nr
          if (npt(i) .GT. 0) then
            twoPtFun(i)=twoPtFun(i)/npt(i)
            error(i)=sqrt(1.0/npt(i))           ! Poisson error
          else
            error(i)=0.0
          endif
          write(20,'(3E15.6)')r(i),twoPtFun(i),error(i)
        enddo
        close(20)

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
