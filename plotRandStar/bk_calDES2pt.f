c ----------- calculate the 2pt function of residual -------------------
c input x, y, e1, e2 from PCA reconstruction and data
c       this requires exposure id
c       calculate the differences in ellipticities
c       assign each pairs to bins in r 
c       take the average of e(r1)*e(r2) in each r bin
c number of r bins and dr is hardwired
c output the 2pt function (r,2pt)
c-----------------------------------------------------------------------
	program main

        integer nstarMax,nstar,i,j,IDexp,ichar
        parameter (nstarMax=100000,IDexp=100)
        real*8 x(nstarMax),y(nstarMax),
     c         e1data(nstarMax),e2data(nstarMax),
     c         de1(nstarMax),de2(nstarMax),x2,y2

        integer nChips,mc,nc,nShapeLet,im,in,iExp,iChip,ngrid,igrid,nr
        parameter (nChips=62,mc=20,nc=40,nShapeLet=2,nr=20)
        parameter (ngrid=nChips*mc*nc)

        integer npt(nr),ibin

        real*8 dr,rmin,rmax,r(nr),rval,twoPtFun(nr),error(nr),e1e2
        character*50 inputStarRec, inputE12recon, outputFile,
     c               aChar, inputMask

        real*8 e1recon(ngrid),e2recon(ngrid)
        real*8 dataMask(nstarMax),avg     ! dataMask is not used

        inputStarRec="starXYe12.dat"
        inputE12recon="/data/mzm/DES_20x40_400/reconEM"
        ! inputMask="../results/MPI/50pct/dataMask"

        outputFile="2ptFunc.dat"

        rmin=0.0
        rmax=2.5

        dr=(rmax-rmin)/nr
        do i = 1, nr
          r(i)=rmin+(i-1/2)*dr
          npt(i)=0
          twoPtFun(i)=0.0
        enddo

        open(21,file=inputE12recon)      ! read in reconstruction
        do i = 1, IDexp-1
          read(21,*)
        enddo
        read(21,*)(e1recon(i),i=1,ngrid), (e2recon(j),j=1,ngrid)
        close(21)

        do i = 1, ngrid
           e1recon(i)=e1recon(i)/sqrt(2.0)
           e2recon(i)=-e2recon(i)/sqrt(2.0)
        enddo

        open(20,file=inputStarRec)
        nstar=0
111     continue
          nstar=nstar+1
          read(20,*,end=100)iExp,iChip,im,in,x2,y2,
     c            e1data(nstar),e2data(nstar),x(nstar),y(nstar)
          e1data(nstar)=e1data(nstar)/sqrt(2.0)
          e2data(nstar)=-e2data(nstar)/sqrt(2.0)
          igrid=(iChip-1)*mc*nc + (in-1)*mc + im
          de1(nstar)=e1recon(igrid)-e1data(nstar)
          de2(nstar)=e2recon(igrid)-e2data(nstar)
        goto 111

100     continue
        close(20)
        nstar=nstar-1
        write(*,'(A,I5)')"      nstar = ",nstar

        outputFile="2ptFuncResidual.dat"
        call cal2pt(x,y,de1,de2,nstar,dataMask,0,
     c                    rmin,rmax,dr,r,nr,outputFile)

        outputFile="2ptFuncData.dat"
        call cal2pt(x,y,e1data,e2data,nstar,dataMask,0,
     c                    rmin,rmax,dr,r,nr,outputFile)

!        outputFile="2ptFuncRecon.dat"   ! recon has different x,y
!        call cal2pt(x,y,e1recon,e2recon,nstar,dataMask,0,
!     c                    rmin,rmax,dr,r,nr,outputFile)

        stop
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
            ibin=(rval-rmin)/dr
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
