c-----------------------------------------------------------------------
c bin the input data
c     icol=1: histogram data
c     icol=2: histogram missing data
c Output is histo.dat
c changed from fixed number of bins to fixed binSize
c-----------------------------------------------------------------------
	program main
        integer NgMax,Ng,i,j,nzbinMax,nzbin,icol,imin,imax,imed
        parameter (NgMax=10**6,nzbinMax=3000,icol=1)
        integer ncount(nzbinMax),nSum
        real*8 zs(NgMax),tempVal,zmin,zmax,dz,zmean,zmed

        parameter (dz=5.0E-5)

        character*100 inputFile,outputFile

        inputFile="SNAPreconErr.dat"
        outputFile="histo.dat"

        open(20,file=inputFile)
        do i = 1, NgMax
          read(20,*,end=50)(tempVal, j=1,icol-1), zs(i)
        enddo
        write(6,*)"    NgMax is exceeded!!!"
        stop

50      close(20)
        Ng=i-1

        zmin=100.d0
        zmax=0.d0
        zmean=0.d0
        do i = 1, Ng
          zmean=zmean+zs(i)
          if (zs(i) .lt. zmin) then
            zmin=zs(i)
            imin=i
          endif
          if (zs(i) .gt. zmax) then
            zmax=zs(i)
            imax=i
          endif
        enddo
        write(6,'(2(A,f12.6),A)')"   z range = {",zmin,",",zmax,"}"
        write(6,'(2(A,I6))')"   Imin = ",imin,"  Imax = ",imax
        write(6,'(A,f12.6)')"   z mean = ", zmean/Ng

        nzbin=(zmax-zmin)/dz

        do i = 1, nzbin
          ncount(i)=0
        enddo

        do i = 1, Ng
          j=INT((zs(i)-zmin)/dz)+1
          if (j .lt. nzbin) then
            ncount(j)=ncount(j)+1
          else
            ncount(nzbin)=ncount(nzbin)+1
          endif
        enddo

        open(21,file=outputFile)
        do i = 1, nzbin
          write(21,101)zmin+(i-0.5)*dz,ncount(i)
        enddo
        close(21)
101     format(1x,f10.5,2x,I8)

        nSum=0
        do i = 1, nzbin
          nSum=nSum+ncount(i)
          if (nSum .GE. Ng/2) goto 102
        enddo
102     continue
        imed=i-1
        write(6,'(A,f12.6)')"   z median = ", zmin+(imed-1)*dz

        end 
c-----------------------------------------------------------------------
