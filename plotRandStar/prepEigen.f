c ----------- prepare data for PSF whisker plot ------------------------
c input x, y, e1, e2 of the eigen PSF pattern
c output x-dx/2, y-dy/2, dx, dy
c gnuplot plots the output by
c   plot 'file.dat' using 1:2:3:4 with vectors head filled lt 2
c-----------------------------------------------------------------------
	program main

        integer nstar,i,IDeigen,ichar
        ! parameter (nstar=900,IDeigen=4)
        parameter (nstarMax=100000,IDeigen=4)
        real*8 x(nstarMax),y(nstarMax),e1(nstarMax),e2(nstarMax),
     c         dx,dy,phi,pi,scale,x2,y2
        character*50 inputXY, inputE12, outputFile, aChar

        integer nstarP
        real*8 dataMask(nstarMax),avg

        call numToChar(IDeigen,aChar,ichar)

        inputXY="../results/gridXY2"
        inputE12="../results/eigenVecSVD"
        ! inputXY="../results/new_BS/randStar/gridXY2"
        ! inputE12="../results/new_BS/randStar/eigenVecEM"

        outputFile="whiskerEigen.dat"

        pi=3.14159265
        scale=0.2

        open(20,file=inputXY)
        do i = 1, nstarMax
          read(20,*,end=100)x(i),y(i),x2,y2
          x(i)=(x(i)+x2)/2.0
          y(i)=(y(i)+y2)/2.0
        enddo
        close(20)
        write(*,'(A)')"Error: nstarMax is reached!"
        stop

100     continue
        close(20)
        nstar=i-1
        write(*,'(A,I5)')"      nstar = ",nstar

        open(21,file=inputE12)
        do i = 1, IDeigen-1
          read(21,*)
        enddo
        read(21,*)(e1(i),i=1,nstar), (e2(j),j=1,nstar)
        close(21)

        nstarP=nstar
        call genWhisker(x,y,e1,e2,nstarP,dataMask,0,scale,
     c                        outputFile,avg)
        write(*,'(A,f12.6)')"      average e ",avg

        stop
        end

