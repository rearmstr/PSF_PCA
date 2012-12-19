! This is used to arrange the 8 BCS CCD chips on the focal plane.
!            1 2 3 4
!            5 6 7 8
! Shift the pixel values in gridXY2 to put CCDs into place.
!
	program main

        integer nChips, Nx, Ny, dx, dy, nCell
        parameter (nChips = 8)  ! total number of CCDs
        parameter (Nx = 10)     ! number of grid along x
        parameter (Ny = 20)     ! number of grid along y
        parameter (dx = 100)    ! chip separation along x
        parameter (dy = 100)    ! chip separation along y

        integer xShift(nChips), yShift(nChips)
        real*8 x1, y1, x2, y2

        character*50 inputXY, outputXY

        nCell = Nx*Ny           ! cells per CCD chip
        inputXY="gridXY2_BCS"
        outputXY="gridXY2_BCS_shifted"

        xShift(5) = 0
        yShift(5) = 4096 + dy

        xShift(6) = 2048 + dx
        yShift(6) = 4096 + dy

        xShift(7) = 2*(2048 + dx)
        yShift(7) = 4096 + dy

        xShift(8) = 3*(2048 + dx)
        yShift(8) = 4096 + dy

        xShift(1) = 0
        yShift(1) = 0

        xShift(2) = 2048 + dx
        yShift(2) = 0

        xShift(3) = 2*(2048 + dx)
        yShift(3) = 0

        xShift(4) = 3*(2048 + dx)
        yShift(4) = 0

        open(20, file=inputXY)
        open(21, file=outputXY)

        do i = 1, nChips
          do j = 1, nCell
            read(20,*)x1,y1,x2,y2
            write(21,100)x1+xShift(i),y1+yShift(i),
     c                   x2+xShift(i),y2+yShift(i)
          enddo
        enddo
100     format(1x,4f12.4)

        close(20)
        close(21)
c
c generate LL, UL, UR, LR corner of the chips
c
        open(22, file="LL.dat")
        x1=0
        y1=0
        do i = 1, nChips
          write(22,101)i,x1,y1,x1+xShift(i),y1+yShift(i)
        enddo
        close(22)
101     format(1x,I4,4f12.4)

        open(22, file="UL.dat")
        x1=0
        y1=4096
        do i = 1, nChips
          write(22,101)i,x1,y1,x1+xShift(i),y1+yShift(i)
        enddo
        close(22)

        open(22, file="UR.dat")
        x1=2048
        y1=4096
        do i = 1, nChips
          write(22,101)i,x1,y1,x1+xShift(i),y1+yShift(i)
        enddo
        close(22)

        open(22, file="LR.dat")
        x1=2048
        y1=0
        do i = 1, nChips
          write(22,101)i,x1,y1,x1+xShift(i),y1+yShift(i)
        enddo
        close(22)

        end
