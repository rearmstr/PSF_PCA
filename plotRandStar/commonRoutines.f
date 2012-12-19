c----------------------List of subroutines -----------------------------
c subroutine genWhisker(x,y,e1,e2,nstar,dataMask,ic,scale,outputFile,avg)
c subroutine numToChar(aNum,aChar,iChar)
c
c-----------------------------------------------------------------------
        subroutine genWhisker(x,y,e1,e2,nstar,dataMask,ic,scale,
     c                        outputFile,avg)
c
c take x,y and e1,e2, generate whisker plot data file in outputFile
c ic controls the use of dataMask (ic=0 no mask; ic=1 use mask;
c    ic=-1 reverse mask)
c scale is the magnification factor for the length of whisker
c avg is the average ellipticity
c-----------------------------------------------------------------------
        integer nstar,ic,i,ncount
        real*8 e1(nstar),e2(nstar),dataMask(nstar),avg,e,pi,scale,phi,
     c         dx,dy,x(nstar),y(nstar)
        character*32 outputFile

        pi=3.14159265
        if (ic .eq. 0) then
          ncount=nstar               ! ncount count unmasked stars
        else
          ncount=0
        endif

        open(23,file=outputFile)
        avg=0.0
        do i = 1, nstar
          e=dsqrt(e1(i)**2+e2(i)**2)
          if (ic .eq. 1) then
            e=e*dataMask(i)          ! e=0 when masked
            ncount=ncount+dataMask(i)
          endif
          if (ic .eq. -1) then
            e=e*(1-dataMask(i))     ! reverse mask
            ncount=ncount+1-dataMask(i)
          endif

          avg=avg+e

          if (dabs(e1(i)) .LT. 1.0E-8) then   ! avoid dividing by zero
            if (e2(i) .GE. 0.0) then
              phi=pi/4.0
            else
              phi=3.0*pi/4.0
            endif
          else                                ! calc phi
            phi=atan(e2(i)/e1(i)) / 2.0
          endif

          if (e1(i)*e2(i) .GT. 0.d0) then     ! adjust phi
            if (e1(i) .LT. 0.d0) phi=phi + pi/2.0
          else
            if (e1(i) .LT. 0.d0) then
              phi=phi + pi/2.0
            else
              phi=phi + pi
            endif
          endif

          dx=e*dcos(phi) * scale             ! whisker's x & y length
          dy=e*dsin(phi) * scale

          write(23,'(5E15.6)')x(i)-dx/2.0, y(i)-dy/2.0, dx, dy, e
        enddo
        close(23)

        avg=avg/ncount

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
