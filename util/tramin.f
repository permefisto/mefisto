      subroutine tramin( nb, dtab )
C     dtab = dtab - min dtab
c
      double precision dtab(nb), min
C
      min = dtab(1)
      do i=2,nb
         if( dtab(i) .lt. min ) min=dtab(i)
      enddo
      print *,'tramin: min=',min
c
      do i=1,nb
         dtab(i)=dtab(i)-min
      enddo
c
      return
      end
