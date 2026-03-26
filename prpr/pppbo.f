      program pppbo
      character*5 str
      common /n/ a(4)
      common  str,m(1000000)
      real    r(1)
      equivalence (m(1),r(1))
C
      do i=1,100
         m(i) = i
      enddo
C
      a(1) = 11.
      a(2) = 22.
      a(3) = 33.
C
      call bar( 99 )
      call foo()
C
      print *,'pppbo:', a
      print *,'pppbo:', (m(k),k=1,15)
      print *,'pppbo:', (r(k),k=16,25)
      STOP
      end
