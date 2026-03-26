      program gc
C234567--------------------------------------------------------------012
C     Resolution d'un systeme lineaire par gradient conjugue
C     Attention la matrice doit etre symetrique definie positive
C234567--------------------------------------------------------------012
      parameter (nmax=10)
      double precision a(1:nmax,1:nmax), b(nmax), x(nmax)
 1    print *,'Nombre de lignes du systeme lineaire a resoudre n='
      read *,n
      if( n .le. 0 .or. n .gt. nmax ) goto 1
C
C     entree des coefficients de chaque ligne de la matrice A
C     et du second membre b du systeme lineaire
      CALL lecAxb( n, A, b )
C
C     resolution du systeme lineaire Ax=b
      CALL gcAxb( n, A, b, x, ierr )
C     Erreur?
      if( ierr .ne. 0 ) then
         print *,'Matrice non definie positive'
         stop
      endif
C
C     Affichage du resultat
      print *,'Le vecteur solution x est'
      print *,('  x(',i,')=',x(i), i=1,n )
      print *
      stop
      end
      SUBROUTINE lecAxb( n, A, b )
      double precision a(1:n, 1:n), b(n)
C     entree des coefficients de chaque ligne de la matrice A
C     et du second membre b
      do 20 i=1,n,1
         do 10 j=1,i
            print *,'A(',i,',',j,')='
            read *,a(i,j)
            a(j,i)=a(i,j)
 10      continue
 20   continue
      do 30 i=1,n
         print *,'b(',i,')='
         read *,b(i)
 30   continue
      return
      end
      subroutine gcAxb( n, A, b, x, ierr )
C     methode du gradient conjugue pour resoudre A x = b
      parameter (nmax=10)
      double precision A(n,n), b(n), x(n)
      double precision r(nmax), v(nmax), Av(nmax)
      double precision muk, lk1, r0r0, rkrk, rk1rk1, Avkvk
      double precision prosca
C
C     X0 = 0
      do 10 i=1,n
         x(i) = 0d0
 10   continue
C
C     R0 = V0 = b - A X0
      do 20 i=1,n
         r(i)=b(i)
         v(i)=b(i)
 20   continue
      r0r0 = prosca( n, r, r )
      rkrk = r0r0
      k=0
C
C     les iterations
 50   call matvec( n, A, v, Av )
      Avkvk = prosca( n, Av, v )
C     rkrk  = prosca( n, r, r )
      muk = rkrk / Avkvk
C     xk+1=xk + muk vk
      CALL cl2vec( n, x, muk, v, x )
C     rk+1 = rk - muk A vk
      CALL cl2vec( n, r, -muk, Av, r )
      rk1rk1 = prosca( n, r, r )
      print *,'iteration ',k+1,'  ||rk+1||**2=',rk1rk1
      lk1 = rk1rk1 / rkrk
C     vk+1 = rk+1 + lk1 vk
      CALL cl2vec( n, r, lk1, v, v )
C     test d'arret
      if( rk1rk1 .gt. 1d-4 * r0r0 ) then
         k = k + 1
         rkrk = rk1rk1
         if( k .gt. n ) then
            ierr = 1
            return
         endif
         goto 50
      endif
      ierr = 0
      return
      end
      double precision function prosca( n, u, v )
C     produit scalaire de 2 vecteurs de n composantes
      double precision u(n),v(n)
      prosca = 0d0
      do 10 i=1,n
         prosca = prosca + u(i) * v(i)
 10   continue
      return
      end
      SUBROUTINE cl2vec( n, u, l, v, w )
C     vecteurs de n composantes et w = u + l v
      double precision u(n), l, v(n), w(n)
      do 10 i=1,n
         w(i) = u(i) + l * v(i)
 10   continue
      return
      end
      SUBROUTINE matvec( n, A, u, v )
C     produit matrice vecteur de n composantes
      double precision A(n,n), u(n), v(n), s
      do 10 i=1,n
         s = 0d0
         do 5 j=1,n
            s = s + A(i,j) * u(j)
 5       continue
         v(i) = s
 10   continue
      return
      end
C-------------------------------------------------------------------------
C perronne@PC:~/COURS_FORTRAN -bash-> f77 gc.f -o gc
C perronne@PC:~/COURS_FORTRAN -bash-> gc
C Nombre de lignes du systeme lineaire a resoudre n=
C 1
C  A( 1, 1)=
C 5
C  b( 1)=
C 5
C  iteration  1  ||rk+1||**2=  7.70371978E-32
C  Le vecteur solution x est
C    x( 1)=  1.
C perronne@PC:~/COURS_FORTRAN -bash-> gc
C  Nombre de lignes du systeme lineaire a resoudre n=
C 2
C  A( 1, 1)=
C 3
C  A( 2, 1)=
C -1
C  A( 2, 2)=
C 2
C  b( 1)=
C 1
C  b( 2)=
C 3
C  iteration  1  ||rk+1||**2=  1.11111111
C  iteration  2  ||rk+1||**2=  3.50658348E-32
C  Le vecteur solution x est
C    x( 1)=  1.  x( 2)=  2.
C perronne@PC:~/COURS_FORTRAN -bash-> gc
C  Nombre de lignes du systeme lineaire a resoudre n=
C 3
C  A( 1, 1)=
C 5
C  A( 2, 1)=
C -1
C  A( 2, 2)=
C 4
C  A( 3, 1)=
C -2
C  A( 3, 2)=
C -3
C  A( 3, 3)=
C 6
C  b( 1)=
C -3
C  b( 2)=
C -2
C  b( 3)=
C 10
C  iteration  1  ||rk+1||**2=  8.14313298
C  iteration  2  ||rk+1||**2=  3.99079872
C  iteration  3  ||rk+1||**2=  6.21766331E-30
C  Le vecteur solution x est
C    x( 1)=  1.  x( 2)=  2.  x( 3)=  3.
C-------------------------------------------------------------------------
