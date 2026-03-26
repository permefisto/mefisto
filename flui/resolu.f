      subroutine resolu(ca,dca,icpl,icpd,icpc,z,r)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                    s.p. resolu
c                    -----------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but : ce sp descend le systeme triangulaire inferieur:
c -----                     (( l ca )) (z) = (r)
c       remonte le systeme triangulaire superieur:
c                           (( u ca )) (z) = (z)
c
c       version generale : matrice a structure quelconque
c
c   parametres d'entree:
c   -------------------
c   ca,dca         : la matrice de conditionnement
c   icpl,icpd,icpc : les pointeurs associes
c   r              : second membre
c
c   parametre de sortie:
c   -------------------
c   z              : vecteur solution
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   programmeur : P. Joly, Laboratoire Jacques-Louis Lions
c                 Universite Pierre et Marie Curie - CNRS (UMR 7598)
c   email       : joly@ann.jussieu.fr - phone  : (33) 1 44 27 44 10
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension ca(nzca),dca(n),z(n),r(n)
      dimension icpl(0:n),icpd(n),icpc(nzca)
      common/rang/n,nza,nzca
c
c
c     la descente  ( l ca )  (z) = (r)
c     --------------------------------
c
      do i=1,n
         s=r(i)
         do k=icpl(i-1)+1,icpd(i)
            j=icpc(k)
            s=s-ca(k)*z(j)
         enddo
         z(i)=s
      enddo
c
c     la remontee   ( u ca )  (z) = (r)
c     ---------------------------------
c
      do i=n,1,-1
         s=z(i)
         do k=icpd(i)+1,icpl(i)
            j=icpc(k)
            s=s-ca(k)*z(j)
         enddo
         z(i)=s*dca(i)
      enddo
c
      return
      end
