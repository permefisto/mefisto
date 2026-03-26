      subroutine prosca(n,x,y,xpy)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                     sp prosca
c                     ---------
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but : calculer le produit scalaire de deux vecteurs
c ----
c
c parametres d'entree :
c ---------------------
c x et y : vecteurs dont le produit est a calculer
c
c parametre de sortie :
c ---------------------
c xpy    : le produit scalaire
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   programmeur : P. Joly, Laboratoire Jacques-Louis Lions
c                 Universite Pierre et Marie Curie - CNRS (UMR 7598)
c   email       : joly@ann.jussieu.fr - phone  : (33) 1 44 27 44 10
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension x(n),y(n)
c
      s=0.d0
      do i=1,n
         s=s+x(i)*y(i)
      enddo
      xpy=s
c
      return
      end
