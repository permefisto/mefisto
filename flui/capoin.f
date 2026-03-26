      subroutine capoin(icpl,icpd,icpc,icplp,icpcp,ind,lca,niv)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                    s.p. capoin
c                    -----------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but : ce sp calcule les pointeurs icplp,icpcp
c ----  correspondants au 'remplissage' de la
c       matrice en cours de factorisation
c
c       version generale : matrice a structure quelconque
c
c   parametres d'entree:
c   -------------------
c   icpl,icpd,icpc : les pointeurs associes
c   ind            : liste des coefficients
c   lca            : longieur maximum du tableau ca
c   niv            : niveau de la factorisation
c
c   parametre de sortie:
c   -------------------
c   icplp,icpcp    : les pointeurs cherches
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   programmeur : P. Joly, Laboratoire Jacques-Louis Lions
c                 Universite Pierre et Marie Curie - CNRS (UMR 7598)
c   email       : joly@ann.jussieu.fr - phone  : (33) 1 44 27 44 10
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension icpl(0:n),icpd(n),icpc(nzca)
      dimension icplp(0:n),icpcp(nzca),ind(n)
      common/rang/n,nza,nzca
c
      do i=1,n
         ind(i)=0
      enddo
      ic1=0
      ic2=0
      k1=1
      do i=1,n
         k2=icpl(i)
         do k=k1,k2
            j=icpc(k)
            ind(j)=i
         enddo
         ind(i)=i
         do k=k1,icpd(i)
            j=icpc(k)
            do l=icpd(j)+1,icpl(j)
               jj=icpc(l)
               if(ind(jj).ne.i) then
                  ic1=ic1+1
                  if(ic1.gt.lca) then
                     print 1000,niv,lca
                     return
                  endif
                  icpcp(ic1)=jj
                  ind(jj)=i
               end if
            enddo
         enddo
         call tri(icpcp(ic2+1),ic1-ic2)
         icplp(i)=ic1
         ic2=ic1
         k1=k2+1
      enddo
      icplp(0)=0
      k1=icpl(n)
      kp1=icplp(n)
      nzero=k1+kp1
      if(nzero.gt.lca) then
         print 1000,niv,lca
         return
      endif
      k=nzero
      do i=n,1,-1
         icpl(i)=k
         kp2=icplp(i-1)
         k2=icpl(i-1)
2        continue
         if(k1.gt.k2) then
            if(kp1.gt.kp2) then
            if(icpc(k1).lt.icpcp(kp1)) then
               icpc(k)=icpcp(kp1)
               kp1=kp1-1
            else
               icpc(k)=icpc(k1)
               k1=k1-1
            end if
            else
               icpc(k)=icpc(k1)
               k1=k1-1
            end if
         else
            if(kp1.gt.kp2) then
               icpc(k)=icpcp(kp1)
               kp1=kp1-1
            else
               go to 1
            end if
         end if
          k=k-1
         go to 2
1        continue
      enddo
      nzca=nzero
c
      return
1000  format(' niveau',i4,' arret des calculs place memoire',i9,
     s       ' insuffisante')
      end
