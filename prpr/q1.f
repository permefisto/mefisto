      subroutine q1quad( xyrect, stquad,
     %                   xyr,
     %                   xy, ierr )
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but :   calcul des 2 coordonnees xyr(2) dans xyrect=[xmin,xmax]x[ymin,ymax]
c -----   image par q1: Rectangle-->stquad (appartenant a Q1**2)
c         par une resolution directe due a nicolas Thenault
c
c entrees:
c --------
c xyrect : 1:2 abscisse-ordonnee 1:2 min max du rectangle
c stquad : les 2 coordonnees des 4 sommets du quadrangle image
c xyr    : 2 coordonnees dans Rectangle d'image a calculer
c
c sorties:
c --------
c xy     : 2vcoordonnees du point image dans le quadrangle stquad
c ierr   : 0 si calcul sans erreur, 1 si rectangle degenere
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c auteur : Perronnet alain   LaboratoireJ.-L. Lions UPMC Paris mars 2006
c234567..............................................................012
      double precision  xyrect(1:2,1:2), stquad(1:2,1:4), xyr(2), xy(2)
      double precision  dist(2), den, a, b, c, d
c
      do 10 i=1,2
         dist(i) = xyrect(i,2) - xyrect(i,1)
         if( dist(i) .le. 0d0 ) then
            ierr = 1
            xy(1) = 0d0
            xy(2) = 0d0
            return
         endif
 10   continue
c
c     le denominateur
      den = dist(1) * dist(2)
c
c     les autres coefficients
      a = xyrect(1,2) - xyr(1)
      b = xyr(1)      - xyrect(1,1)
      c = xyrect(2,2) - xyr(2)
      d = xyr(2)      - xyrect(2,1)
      do 20 i=1,2
         xy(i) = ( stquad(i,1) * a * c
     %           + stquad(i,2) * b * c
     %           + stquad(i,3) * b * d
     %           + stquad(i,4) * a * d
     %           ) / den
 20   continue
      ierr = 0
c
      return
      end


      subroutine q1qinv( xyrect, stquad,
     %                   xy,
     %                   xyr, ierr )
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but :   calcul des 2 coordonnees (xr,yr) dans Rectangle=[xmin,xmax]x[ymin,ymax]
c -----   image par q1: Rectangle-->stquad (appartenant a Q1**2)
c         par une resolution directe due a nicolas Thenault
c
c entrees:
c --------
c xyrect : 1:2 abscisses 1:2 min max du rectangle
c stquad : les 2 coordonnees des 4 sommets du quadrangle image
c xy     : 2 coordonnees dans le quadrangle stquad
c
c sorties:
c --------
c xyr    : 2 coordonnees du point dans le rectangle
c ierr   : 0 si calcul sans erreur, 1 si quadrangle degenere
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c auteurs: Thenault Tulenew  analyse numerique paris        janvier 1998
c modifs : Perronnet alain   LaboratoireJ.-L. Lions UPMC Paris mars 2006
c234567..............................................................012
      double precision zero
      parameter       (zero=1d-10)
c     ATTENTION: UNE MAUVAISE VALEUR TROP PETITE => ERREURS dans Q1-1!
c
      double precision xy(2),xyr(2),xyrect(2,2),stquad(2,4)
      double precision dist(2),xr,yr
      double precision a,b,c,d,alpha,beta,gamma,delta,x0,y0,t(2),u,v,w
c
      ierr = 0
c
      a = stquad(1,1)
      b = stquad(1,2) - stquad(1,1)
      c = stquad(1,4) - stquad(1,1)
      d = stquad(1,1) - stquad(1,2) + stquad(1,3) - stquad(1,4)
c
      alpha = stquad(2,1)
      beta  = stquad(2,2) - stquad(2,1)
      gamma = stquad(2,4) - stquad(2,1)
      delta = stquad(2,1) - stquad(2,2) + stquad(2,3) - stquad(2,4)
c
      u = beta  * c - b * gamma
      if( abs(u) .le. zero ) then
c        quadrangle degenere
         ierr = 1
         xyr(1) = xyrect(1,1)
         xyr(2) = xyrect(2,1)
         return
      endif
      v = delta * c - d * gamma
      w = b * delta - beta * d
c
      x0 = c * (xy(2)-alpha) - gamma * (xy(1)-a)
      y0 = b * (xy(2)-alpha) - beta  * (xy(1)-a)
c
      a = v  * w
      b = u  * u - w * x0 - v * y0
      c = x0 * y0
c
      if( abs(a) .gt. zero ) then
c
         delta = sqrt( b*b-4d0*a*c )
         if( b .ge. 0d0 ) then
            t(2) = -b - delta
         else
            t(2) = -b + delta
         endif
c        la racine de plus grande valeur absolue
c       (elle donne le plus souvent le point exterieur au carre unite
c        donc a tester en second pour reduire les calculs)
         t(2) = t(2) / ( 2d0 * a )
c        calcul de la seconde racine a partir de la somme => plus stable
         t(1) = - b/a - t(2)
c
         do 10 i=1,2
c
c           la solution i donne t elle un point interne au carre unite?
            xr = ( x0 - v * t(i) ) / u
            yr = ( w * t(i) - y0 ) / u
c
            if( 0d0 .le. xr .and. xr .le. 1d0 ) then
               if( 0d0 .le. yr .and. yr .le. 1d0 ) goto 8000
            endif
c
c           le point (xr,yr) n'est pas dans le carre unite
c           cela peut etre du aux erreurs d'arrondi
c           => choix par le minimum de la distance aux bords du carre
            dist(i) = max( 0d0, -xr, xr-1d0, -yr, yr-1d0 )
c
 10      continue
c
         if( dist(1) .gt. dist(2) ) then
c           Q1(xr,yr) pour la racine 2 est plus proche de xy(1),xy(2)
c           xr yr sont deja calcules
            goto 8000
         endif
c
      else if ( abs(b) .gt. zero ) then
         t(1) = - c / b
      else
         t(1) = 0d0
      endif
c
c     les 2 coordonnees du point dans le carre unite
      xr = ( x0 - v * t(1) ) / u
      yr = ( w * t(1) - y0 ) / u
c
c     les 2 coordonnees du point dans le rectangle [xmin,xmax]x[ymin,ymax]
 8000 xyr(1) = xr * ( xyrect(1,2) - xyrect(1,1) ) + xyrect(1,1)
      xyr(2) = yr * ( xyrect(2,2) - xyrect(2,1) ) + xyrect(2,1)
c
      return
      end


      subroutine sotr3l( longcote, x, y, ierr )
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but :    calcul des 2 coordonnees du sommet superieur du triangle
c -----    defini par la longueur de ses 3 cotes
c
c entrees:
c --------
c longcote : longueur des 3 aretes des 3 cotes du triangle
c          longcote(1) >= longcote(2) et longcote(1) >= longcote(3)
c
c sorties:
c --------
c x, y   : coordonnees du sommet superieur du triangle
c          de sommet 1 l'origine, de sommet 2 (longcote(1),0)
c ierr   : 0 si le triangle est constructible
c          1 sinon
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c auteur : alain perronnet  analyse numerique paris       septembre 1993
c234567..............................................................012
      double precision longcote( 3 ), x, y
c
      x = (longcote(1)**2+longcote(3)**2-longcote(2)**2)
     %  / (2.0*longcote(1))
      y = longcote(3)**2 - x**2
      if( y .le. 0 ) then
c        impossibilite de construire le triangle
         ierr = 1
      else
         y = sqrt( longcote(3)**2 - x**2 )
         ierr = 0
      endif
      end

      subroutine soqu4l( longcote, x3, y3, x4, y4, ierr )
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but :  calcul des 2 coordonnees des 2 sommets superieurs du quadrangle
c -----  defini par la longueur de ses 4 cotes
c
c entrees:
c --------
c longcote: longueur des 3 aretes des 4 cotes du quadrangle
c           longcote(1) >= longcote(2),longcote(3),longcote(4)
c           longcote(1) <  longcote(2)+longcote(3)+longcote(4)
c
c sorties:
c --------
c x3,y3, x4,y4 : coordonnees des sommets superieurs du quadrangle
c                de sommet 1 l'origine, de sommet 2 (longcote(1),0)
c ierr   : 0 si x3,y3, x4,y4 sont calcules normalement; >0 sinon
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c auteur : alain Perronnet analyse numerique Paris          octobre 1993
c modifs : alain Perronnet Laboratoire J.-L. Lions UPMC Paris  mars 2006
c234567..............................................................012
      parameter (nbp=10)
      double precision longcote(4), x3,y3, x4,y4
      double precision c3, c32, a, anginc, cosa, sina, d
c
c     verification sur l'existence d'une solution
      a = 1d-3 * longcote(1)
      if( longcote(1) .ge. longcote(2)+longcote(3)+longcote(4) ) then
c        cote 1 trop long pour les 3 autres cotes
c        le quadrangle est le rectangle de cotes 1 3 longcote(1)
c        et de cotes 2 4 de max(longcote(2), longcote(4))
         ierr = 1
         x3 = longcote(1)
         x4 = 0d0
         d = max( longcote(2), longcote(4) )
         y3 = d
         y4 = d
         return
      else if( longcote(2) .le. a ) then
         ierr = 2
         x3 = longcote(1)
         y3 = a
         x4 = 0d0
         y4 = longcote(4)
         return
      else if( longcote(3) .le. a ) then
         c3 = longcote(3)
         longcote(3) = longcote(4)
         call sotr3l( longcote, x3, y3, ierr )
         longcote(3) = c3
         if( ierr .eq. 0 ) then
c           triangle S1 S2 S3
            ierr = 3
            x4 = x3 - a
            y4 = y3 - a
         else
c           le quadrangle est le rectangle de cotes 1 3 longcote(1)
c           et de cotes 2 4 de max(longcote(2), longcote(4))
            ierr = 5
            x3 = longcote(1)
            x4 = 0d0
            d = max( longcote(2), longcote(4) )
            y3 = d
            y4 = d
         endif
         return
      else if( longcote(4) .le. a ) then
         ierr = 4
         x3 = longcote(1)
         y3 = longcote(2)
         x4 = 0d0
         y4 = a
         return
      endif
c
c     recherche de 2 points des cercles centres
c     en s2 de rayon l2, en s1 de rayon l4
      c32    = longcote(3) ** 2
      a      = atan(1d0) * 2d0
      anginc = a / nbp
c
c     iterations pour affiner l'angle aux sommets 1 et 2
      do 20 iter = 1,5
c        boucle sur les angles
 10      a = a - anginc
         cosa = cos( a )
         sina = sin( a )
         x4 = longcote(4) * cosa
         y4 = longcote(4) * sina
         x3 = longcote(1) - longcote(2) * cosa
         y3 =               longcote(2) * sina
         c3 = (x4-x3)**2 + (y4-y3)**2
         if( c3 .gt. c32 ) goto 10
c
c        diminution de l'increment de l'angle
         a      = a + anginc
         anginc = anginc / nbp
 20   continue
c
      ierr = 0
      return
      end

      subroutine quadplanlongcote( LCote, StQuad, ierr )
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c but : construire les 4 sommets d'un quadrangle plan respectant 
c ----- les 4 longueurs donnees de ses cotes
c
c entree:
c---------
c LCote: longueur des 4 cotes
c
c sorties:
c --------
c StQuad : Les 2 coordonnees des 4 sommets (x1,y1,x2,y2,...,y4)
c ierr   : 0 si StQuad est construit, 1 sinon
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c auteur : perronnet alain Laboratoire J.-L. Lions UPMC Paris  mars 2006
c2345x7..............................................................012
      double precision  LCote(4), StQuad(8), LongCote(4)
      double precision  dmax, x,y, x3,y3, x4,y4, d, e1f1, e2f1, sqrt
c
      ierr = 0
c
c     le cote le plus long devient le cote 1 du quadrangle
      dmax = max( LCote(1), LCote(2), LCote(3), LCote(4) )
      do 10 imax=1,4
         if( dmax .eq. LCote(imax) ) goto 20
 10   continue
c
c     le cote le plus long est le imax-eme cote
 20   LongCote(1) = LCote(imax)
      do 30 i=1,3
         j = imax + i
         if( j .gt. 4 ) j=j-4
         LongCote(i+1) = LCote(j)
 30   continue
      print *,'LCote=',LCote
      print *,'LongCote=',LongCote
c
c     construction des coordonnees des sommets superieurs du quadrangle
c     dans le repere e1=Le PLUS LONG des 4 COTES, e2= Rotation +Pi/2 de e1
      call soqu4l( LongCote, x3, y3, x4, y4, ierr )
      print *,'Sortie soqu4l: ierr=',ierr,' x3=',x3,' y3=',y3,
     %        ' x4=',x4,' y4=',y4
c
c     le premier sommet de StQuad est (0,0) dans le repere (f1,f2)
      StQuad(1) = 0d0
      StQuad(2) = 0d0
c
c     l'abscisse dans le repere (f1,f2) du second sommet est LCote(1)
      dmax = 1d-3 * dmax
      if( LCote(1) .gt. dmax ) then
         x = LCote(1)
      else
         x = dmax
      endif
      StQuad(3) = x
c     l'ordonnee du second sommet est nulle dans le repere (f1,f2)
      StQuad(4) = 0d0
c
      goto( 100, 200, 300, 400 ), imax
c
c     le cote le plus long est le premier  -------------------
 100  StQuad(5) = x3
      StQuad(6) = y3
c
      StQuad(7) = x4
      StQuad(8) = y4
      return
c
c     le cote le plus long est le second  ---------------------
c     f1 = S4S1 / ||S4S1|| et S1=(0,0)
c     le repere e1,e2 devient le repere f1,f2
 200  d = sqrt( x4 * x4 + y4 * y4 )
c
c     la matrice de rotation (ei,fj)
      e1f1 = -x4 / d
      e2f1 = -y4 / d
c     f2 est la rotation de Pi/2 de f1
c     e1f2 = -e2f1
c     e2f2 =  e1f1
c
c     le 3-eme sommet dans le repere (f1,f2) est l'ancien second (LCote(2),0)
      x = LCote(2) - x4
      y = - y4
      StQuad(5) =  x * e1f1 + y * e2f1
      StQuad(6) = -x * e2f1 + y * e1f1
c
c     le 4-eme sommet dans le repere (f1,f2) est l'ancien troisieme (x3,y3)
      x = x3 - x4
      y = y3 - y4
      StQuad(7) =  x * e1f1 + y * e2f1
      StQuad(8) = -x * e2f1 + y * e1f1
      return
c
c     le cote le plus long est le troisieme ----------------------
c     f1 = S3S4 / ||S3S4|| et S1=(0,0)
c     le repere e1,e2 devient le repere f1,f2
 300  d = sqrt( (x4-x3)**2 + (y4-y3)**2 )
c
c     la matrice de rotation (ei,fj)
      e1f1 = (x4-x3) / d
      e2f1 = (y4-y3) / d
c     f2 est la rotation de Pi/2 de f1
c     e1f2 = -e2f1
c     e2f2 =  e1f1
c
c     le 3-eme sommet dans le repere (f1,f2) est l'ancien premier (0,0)
      x = - x3
      y = - y3
      StQuad(5) =  x * e1f1 + y * e2f1
      StQuad(6) = -x * e2f1 + y * e1f1
c
c     le 4-eme sommet dans le repere (f1,f2) est l'ancien second (LCote(3),0)
      x = LCote(3) - x3
      y = - y3
      StQuad(7) =  x * e1f1 + y * e2f1
      StQuad(8) = -x * e2f1 + y * e1f1
      return
c
c     le cote le plus long est le quatrieme ------------------------
c     f1 = S2S3 / ||S2S3||
c     le repere e1,e2 devient le repere f1,f2
 400  x = x3 - LCote(4)
      d = sqrt( x * x + y3 * y3 )
c
c     la matrice de rotation (ei,fj)
      e1f1 = x  / d
      e2f1 = y3 / d
c     f2 est la rotation de Pi/2 de f1
c     e1f2 = -e2f1
c     e2f2 =  e1f1
c
c     le 3-eme sommet dans le repere (f1,f2) est l'ancien quatrieme (x4,y4)
      x = x4 - LCote(4)
      StQuad(5) =  x * e1f1 + y4 * e2f1
      StQuad(6) = -x * e2f1 + y4 * e1f1
c
c     le 4-eme sommet dans le repere (f1,f2) est l'ancien premier (0,0)
      x = - LCote(4)
      StQuad(7) =  x * e1f1
      StQuad(8) = -x * e2f1
      return
      end

      program q1q1
      double precision  LCote(4), StQuad(8), x,y, x0,y0
      double precision  xyrect(1:2,1:2), xyr(2), xy(2)
c
c     rectangle des parametres
      xyrect(1,1) = 0d0
      xyrect(2,1) =-2d0
      xyrect(1,2) = 2d0
      xyrect(2,2) = 0d0
c
      do 10 i=1,4
         print *,'Longueur du cote',i,'=?'
         read *,LCote(i)
 10   continue
c
      call quadplanlongcote( LCote, StQuad, ierr )
c
      do 20 i=1,4
         x0 = StQuad(2*i-1)
         y0 = StQuad(2*i)
         if( i .eq. 4 ) then
            i1 = 1
         else
            i1 = i + 1
         endif
         x = StQuad(2*i1-1)
         y = StQuad(2*i1)

         print 10020,i,x0,i,y0,i,sqrt((x-x0)**2+(y-y0)**2)
10020    format('X',i1,'=',g16.7,' Y',i1,'=',g16.7,
     %        ' Long Cote',i1,'=',g16.7)
 20   continue
c
      xyr(1) = -3.141493957416214d-05
      xyr(2) = 1.9999371679001987d0
      call q1quad( xyrect, StQuad, xyr,
     %             xy, ierr )
      print *
      print *,'Q1  (',xyr(1),',',xyr(2),' )=(',xy(1),',',xy(2),' )'
      call q1qinv( xyrect, stquad, xy,
     %             xyr, ierr )
      print *,'Q1-1(',xy(1),',',xy(2),' )=(',xyr(1),',',xyr(2),' )'
c
      xyr(1) = 1d0
      xyr(2) = 0.5d0
      call q1quad( xyrect, StQuad, xyr,
     %             xy, ierr )
      print *
      print *,'Q1  (',xyr(1),',',xyr(2),' )=(',xy(1),',',xy(2),' )'
      call q1qinv( xyrect, stquad, xy,
     %             xyr, ierr )
      print *,'Q1-1(',xy(1),',',xy(2),' )=(',xyr(1),',',xyr(2),' )'
c
      xyr(1) = (xyrect(1,1)+xyrect(1,2))/2
      xyr(2) = (xyrect(2,1)+xyrect(2,2))/2
      call q1quad( xyrect, StQuad, xyr,
     %             xy, ierr )
      print *
      print *,'Q1  (',xyr(1),',',xyr(2),' )=(',xy(1),',',xy(2),' )'
      call q1qinv( xyrect, stquad, xy,
     %             xyr, ierr )
      print *,'Q1-1(',xy(1),',',xy(2),' )=(',xyr(1),',',xyr(2),' )'
c
      xyr(1) = xyrect(1,2)
      xyr(2) = xyrect(2,2)
      call q1quad( xyrect, StQuad, xyr,
     %             xy, ierr )
      print *
      print *,'Q1  (',xyr(1),',',xyr(2),' )=(',xy(1),',',xy(2),' )'
      call q1qinv( xyrect, stquad, xy,
     %             xyr, ierr )
      print *,'Q1-1(',xy(1),',',xy(2),' )=(',xyr(1),',',xyr(2),' )'
c
      stop
      end


c -bash-> q1
c Longueur du cote 1=?0
c Longueur du cote 2=?3
c Longueur du cote 3=?7
c Longueur du cote 4=?5
c LCote=  0.  3.  7.  5.
c LongCote=  7.  5.  0.  3.
c Sortie soqu4l: ierr= 3 x3=2.35714286 y3=1.85576872 x4=2.35014286 y4=1.84876872
c X1=    0.000000     Y1=    0.000000     Long Cote1=   0.7000000E-02
c X2=   0.7000000E-02 Y2=    0.000000     Long Cote2=    2.993049    
c X3=    2.978978     Y3=  -0.3545251     Long Cote3=    7.000000    
c X4=   -1.970769     Y4=    4.595222     Long Cote4=    5.000000    
c
c Q1  (  0.,  0. )=(  0.,  0. )
c Q1-1(  0.,  0. )=( -2.15312999E-16,  2.87856587E-19 )
c
c Q1  (  1.,  0.5 )=( -0.182020564,  0.839446388 )
c Q1-1( -0.182020564,  0.839446388 )=(  1.,  0.5 )
c
c Q1  (  2.,  1. )=(  0.253802306,  1.06017434 )
c Q1-1(  0.253802306,  1.06017434 )=(  2.,  1. )
c
c Q1  (  4.,  2. )=(  2.97897835, -0.354525051 )
c Q1-1(  2.97897835, -0.354525051 )=(  4.,  2. )
c -bash-> 
