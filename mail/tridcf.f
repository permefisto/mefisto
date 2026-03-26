      SUBROUTINE TRIDCF( NBCF0,  nbstpe, nostpe, PXYD,   NOARST,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, N1ARTR, NOARTR,
     %                   MXARCF, N1ARCF, NOARCF, LARMIN,
     %                   NBTRCF, NOTRCF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRIANGULATION DIRECTE DE NBCF0 CONTOURS FERMES (CF)
C -----    DEFINIS PAR LA LISTE CIRCULAIRE DE LEURS ARETES PERIPHERIQUES
C          AVEC INTEGRATION DE NBSTPE SOMMETS ISOLES A L'UN DES CF INITIAUX
C
C ENTREES:
C --------
C NBCF0  : NOMBRE INITIAL DE CF A TRIANGULER
c nbstpe : nombre de sommets isoles a l'interieur des cf et
c          a devenir sommets de la triangulation
c nostpe : numero dans pxyd des nbstpe sommets isoles
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE
C MXSOAR : NOMBRE MAXIMAL D'ARETES STOCKABLES DANS LE TABLEAU NOSOAR
C          ATTENTION: MXSOAR>3*MXSOMM OBLIGATOIRE!
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C MXARCF : NOMBRE MAXIMAL D'ARETES DECLARABLES DANS NOARCF,N1ARCF,LARMIN,NOTRCF
C
C MODIFIES:
C ---------
C NOARST : NOARST(I) NUMERO D'UNE ARETE DE SOMMET I
C N1SOAR : NO DE L'EVENTUELLE PREMIERE ARETE LIBRE DANS LE TABLEAU NOSOAR
C          CHAINAGE DES VIDES SUIVANT EN 3 ET PRECEDANT EN 2 DE NOSOAR
C NOSOAR : NUMERO DES 2 SOMMETS, NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          HACHAGE DES ARETES = NOSOAR(1)+NOSOAR(2)*2
C          AVEC MXSOAR>=3*MXSOMM
C          UNE ARETE I DE NOSOAR EST VIDE <=> NOSOAR(1,I)=0 ET
C          NOSOAR(2,ARETE VIDE)=L'ARETE VIDE QUI PRECEDE
C          NOSOAR(3,ARETE VIDE)=L'ARETE VIDE QUI SUIT
C
C N1ARTR : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOARTR
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOARTR(2,.)
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1 = 0 SI TRIANGLE VIDE => ARETE2 = TRIANGLE VIDE SUIVANT
C
C N1ARCF : NUMERO DE LA PREMIERE ARETE DE CHACUN DES NBCF0 CF
C          N1ARCF(0)   NO DE LA PREMIERE ARETE VIDE DU TABLEAU NOARCF
C          NOARCF(2,I) NO DE L'ARETE SUIVANTE
C NOARCF : (1,*) NUMERO DU SOMMET
C          (2,*) NUMERO DE L'ARETE SUIVANTE DU CF DANS LE TABLEAU NOARCF
C          (3,*) NUMERO DE L'ARETE DANS LE TABLEAU NOSOAR
C
C AUXILIAIRES :
C -------------
C LARMIN : TABLEAU (MXARCF)   AUXILIAIRE
C          STOCKER LA LISTE DES NUMEROS DES MEILLEURES ARETES
C          LORS DE LA SELECTION DU MEILLEUR SOMMET DU CF A TRIANGULER
C          CF LE SP TRCHTD
C
C SORTIE :
C --------
C NBTRCF : NOMBRE DE  TRIANGLES DES NBCF0 CF
C NOTRCF : NUMERO DES TRIANGLES DES NBCF0 CF DANS LE TABLEAU NOARTR
C IERR   : 0 SI PAS D'ERREUR
C          2 SATURATION DE L'UN DES DES TABLEAUX NOSOAR, NOARTR, ...
C          3 SI CONTOUR FERME REDUIT A MOINS DE 3 ARETES
C          4 SATURATION DU TABLEAU NOTRCF
C          5 UNE ARETE APPARTENANT A 3 TRIANGLES DIFFERENTS
C         51 SATURATION DU TABLEAU NOSOAR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    MARS    1997
C MODIFS : ALAIN PERRONNET LABORATOIRE JL LIONS UPMC PARIS  OCTOBRE 2006
C....................................................................012
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           nostpe(nbstpe),
     %                  NOARTR(MOARTR,*),
     %                  NOSOAR(MOSOAR,MXSOAR),
     %                  NOARST(*),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(3,MXARCF),
     %                  LARMIN(MXARCF),
     %                  NOTRCF(MXARCF)
C
      integer           nosotr(3)
      double precision  d, diptdr, surtd2, dmin
C
C     DEPART AVEC NBCF0 CF A TRIANGULER
      NBCF = NBCF0
C
C     LE NOMBRE DE TRIANGLES FORMES DANS L'ENSEMBLE DES CF
      NBTRCF = 0
c
c     le nombre restant de sommets isoles a integrer au cf
      nbstp = nbstpe
c
 1    if( nbstp .le. 0 ) goto 10
c
c     il existe au moins un sommet isole
c     recherche d'un cf dont la premiere arete forme un triangle
c     d'aire>0 avec un sommet isole et recherche du sommet isole
c     le plus proche de cette arete
c     ==========================================================
      imin = 0
      dmin = 1d123
      do 6 ncf=1,nbcf
c        le cf en haut de pile a pour arete avant la premiere arete
         na1 = n1arcf( ncf )
         na2 = na1
c        recherche de l'arete qui precede la premiere arete
 2       if( noarcf( 2, na2 ) .ne. na1 ) then
            na2 = noarcf( 2, na2 )
            goto 2
         endif
c        l'arete na0 dans noarcf qui precede n1arcf( ncf )
         na0 = na2
c        la premiere arete du cf
         na1   = noarcf( 2, na0 )
c        son numero dans nosoar
         noar1 = noarcf( 3, na1 )
c        l'arete suivante
         na2   = noarcf( 2, na1 )
c        le no pxyd des 2 sommets de l'arete na1
         ns1   = noarcf( 1, na1 )
         ns2   = noarcf( 1, na2  )
         do 3 i=1,nbstpe
c           le sommet isole ns3
            ns3 = nostpe( i )
            if( ns3 .le. 0 ) goto 3
c           aire du triangle arete na1 et sommet ns3
            d = surtd2( pxyd(1,ns1), pxyd(1,ns2), pxyd(1,ns3) )
            if( d .gt. 0d0 ) then
c              distance de ce sommet ns3 a l'arete na1
               d = diptdr( pxyd(1,ns3),  pxyd(1,ns1), pxyd(1,ns2) )
               if( d .lt. dmin ) then
                  dmin = d
                  imin = i
               endif
            endif
 3       continue
         if( imin .gt. 0 ) then
c           le sommet imin de nostpe est a distance minimale de
c           la premiere arete du cf de numero ncf
c           la formation de l'arete ns2-ns3 dans le tableau nosoar
            call fasoar( ns2, ns3, -1, -1,  0,
     %                   mosoar, mxsoar, n1soar, nosoar, noarst,
     %                   noar2,  ierr )
            if( ierr .ne. 0 ) goto 9900
c           la formation de l'arete ns3-ns1 dans le tableau nosoar
            call fasoar( ns3, ns1, -1, -1,  0,
     %                   mosoar, mxsoar, n1soar, nosoar, noarst,
     %                   noar3,  ierr )
            if( ierr .ne. 0 ) goto 9900
c
c           ajout dans noartr du triangle de sommets ns1 ns2 ns3
c           et d'aretes na1, noar2, noar3 dans nosoar
            call trcf3a( ns1,   ns2,   ns3,
     %                   noar1, noar2, noar3,
     %                   mosoar, nosoar,
     %                   moartr, n1artr, noartr,
     %                   nt )
            if( nt .le. 0 ) then
               ierr = 7
               return
            endif
            if( nbtrcf .ge. mxarcf ) then
               write(imprim,*) 'tridcf: saturation du tableau notrcf'
               ierr = 8
               return
            endif
            nbtrcf = nbtrcf + 1
            notrcf( nbtrcf ) = nt
c
c           modification du cf. creation d'une arete dans noarcf
            na12 = n1arcf(0)
            if( na12 .le. 0 ) then
               write(imprim,*) 'tridcf: saturation du tableau noarcf'
               ierr = 10
               return
            endif
c           la 1-ere arete vide de noarcf est mise a jour
            n1arcf(0) = noarcf( 2, na12 )
c
c           l'arete suivante de na0
            noarcf( 1, na1 ) = ns1
            noarcf( 2, na1 ) = na12
            noarcf( 3, na1 ) = noar3
c           l'arete suivante de na1
            noarcf( 1, na12 ) = ns3
            noarcf( 2, na12 ) = na2
            noarcf( 3, na12 ) = noar2
c
c           un sommet isole traite
            nbstp = nbstp - 1
            nostpe( imin ) = - nostpe( imin )
            goto 1
         endif
c
 6    continue
c
      if( imin .eq. 0 ) then
         write(imprim,*) 'tridcf: Il reste',nbstp,
     %                   ' sommets isoles non triangules'
         write(imprim,*) 'Ameliorer l''algorithme'
         CALL XVPAUSE
         ierr = 9
         return
      endif
C
C     TANT QU'IL EXISTE UN CF A TRIANGULER FAIRE
C     LA TRIANGULATION DIRECTE DU CF
C     ==========================================
 10   IF( NBCF .GT. 0 ) THEN
C
C        LE CF EN HAUT DE PILE A POUR PREMIERE ARETE
         NA01 = N1ARCF( NBCF )
         NA1  = NOARCF( 2, NA01 )
C
C        CHOIX DU SOMMET DU CF A RELIER A L'ARETE NA1
C        --------------------------------------------
         CALL TRCHTD( PXYD, NA01, NA1, NOARCF,
     %                NA03, NA3,  LARMIN )
         IF( NA3 .EQ. 0 ) THEN
            IERR = 3
            RETURN
         ENDIF
C
C        L'ARETE SUIVANTE DE NA1
         NA02 = NA1
         NA2  = NOARCF( 2, NA1 )
C
C        FORMATION DU TRIANGLE ARETE NA1 - SOMMET NOARCF(1,NA3)
C        ------------------------------------------------------
         CALL TRCF3S( NBCF,   NA01, NA1, NA02, NA2, NA03, NA3,
     %                MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                MOARTR, N1ARTR, NOARTR, NOARST,
     %                MXARCF, N1ARCF, NOARCF, NT )
         IF( NT .LE. 0 ) THEN
C           ERREUR ou SATURATION DU TABLEAU NOARTR OU NOARCF OU N1ARCF
            IERR = 2
            RETURN
         ENDIF
C
C        AJOUT DU TRIANGLE CREE A SA PILE
         IF( NBTRCF .GE. MXARCF ) THEN
            WRITE(IMPRIM,*) 'TRIDCF: SATURATION DU TABLEAU NOTRCF'
            IERR = 4
            RETURN
         ENDIF
         NBTRCF = NBTRCF + 1
         NOTRCF( NBTRCF ) = NT
         GOTO 10
      ENDIF
C
C     MISE A JOUR DU CHAINAGE DES TRIANGLES DES ARETES
C     ================================================
      DO 30 NTP0 = 1, NBTRCF
C
C        LE NUMERO DU TRIANGLE AJOUTE DANS LE TABLEAU NOARTR
         NT0 = NOTRCF( NTP0 )
C
CCCC        AIRE SIGNEE DU TRIANGLE NT0
CCCC        LE NUMERO DES 3 SOMMETS DU TRIANGLE NT
CCC         CALL NUSOTR( NT0, MOSOAR, NOSOAR, MOARTR, NOARTR,
CCC     %                NOSOTR )
CCC         D = SURTD2( PXYD(1,NOSOTR(1)), PXYD(1,NOSOTR(2)),
CCC     %               PXYD(1,NOSOTR(3)) )
CCC         IF( D .LE. 0 ) THEN
CCCC
CCCC           UN TRIANGLE D'AIRE NEGATIVE DE PLUS
CCC            WRITE(IMPRIM,*) 'TRIANGLE ',NT0,' ST:',NOSOTR,
CCC     %                      ' D AIRE ',D,'<=0'
CCC            CALL XVPAUSE
CCC         ENDIF
C
C        TRACE DU TRIANGLE NT0
         CALL MTTRTR( PXYD, NT0, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                NCTURQ, NCBLAN )
C
C        BOUCLE SUR LES 3 ARETES DU TRIANGLE NT0
         DO 20 I=1,3
C
C           LE NUMERO DE L'ARETE I DU TRIANGLE DANS LE TABLEAU NOSOAR
            NOAR = ABS( NOARTR(I,NT0) )
C
C           CE TRIANGLE EST IL DEJA CHAINE DANS CETTE ARETE?
            NT1 = NOSOAR(4,NOAR)
            NT2 = NOSOAR(5,NOAR)
            IF( NT1 .EQ. NT0 .OR. NT2 .EQ. NT0 ) GOTO 20
C
C           AJOUT DE CE TRIANGLE NT0 A L'ARETE NOAR
            IF( NT1 .LE. 0 ) THEN
C               LE TRIANGLE EST AJOUTE A L'ARETE
                NOSOAR( 4, NOAR ) = NT0
            ELSE IF( NT2 .LE. 0 ) THEN
C               LE TRIANGLE EST AJOUTE A L'ARETE
                NOSOAR( 5, NOAR ) = NT0
            ELSE
C              L'ARETE APPARTIENT A 2 TRIANGLES DIFFERENTS DE NT0
C              ANOMALIE. CHAINAGE DES TRIANGLES DES ARETES DEFECTUEUX
C              A CORRIGER
               WRITE(IMPRIM,*) 'TRIDCF: ERREUR 1 ARETE DANS 3 TRIANGLES'
               WRITE(IMPRIM,*) 'TRIDCF: ARETE NOSOAR(',NOAR,')=',
     %                          (NOSOAR(K,NOAR),K=1,MOSOAR)
               CALL NUSOTR( NT0, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR)
               WRITE(IMPRIM,*) 'TRIDCF: TRIANGLE NT0=',NT0,' ST:',
     %                          (NOSOTR(K),K=1,3)
               CALL NUSOTR( NT1, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR)
               WRITE(IMPRIM,*) 'TRIDCF: TRIANGLE NT1=',NT1,' ST:',
     %                          (NOSOTR(K),K=1,3)
               CALL NUSOTR( NT2, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR)
               WRITE(IMPRIM,*) 'TRIDCF: TRIANGLE NT2=',NT2,' ST:',
     %                          (NOSOTR(K),K=1,3)
ccc               CALL XVPAUSE
               IERR = 5
               RETURN
            ENDIF
C
 20      CONTINUE
c
 30   continue
      return
c
c     erreur tableau nosoar sature
 9900 write(imprim,*) 'saturation du tableau nosoar'
      ierr = 51
      return
      end
