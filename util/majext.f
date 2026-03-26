      SUBROUTINE MAJEXT( MNSOMM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   MISE A JOUR DU CADRE MINMAX DES NBCOOR COORDONNEES DES SOMMETS
C ----

C ENTREE  :
C ---------
C MNSOMM  : ADRESSE MCN DU TABLEAU 'XYZSOMMET' A TRAITER

C SORTIES :
C ---------
C DANS LE COMMON /XYZEXT/
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2005
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))

      COMMON / MSSFTA / MOTSSF(27), MOTSMN, NTADAM
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"
C.......................................................................

C     PROTECTION
      IF( MNSOMM .LE. 0 .OR. MNSOMM .GT. MOTSMN ) RETURN

C     LE NOMBRE DE XYZ
      NBSOM = MCN( MNSOMM + WNBSOM )

C     LE NOMBRE DE COORDONNEES D'UN SOMMET
      NBCOOR = MCN( MNSOMM + WBCOOR )

C     L'ADRESSE -1 DES COORDONNEES
      MN = MNSOMM + WYZSOM - 1

      IF( INIEXT .EQ. 0 .AND. NBSOM .GT. 1 ) THEN

C        COOEXT N'EST PAS INITIALISE => MIN MAX DE CE PLSV
         DO I=1,NBCOOR
            R = RMCN( MN + I )
            COOEXT(I,1) = R
            COOEXT(I,2) = R
         ENDDO
         MN = MN + NBCOOR

         DO J=2,NBSOM
            DO I=1,NBCOOR
C              LA COORDONNEE I DU POINT J
               R = RMCN( MN + I )
               IF( R .LT. COOEXT(I,1) ) THEN
                  COOEXT(I,1) = R
               ENDIF
               IF( R .GT. COOEXT(I,2) ) THEN
                  COOEXT(I,2) = R
               ENDIF
            ENDDO
            MN = MN + NBCOOR
         ENDDO
C        MODIFICATION DE COOEXT
         INIEXT = 2

      ELSE

C        LA BOUCLE SUR LES NBCOOR COORDONNEES DES NBSOM SOMMETS
C        REINITIALISATION DE COOEXT
         INIEXT = 1
         DO J=1,NBSOM
            DO I=1,NBCOOR
C              LA COORDONNEE I DU POINT J
               R = RMCN( MN + I )

C              MISE A JOUR DE LA COORDONNEE MIN
               IF( R .LT. COOEXT(I,1) ) THEN
                  COOEXT(I,1) = R
C                 MODIFICATION DE COOEXT
                  INIEXT = 2
               ENDIF
C              MISE A JOUR DE LA COORDONNEE MAX
               IF( R .GT. COOEXT(I,2) ) THEN
                  COOEXT(I,2) = R
C                 MODIFICATION DE COOEXT
                  INIEXT = 2
               ENDIF
            ENDDO
            MN = MN + NBCOOR
         ENDDO
      ENDIF

cccC     VALEURS PAR DEFAUT DES VARIABLES DE $MEFISTO/incl/xyzext.inc
ccc      DO I=1,NBCOOR
ccc         R = ( COOEXT(I,2) - COOEXT(I,1) ) / 20  =>  OVERFLOW!
ccc         COOEXT(I,1) = COOEXT(I,1) - R
ccc         COOEXT(I,2) = COOEXT(I,2) + R
ccc      ENDDO

      DO I=NBCOOR+1,6
         COOEXT(I,1) = 0
         COOEXT(I,2) = 0
      ENDDO

      MOAXYZ = 0
      XYZAMPLI(1) = 1
      XYZAMPLI(2) = 1
      XYZAMPLI(3) = 1
      XYZAMPLI(4) = 1

cccC     AFFICHAGE DES COORDONNEES MIN et MAX
ccc      IF( INIEXT .EQ. 2 ) THEN
ccc         WRITE(IMPRIM,*)
ccc         DO 90 I=1,NBCOOR
ccc            WRITE(IMPRIM,10090) I, COOEXT(I,1), COOEXT(I,2)
ccc10090       FORMAT('majext: COORD ',I1,' MIN=',G15.7,' MAX=',G15.7)
ccc 90      ENDDO
ccc      ENDIF

      RETURN
      END
