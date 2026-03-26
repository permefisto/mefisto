      SUBROUTINE AD2T2T( I,      NT,
     %                   NOTRIA, NOTRSO, PXYD, NUISOP,
     %                   NT1,    NBCHGT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECHANGE EVENTUEL DES DIAGONALES DE 2 TRIANGLES ADJACENTS PAR
C -----    L'ARETE I DU TRIANGLE NT ET SI QUARANGLE CONVEXE
C          POUR RENDRE LA TRIANGULATION DE TYPE DELAUNAY
C
C ENTREES:
C --------
C I      : NUMERO DE 1 A 3 DE L'ARETE
C NT     : NUMERO NOTRIA DU TRIANGLE A CONSIDERER
C PXYD   : X Y DISTANCE SOUHAITEE DES SOMMETS
C NOTRIA : LISTE CHAINEE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                          ADJACENT PAR L'ARETE i
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C PXYD   : X Y DISTANCE SOUHAITEE DES SOMMETS
C NUISOP : NUMERO DE L'ISO DU POINT ET 0 SINON
C
C SORTIE :
C --------
C NT1    : NUMERO NOTRIA DU TRIANGLE OPPOSE PAR L'ARETE I DU TRIANGLE NT
C NBCHGT : NOMBRE DE CHANGEMENT DE DIAGONALE
C          0 PAS D'ECHANGE
C          1 SI LA DIAGONALE A ETE CHANGEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1995
C....................................................................012
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION
C
      INTEGER           NOTRIA(6,*),
     %                  NOTRSO(*),
     %                  NUISOP(*)
      DOUBLE PRECISION  PXYD(3,*),
     %                  CETRIA(1:3)
      DOUBLE PRECISION  S, XG, YG
C
      NBCHGT = 0
C
C     LE TRIANGLE OPPOSE A L'ARETE I
      NT1 = NOTRIA(I+3,NT)
      IF( NT1 .LE. 0 ) GOTO 9999
      IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 9999
C
C     L'ARETE I DE NT EST COMMUNE A 2 TRIANGLES
      NS1 = NOTRIA(I,NT)
C
C     LE 2-EME SOMMET DU TRIANGLE NT
      IF( I .LT. 3 ) THEN
         I1 = I + 1
      ELSE
         I1 = 1
      ENDIF
      NS2 = NOTRIA(I1,NT)
C
      IF( NUISOP(NS1) .GT. 0 .AND. NUISOP(NS2) .GT. 0 .AND.
     %    NUISOP(NS1) .EQ. NUISOP(NS2) ) GOTO 9999
C
C     LE 3-EME SOMMET DU TRIANGLE NT
      IF( I .GT. 1 ) THEN
         I2 = I - 1
      ELSE
         I2 = 3
      ENDIF
      NS3 = NOTRIA(I2,NT)
C
C     LE NUMERO D'ARETE DE NS1-NS2 DANS NT1
      IF( NOTRIA(1,NT1) .EQ. NS2 ) THEN
         NA = 1
      ELSE IF( NOTRIA(2,NT1) .EQ. NS2 ) THEN
         NA = 2
      ELSE IF( NOTRIA(3,NT1) .EQ. NS2 ) THEN
         NA = 3
      ELSE
         WRITE(IMPRIM,*) 'AD2T2T: SOMMET ',NS2,' NON DANS TR ',NT1
         GOTO 9999
      ENDIF
      IF( NA .LT. 3 ) THEN
         NA1 = NA + 1
      ELSE
         NA1 = 1
      ENDIF
      IF( NA1 .LT. 3 ) THEN
         NA2 = NA1 + 1
      ELSE
         NA2 = 1
      ENDIF
      NS4 = NOTRIA(NA2,NT1)
C
C     LE SOMMET NS4 EST IL DANS LE CERCLE CIRCONSCRIT DE NT ?
C     -------------------------------------------------------
C     CETRIA : COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT ET
C              CARRE DU RAYON
      NA = -1
      CALL CENCED( PXYD(1,NS2), PXYD(1,NS3), PXYD(1,NS4),
     %             CETRIA, NA )
      IF( NA .NE. 0 ) GOTO 9999
C
      XG = CETRIA(1) - PXYD(1,NS4)
      YG = CETRIA(2) - PXYD(2,NS4)
      IF( XG*XG + YG*YG .LE. CETRIA(3)*0.9999D0 ) THEN
C
C        OUI: NS4 EST DANS LE CERCLE CIRCONSCRIT A NT
C        => NS3 EST AUSSI DANS LE CERCLE CIRCONSCRIT DE NT1
C        LA CONDITION AU DESSUS SUFFIT CAR
C        LA SOMME DES 2 ANGLES OPPOSES EST < 180 DEGRES
C
C        LE QUADRANGLE FORME PAR LES 2 TRIANGLES EST IL CONVEXE?
C        -------------------------------------------------------
C        2 FOIS LA SURFACE DU TRIANGLE 234
         S =(PXYD(1,NS3)-PXYD(1,NS2))*(PXYD(2,NS4)-PXYD(2,NS2))
     %     -(PXYD(2,NS3)-PXYD(2,NS2))*(PXYD(1,NS4)-PXYD(1,NS2))
         IF( S .LE. 0 ) GOTO 9999
C
C        2 FOIS LA SURFACE DU TRIANGLE 143
         S =(PXYD(1,NS4)-PXYD(1,NS1))*(PXYD(2,NS3)-PXYD(2,NS1))
     %     -(PXYD(2,NS4)-PXYD(2,NS1))*(PXYD(1,NS3)-PXYD(1,NS1))
         IF( S .LE. 0 ) GOTO 9999
C
C        ICI : LE DECOUPAGE ACTUEL DU QUADRANGLE CONVEXE N'EST PAS DELAUNAY
         GOTO 20
      ENDIF
C
C     ICI LE DECOUPAGE NT NT1 EST DEJA DE TYPE DELAUNAY
      GOTO 9999
C
C     ICI NT NT1 DOIVENT ETRE REMPLACES SUR EUX MEMES POUR ETRE DELAUNAY
C     ==================================================================
 20   NBCHGT = 1
C
      IF( TRATRI ) THEN
C        TRACE DES TRIANGLES A MODIFIER
         CALL  DVTRTR( PXYD , NOTRIA , NT  , NCJAUN , NCNOIR )
         CALL  DVTRTR( PXYD , NOTRIA , NT1 , NCJAUN , NCNOIR )
      ENDIF
C
C     MISE A JOUR DES 2 SOMMETS
      NOTRIA( I,NT ) = NS4
      NOTRIA(NA,NT1) = NS3
C
C     LES TRIANGLES OPPOSES
      NTT1 = NOTRIA(I1+3,NT)
      NTT2 = NOTRIA(I2+3,NT)
      NTT3 = NOTRIA(NA1+3,NT1)
      NTT4 = NOTRIA(NA2+3,NT1)
C
      NOTRIA(I +3 ,NT ) = NTT4
      NOTRIA(I1+3 ,NT ) = NTT1
      NOTRIA(I2+3 ,NT ) = NT1
C
      NOTRIA(NA +3,NT1) = NTT2
      NOTRIA(NA1+3,NT1) = NTT3
      NOTRIA(NA2+3,NT1) = NT
C
C     L'INVERSE POUR LES TRIANGLES NTT2 ET NTT4
      IF( NTT4 .GT. 0 ) THEN
C        LE NUMERO D'ARETE DE NS1-NS2 DANS NT1
         IF( NOTRIA(1,NTT4) .EQ. NS2 ) THEN
            NA = 4
         ELSE IF( NOTRIA(2,NTT4) .EQ. NS2 ) THEN
            NA = 5
         ELSE
            NA = 6
         ENDIF
         NOTRIA(NA,NTT4) = NT
      ENDIF
C
      IF( NTT2 .GT. 0 ) THEN
         IF( NOTRIA(1,NTT2) .EQ. NS1 ) THEN
            NA = 4
         ELSE IF( NOTRIA(2,NTT2) .EQ. NS1 ) THEN
            NA = 5
         ELSE
            NA = 6
         ENDIF
         NOTRIA(NA,NTT2) = NT1
      ENDIF
C
C     LE TRIANGLE CONTENANT LES 4 SOMMETS
      NOTRSO(NS1) = NT1
      NOTRSO(NS2) = NT
      NOTRSO(NS3) = NT1
      NOTRSO(NS4) = NT
C
      IF( TRATRI ) THEN
C        TRACE DES 2 TRIANGLES APRES ECHANGE DE LA DIAGONALE
         CALL  DVTRTR( PXYD , NOTRIA , NT  , NCVERT, NCBLAN )
         CALL  DVTRTR( PXYD , NOTRIA , NT1 , NCVERT, NCBLAN )
      ENDIF
C
 9999 RETURN
      END
