      SUBROUTINE SE30TGA( QUALMN, COSMAXPL, MXSOM,  NBSOM,  XYZSOM,
     %                      NT0,    M1TRIA,   NBTRIA, NOTRIA, MODIF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     POUR AMELIORER LA QUALITE D'UNE TRIANGULATION
C -----     TRAITER LE PLUS GRAND ANGLE >100 degres DU TRIANGLE NT0 PAR
C           CREATION DU SOMMET NS5 MILIEU DU SEGMENT DES 2 SOMMETS PROJETES
C           SUR L'ARETE OPPOSEE et
C           DECOUPAGE EN 2 DE NT0 ET DU TRIANGLE NT1 OPPOSE A L'ARETE COMMUNE
C           2 TRIANGLES => 4 TRIANGLES DONC 2 AJOUTES
C           et SEULEMENT SI LA QUALITE MINIMALE DES TRIANGLES EST AMELIOREE

C ENTREES:
C --------
C QUALMN : QUALITE MINIMALE AU DESSOUS DE LAQUELLE LA QUALITE DOIT
C          ETRE AMELIOREE
C COSMAXPL: COSINUS DE L'ANGLE DIEDRE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C           LES 2 TRIANGLES SONT CONSIDERES NON-COPLANAIRES
C NT0    :  NUMERO NOTRIA DU TRIANGLE D'ANGLE TROP GRAND
C NBSOM  :  NOMBRE DE SOMMETS DU MAILLAGE
C M1TRIA :  NOMBRE DE MOTS DE NOTRIA POUR UN TRIANGLE ( 3 ou 4 ou 6 )
C NBTRIA :  NOMBRE DE TRIANGLES DU MAILLAGE

C MODIFIES:
C ---------
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS
C NOTRIA : NUMERO DES 3 SOMMETS ET 3 TRIANGLES ADJACENTS PAR LES ARETES

C SORTIES:
C --------
C NBSOM  : NOMBRE DE SOMMETS   DU MAILLAGE
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE
C MODIF  : 1 SI LE PLUS GRAND ANGLE DU TRIANGLE NT0 A ETE DECOUPE EN 2
C          0 SI LE TRIANGLE NT0 N'A PAS ETE MODIFIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET  Veulettes & St Pierre du Perray    Novembre 2019
C2345X7..............................................................012
      PARAMETER        (MXTRIT=128)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0
      CHARACTER*48      KNM

      INTEGER           NOTRIA(M1TRIA,*), LITRIT(MXTRIT), NOSOTR(3)
      REAL              XYZSOM(3,MXSOM), XYZ(3,3), LONARE(3), VECT(3,2),
     %                  XYZM23(3), XYZPR4(3)

      TRACTE0 = TRACTE
      MODIF   = 0
      NT12    = 0
      NT31    = 0
      NT24    = 0
      NT43    = 0

C     CONVERSION RADIANS DEGRES
      RADEGR = 45. / ATAN( 1. )

C     CALCUL DE LA QUALITE Q0 DU TRIANGLE NT0
      CALL QUATRI( NOTRIA(1,NT0), XYZSOM, Q0 )
      IF( Q0 .GE. QUALMN ) GOTO 9999

C     TRACE DU TRIANGLE NT0 ET DES 3 TRIANGLES ADJACENTS
      LITRIT(1) = NT0
      KNM = 'Debut se30tga: Triangle                     '
      WRITE(KNM(27:33),'(I7)') NT0
      CALL SANSDBL( KNM, NC )
      CALL TRTRIAN( KNM(1:NC), XYZSOM, M1TRIA, MXTRIT,
     %              1, LITRIT, NOTRIA )

C     LE TRIANGLE NT0 EST DE MAUVAISE QUALITE
C     CALCUL DES ANGLES MIN MAX ET PERIMETRE DU TRIANGLE NT0
C     ------------------------------------------------------
      DO J=1,3
C        LES 3 COORDONNEES DU SOMMET J DU TRIANGLE NT0
         NS = NOTRIA( J, NT0 )
         XYZ( 1, J ) = XYZSOM( 1, NS )
         XYZ( 2, J ) = XYZSOM( 2, NS )
         XYZ( 3, J ) = XYZSOM( 3, NS )
      ENDDO

C     LA LONGUEUR DES 3 COTES DU TRIANGLE NT0
      DO J=1,3
         IF( J .EQ. 3 ) THEN
            J1=1
         ELSE
            J1 = J+1
         ENDIF
         LONARE(J) = SQRT( ( XYZ(1,J)-XYZ(1,J1) ) ** 2
     %                   + ( XYZ(2,J)-XYZ(2,J1) ) ** 2
     %                   + ( XYZ(3,J)-XYZ(3,J1) ) ** 2 )
      ENDDO

C     CALCUL DU COSINUS MIN ET MAX DES 3 ANGLES DU TRIANGLE NT0
      COSMA = -2.
      COSMI =  2.
      DO J=1,3

         IF( J .EQ. 1 ) THEN
            J0 = 3
         ELSE
            J0 = J - 1
         ENDIF

         IF( J .EQ. 3 ) THEN
            J1 = 1
         ELSE
            J1 = J + 1
         ENDIF

C        LE COSINUS DE L'ANGLE J
         DO KL=1,3
            VECT(KL,1) = XYZ(KL,J0) - XYZ(KL,J)
            VECT(KL,2) = XYZ(KL,J1) - XYZ(KL,J)
         ENDDO
         COSJ = PROSCR( VECT(1,1), VECT(1,2), 3 )

         IF( LONARE(J0) .NE. 0. .AND. LONARE(J) .NE. 0. ) THEN
            COSJ = COSJ / ( LONARE(J0) * LONARE(J) )
         ELSE
            COSJ = 1.
         ENDIF

         IF( COSJ .LT. COSMI ) THEN
            COSMI = COSJ
C           NUMERO DE L'ANGLE MAX
            IAMAX = J
         ENDIF
         IF( COSJ .GT. COSMA ) THEN
            COSMA = COSJ
C           NUMERO DE L'ANGLE MIN
            IAMIN = J
         ENDIF

      ENDDO

C     TRAITEMENT DES ERREURS D'ARRONDIS DU COSINUS
      IF( COSMA .GT. 1. ) THEN
         COSMA = 1.
      ENDIF
      IF( COSMI .LT. -1. ) THEN
         COSMI = -1.
      ENDIF

C     LES ANGLES Min et MAX en DEGRES du TRIANGLE NT0
      ANGLMAX = ACOS(COSMI) * RADEGR
      ANGLMIN = ACOS(COSMA) * RADEGR

      IF( ANGLMAX .LE. 100.0 ) GOTO 9999


C     UN ANGLE EN DEGRES MAXIMUM DU TRIANGLE NT0 > 100 degres ->
C     NT0 EST DECOUPE EN 2 TRIANGLES A PARTIR DU MILIEU DU SEGMENT
C     DES 2 SOMMETS PROJETES SUR LEUR ARETE OPPOSEE COMMUNE et
C     SEULEMENT SI LA QUALITE MINIMALE DES TRIANGLES EST AUGMENTEE
C     ------------------------------------------------------------
      NSENS = 0

C     PAR PERMUTATION => 1 SEUL CAS IAMAX=1 NUMERO DU SOMMET DE
C                                           TROP GRAND ANGLE
      IF( IAMAX .EQ. 2 ) THEN
         I                = NOTRIA( 1, NT0 )
         NOTRIA( 1, NT0 ) = NOTRIA( 2, NT0 )
         NOTRIA( 2, NT0 ) = NOTRIA( 3, NT0 )
         NOTRIA( 3, NT0 ) = I
         I                = NOTRIA( 4, NT0 )
         NOTRIA( 4, NT0 ) = NOTRIA( 5, NT0 )
         NOTRIA( 5, NT0 ) = NOTRIA( 6, NT0 )
         NOTRIA( 6, NT0 ) = I
      ELSE IF( IAMAX .EQ. 3 ) THEN
         I                = NOTRIA( 1, NT0 )
         NOTRIA( 1, NT0 ) = NOTRIA( 3, NT0 )
         NOTRIA( 3, NT0 ) = NOTRIA( 2, NT0 )
         NOTRIA( 2, NT0 ) = I
         I                = NOTRIA( 4, NT0 )
         NOTRIA( 4, NT0 ) = NOTRIA( 6, NT0 )
         NOTRIA( 6, NT0 ) = NOTRIA( 5, NT0 )
         NOTRIA( 5, NT0 ) = I
      ENDIF

C     CALCUL DE LA QUALITE Q0 DU TRIANGLE NT0
      CALL QUATRI( NOTRIA(1,NT0), XYZSOM, Q0 )
      IF( Q0 .GE. QUALMN ) GOTO 9995

C     AFFICHAGE DU TRIANGLE NT0
      PRINT 10015, NT0,Q0,IAMAX,ANGLMAX,ANGLMIN,LONARE,
     %            (NOTRIA(L,NT0),(XYZ(K,L),K=1,3),L=1,3)

10015 FORMAT(/' se30tga: traitement du TRIANGLE',I7,': QUALITE=',F8.5,
     %   ' avec l''ANGLE ',I1,' MAX de ',F7.2,
     %   ' DEGRES et un ANGLE MIN de ',F7.2,' DEGRES'/
     %   ' LONGUEUR des 3 ARETES=',3G15.6/
     %    3(' se30tga: St',I8,' : X=',G15.6,' Y=',G15.6,' Z=',G15.6/))

C     LES 3 SOMMETS DE NT0
      NS1 = NOTRIA( 1, NT0 )
      NS2 = NOTRIA( 2, NT0 )
      NS3 = NOTRIA( 3, NT0 )

C     LES 3 COORDONNEES DU POINT PROJETE NS1 SUR LA DROITE NS2-NS3
      CALL PTPRDR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %             XYZM23 )

C     LE TRIANGLE NT1 OPPOSE PAR L'ARETE 2 OPPOSE AU SOMMET 1
C     ( DE PLUS GRAND ANGLE ) DE NT0 SOIT ENCORE
C     LE TRIANGLE ADJACENT A L'ARETE NOTRIA(2:3,NT0)
      NT1 = NOTRIA( 5, NT0 )

      IF( NT1 .GT. 0 ) THEN

C        CALCUL DE LA QUALITE Q1 DU TRIANGLE NT1
         CALL QUATRI( NOTRIA(1,NT1), XYZSOM, Q1 )
         Q0Q1MIN = MIN( Q0, Q1 )

C        NS4 3-EME SOMMET DE NT1 N'APPARTENANT PAS A L'ARETE NS2-NS3 DE NT0
         DO I=1,3
            NS4 = NOTRIA(I,NT1)
            IF( NS4 .NE. NS2 .AND. NS4 .NE. NS3 ) GOTO 15
         ENDDO

C        LES 3 COORDONNEES DU POINT PROJETE NS4 SUR LA DROITE NS2-NS3
 15      CALL PTPRDR( XYZSOM(1,NS4), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %                XYZPR4 )

C        LES 3 COORDONNEES DU MILIEU DE L'ARETE DES 2 POINTS
C        PROJETES NS1->PR1, NS4->PR4 SUR L'ARETE NS2-NS3 OPPOSEE
C        AU PLUS GRAND ANGLE DU TRIANGLE NT0 DE SOMMET NS1
         DO I=1,3
            XYZM23(I) = ( XYZM23(I) + XYZPR4(I) ) / 2
         ENDDO

      ELSE

C        PAS DE TRIANGLE OPPOSE NT1 -> CHOIX DU POINT PROJETE POUR
C        DIMINUER LE GRAND ANGLE ET AVOIR 90 DEGRES SEULEMENT AU
C        POINT PROJETE
         Q0Q1MIN = Q0
         NS4     = 0

      ENDIF


C     CREATION DU SOMMET NS5 MILIEU DE L'ARETE 2 DE NT0
C     ET DECOUPAGE EN 2 DE NT0 ET EVENTUELLEMENT DE NT1
C     2T -> 4T TRIANGLES DE SOMMET COMMUN NS5
C     -------------------------------------------------
      IF( NBSOM .GE. MXSOM ) THEN
         PRINT *,'se30tga: TABLEAU XYZSOM(',MXSOM,') SATURE'
         MODIF = 0
         GOTO 9999
      ENDIF

C     CREATION DU SOMMET NS5 MILIEU DU SEGMENT DES SOMMETS PROJETES
      NBSOM = NBSOM + 1
      NS5   = NBSOM
      DO I=1,3
         XYZSOM(I,NS5) = XYZM23(I)
      ENDDO

C     CALCUL DE LA QUALITE DES 2 SOUS-TRIANGLES DE SOMMET NS5
      NOSOTR(1) = NS1
      NOSOTR(2) = NS2
      NOSOTR(3) = NS5
      CALL QUATRI( NOSOTR, XYZSOM, Q2 )
      IF( Q2 .LE. Q0Q1MIN ) THEN
C        LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
         NBSOM = NBSOM - 1
         GOTO 9990
      ENDIF

      NOSOTR(1) = NS1
      NOSOTR(2) = NS5
      NOSOTR(3) = NS3
      CALL QUATRI( NOSOTR, XYZSOM, Q4 )
      IF( Q4 .LE. Q0Q1MIN ) THEN
C        LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
         NBSOM = NBSOM - 1
         GOTO 9990
      ENDIF

C     LE TRIANGLE NT1 OPPOSE PAR L'ARETE 2 OPPOSE AU SOMMET 1
C     ( DE PLUS GRAND ANGLE ) DE NT0 SOIT ENCORE
C     LE TRIANGLE ADJACENT A L'ARETE NOTRIA(2:3,NT0) EST
C     NT1 = NOTRIA( 5, NT0 )
C     LES 2 AUTRES TRIANGLES ADJACENTS AUX ARETES DE NT0
      NT12 = NOTRIA( 4, NT0 )
      NT31 = NOTRIA( 6, NT0 )

      IF( NT1 .LE. 0 ) THEN

C        PAS DE TRIANGLE NT1 OPPOSE A L'ARETE NS2-NS3
C        --------------------------------------------
         NS4 = 0
         NT1 = 0
         NT3 = 0

C        LE NOUVEAU SOUS-TRIANGLE NT2 DE NT0
         NBTRIA = NBTRIA + 1
         NT2    = NBTRIA
 
C        LES 4 TRIANGLES DE SOMMET CENTRAL NS5
         NOTRIA(1,NT0) = NS1
         NOTRIA(2,NT0) = NS2
         NOTRIA(3,NT0) = NS5
         NOTRIA(4,NT0) = NT12
         NOTRIA(5,NT0) = 0
         NOTRIA(6,NT0) = NT2

         NOTRIA(1,NT2) = NS1
         NOTRIA(2,NT2) = NS5
         NOTRIA(3,NT2) = NS3
         NOTRIA(4,NT2) = NT0
         NOTRIA(5,NT2) = 0
         NOTRIA(6,NT2) = NT31

         PRINT*,'se30tga: 1T->2T NT0=',NT0,
     %          ' NS5=',NS5,' NOUVEAU SOMMET MILIEU de l''ARETE',NS2,NS3

C        ESSAI D'ECHANGER LES 2 DIAGONALES POUR AMELIORER LA QUALITE
         CALL EC2DIA( COSMAXPL, NT0, 1, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT2, 3, NOTRIA, NBSOM, XYZSOM,  N )
         GOTO 1000

      ELSE

C        IL EXISTE UN TRIANGLE NT1 OPPOSE A L'ARETE NS2-NS3
C        --------------------------------------------------
C        RECHERCHE DE L'ARETE NS2-NS3 DANS NT1
         DO I=1,3
            IF( I .NE. 3 ) THEN
               I1 = I + 1
            ELSE
               I1 = 1
            ENDIF
            IF( NOTRIA(I ,NT1) .EQ. NS3 .AND.
     %          NOTRIA(I1,NT1) .EQ. NS2 ) THEN
               NSENS = 1
               GOTO 20
            ELSE IF( NOTRIA(I ,NT1) .EQ. NS2 .AND.
     %               NOTRIA(I1,NT1) .EQ. NS3 ) THEN
               NSENS = -1
               GOTO 20
            ENDIF
         ENDDO

C        NS4 3-EME SOMMET DE NT1 N'APPARTENANT PAS A L'ARETE NS2-NS3 DE NT0
 20      IF( I1 .NE. 3 ) THEN
            I2 = I1 + 1
         ELSE
            I2 = 1
         ENDIF
         NS4 = NOTRIA(I2,NT1)

         NOSOTR(1) = NS2
         NOSOTR(2) = NS4
         NOSOTR(3) = NS5
         CALL QUATRI( NOSOTR, XYZSOM, Q3 )
         IF( Q3 .LE. Q0Q1MIN ) THEN
C           LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
            NBSOM = NBSOM - 1
            GOTO 9990
         ENDIF

         NOSOTR(1) = NS5
         NOSOTR(2) = NS4
         NOSOTR(3) = NS3
         CALL QUATRI( NOSOTR, XYZSOM, Q5 )
         IF( Q5 .LE. Q0Q1MIN ) THEN
C           LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
            NBSOM = NBSOM - 1
            GOTO 9990
         ENDIF

         IF( NSENS .GT. 0 ) THEN
            NT24 = NOTRIA(3+I1,NT1)
            NT43 = NOTRIA(3+I2,NT1)
         ELSE
            NT24 = NOTRIA(3+I2,NT1)
            NT43 = NOTRIA(3+I1,NT1)
         ENDIF

C        LES 2 NOUVEAUX SOUS-TRIANGLES DE NT0
         NBTRIA = NBTRIA + 1
         NT2    = NBTRIA

         NBTRIA = NBTRIA + 1
         NT3    = NBTRIA
 
C        LES 4 TRIANGLES DE SOMMET CENTRAL NS5
         NOTRIA(1,NT0) = NS1
         NOTRIA(2,NT0) = NS2
         NOTRIA(3,NT0) = NS5
         NOTRIA(4,NT0) = NT12
         NOTRIA(5,NT0) = NT1
         NOTRIA(6,NT0) = NT2

         NOTRIA(1,NT1) = NS2
         NOTRIA(2,NT1) = NS4
         NOTRIA(3,NT1) = NS5
         NOTRIA(4,NT1) = NT24
         NOTRIA(5,NT1) = NT3
         NOTRIA(6,NT1) = NT0

         NOTRIA(1,NT2) = NS1
         NOTRIA(2,NT2) = NS5
         NOTRIA(3,NT2) = NS3
         NOTRIA(4,NT2) = NT0
         NOTRIA(5,NT2) = NT3
         NOTRIA(6,NT2) = NT31

         NOTRIA(1,NT3) = NS5
         NOTRIA(2,NT3) = NS4
         NOTRIA(3,NT3) = NS3
         NOTRIA(4,NT3) = NT1
         NOTRIA(5,NT3) = NT43
         NOTRIA(6,NT3) = NT2
      
C        MISE A JOUR DANS LES TRIANGLES ADJACENTS DES TRIANGLES MODIFIES
         IF( NT43 .GT. 0 ) THEN
            DO I=4,6
               IF( NOTRIA(I,NT43) .EQ. NT1 ) GOTO 110
            ENDDO
            PRINT*,'PB se30tga: TRIANGLE NT1=',NT1,
     %             ' NON ADJACENT AU TRIANGLE NT43=',NT43
            GOTO 9990
 110        NOTRIA(I,NT43) = NT3
         ENDIF

         IF( NT31 .GT. 0 ) THEN
            DO I=4,6
               IF( NOTRIA(I,NT31) .EQ. NT0 ) GOTO 115
            ENDDO
            PRINT*,'PB se30tga: TRIANGLE NT0=',NT0,
     %           ' NON ADJACENT AU TRIANGLE NT31=',NT31
            GOTO 9990
 115        NOTRIA(I,NT31) = NT2
         ENDIF

         PRINT*,'se30tga: 2T->4T NT0=',NT0,' NT1=',NT1,
     %          ' NS5=',NS5,' NOUVEAU SOMMET MILIEU de l''ARETE',NS2,NS3

C        ESSAI D'ECHANGER LES 2 DIAGONALES POUR AMELIORER LA QUALITE
         CALL EC2DIA( COSMAXPL, NT0, 1, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT1, 1, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT2, 3, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT3, 2, NOTRIA, NBSOM, XYZSOM,  N )
         GOTO 1000

      ENDIF


C     MODIFICATION EFFECTUEE DU TRIANGLE NT0
C     TRACE DES TRIANGLES MODIFIES ET DE LEURS 3 TRIANGLES ADJACENTS
C     ==============================================================
 1000 MODIF = 1

      PRINT*,'se30tga: le TRIANGLE',NT0,
     %       ' St:',(NOTRIA(K,NT0),K=1,3),
     %       ' Angle Min=',ANGLMIN,
     %       ' son PLUS GRAND ANGLE',ANGLMAX,' a ETE DECOUPE'

C     TRACE DES TRIANGLES CREES ET DE LEURS ADJACENTS
      LITRIT(1) = NT0
      LITRIT(2) = NT1
      LITRIT(3) = NT2
      LITRIT(4) = NT3
      N = 4
C     COMPRESSION DU TABLEAU LITRIT PAR SUPPRESSION DES ZEROS
      M = 0
      DO K=1,N
         L = LITRIT(K)
         IF( L .GT. 0 ) THEN
            M = M+1
            LITRIT(M) = L
         ENDIF
      ENDDO
C     TRACE DES M TRIANGLES CREES ET DE LEURS ADJACENTS
      CALL TRTRIAN( 'se30tga modif=1:', XYZSOM, M1TRIA, MXTRIT, M,
     %               LITRIT, NOTRIA )
      GOTO 9999


C     PAS DE MODIFICATION DU TRIANGLE NT0
C     ===================================
 9990 MODIF = 0

      PRINT*,'se30tga: TRIANGLE',NT0,
     %       ' St:',(NOTRIA(K,NT0),K=1,3),
     %       ' Qualite=',Q0,
     %       ' Angle Min=',ANGLMIN,
     %       ' Angle Max=',ANGLMAX,' NON MODIFIE'

C     TRACE DU TRIANGLE NT0 et de ses VOISINS
 9995 LITRIT(1) = NT0
      CALL TRTRIAN( 'se30tga: 0 modif', XYZSOM, M1TRIA, MXTRIT, 1,
     %               LITRIT, NOTRIA )

 9999 TRACTE = TRACTE0
      RETURN
      END
