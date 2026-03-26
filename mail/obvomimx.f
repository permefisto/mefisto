      SUBROUTINE OBVOMIMX( KNMVOL, LADEFI, HEXAVOLU, NBSOCO, NBFACO,
     %                     IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LES COORDONNEES DU POINT MINIMUM ET MAXIMUM
C -----     DES SOMMETS DE L'OBJET VOLUMIQUE PARTITIONNE
C           CALCUL DU NOMBRE DE SOMMETS ET FACES DU CONTOUR DU VOLUME
C           VERIFICATION DE LA FERMETURE DES SURFACES

C ENTREES:
C --------
C KNMVOL :  NOM DU VOLUME A TETRAEDRISER
C LADEFI :  TABLEAU DE DEFINITION DU VOLUME PARTITIONNE
C           CF '~td/d/a_volume__definition'

C SORTIES:
C --------
C HEXAVOLU: MIN ET MAX DES COORDONNEES DES SOMMETS DES VOLUMES
C           C-A-D DE L'HEXADRE ENGLOBANT LE VOLUME
C NBSOCO  : NOMBRE DE POINTS INTERNES EN ENTREE
C           NOMBRE DE POINTS SOMMETS DU VOLUME EN SORTIE
C NBFACO  : NOMBRE DE FACES INITIALES DES SURFACES CONTOUR DU VOLUME
C IERR    : =0 SI PAS D'ERREUR
C           >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1991
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*24      KNMVOL
      REAL              HEXAVOLU(6,2), COIN(6,2)
      INTEGER           LADEFI(0:*)
      CHARACTER*24      KNOM

C     TYPE DU VOLUME
C     ==============
      NUTYVO = LADEFI( WUTYVO )

C     RECHERCHE DES COORDONNEES HEXAVOLU DE L'HEXAEDRE ENGLOBANT
C     ========================================================
      R = RINFO( 'GRAND' )
      HEXAVOLU(1,1) =  R
      HEXAVOLU(2,1) =  R
      HEXAVOLU(3,1) =  R
      HEXAVOLU(4,1) =  0.0
      HEXAVOLU(5,1) =  0.0
      HEXAVOLU(6,1) =  0.0

      HEXAVOLU(1,2) = -R
      HEXAVOLU(2,2) = -R
      HEXAVOLU(3,2) = -R
      HEXAVOLU(4,2) =  0.0
      HEXAVOLU(5,2) =  0.0
      HEXAVOLU(6,2) =  0.0

      NBFACO = 0
      DO N=1,LADEFI( WBVOPA )

C        BOUCLE SUR LES VOLUMES DE LA PARTITION
C        --------------------------------------
C        LE NUMERO DU VOLUME N
         IF( NUTYVO .NE. 20 ) THEN
            NOVO = LADEFI( WUVOPA - 1 + N )
         ELSE
            NOVO = LADEFI( WUVO20 - 1 + N )
         ENDIF

C        LE TABLEAU LEXIQUE DE CE VOLUME
         CALL LXNLOU( NTVOLU, NOVO, NTLXVO, MNDFVO )
C        LE TABLEAU 'DEFINITION' DE CE VOLUME
         CALL LXTSOU( NTLXVO, 'DEFINITION', NTDFVO, MNDFVO )
C        LE TYPE DEFINITION DE CE VOLUME
         NOTYVO = MCN( MNDFVO + WUTYVO )
         IF( NOTYVO .NE. 8 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'ERREUR: NON VOLUME D''UN VOLUME PARTITIONNE'
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF

C        LE TABLEAU DES SURFACES FERMEES DU CONTOUR EST COMPLETE
         NBSF1V = MCN( MNDFVO + WBSF1V )
         MNSV   = MNDFVO + WUSF1V - 1

         DO 10 IS=1,NBSF1V

C           BOUCLE SUR LES SURFACES FERMEES DU CONTOUR
C           ------------------------------------------
C           LE NUMERO DE LA SURFACE FERMEE IS DU CONTOUR
            NOSU = MCN( MNSV + IS )
C           LE LEXIQUE DE LA SURFACE
            CALL LXNLOU( NTSURF, NOSU, NTLXSU, MNLXSU )
            IF( NTLXSU .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:10),'(I10)') NOSU
               KERR(1)='obvomimx: SURFACE INCONNUE'//KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 2
               GOTO 10
            ENDIF

C           VERIFICATION : LA SURFACE EST-ELLE FERMEE ?
C           C-A-D CHAQUE ARETE APPARTIENT ELLE A 2 TRIANGLES
            CALL OBJFER( 3, NOSU, 1, NOFERM )
            IF( NOFERM .LE. 0 ) THEN
               CALL NMOBNU( 'SURFACE', NOSU, KNOM )
               NBLGRC(NRERR) = 1
               KERR(1) ='obvomimx: SURFACE NON FERMEE ' // KNOM
               CALL LEREUR
               IERR = 3
               GOTO 10
            ENDIF

C           LE TABLEAU 'NSEF' DE CETTE SURFACE FERMEE
            CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )
C           TOUTES LES VERIFICATIONS SUR LES NSEF ONT ETE
C           FAITES DANS CALL OBJFER

C           LE TYPE STRUCTURE OU NON DE LA SURFACE FERMEE
            NUTYMA = MCN( MNFASU + WUTYMA )
C           LE NOMBRE DE FACES DE LA SURFACE
            IF( NUTYMA .EQ. 0 ) THEN
C              SURFACE NON STRUCTUREE
               NBFASU = MCN( MNFASU + WBEFOB )
            ELSE IF( NUTYMA .EQ. 3 ) THEN
C              TRIANGLE STRUCTURE
               NBFASU = MCN( MNFASU + WBARTR ) ** 2
            ELSE
               NBFASU = 0
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:10),'(I10)') NOSU
               KERR(1) = 'obvomimx: SURFACE' // KERR(MXLGER)(1:10)
     %                // ' NON TRIANGULEE'
               CALL LEREUR
               IERR = 4
            ENDIF
C           LE NOMBRE DE FACES DES VOLUMES A TRAITER
            NBFACO = NBFACO + NBFASU

C           LE NOMBRE DE SOMMETS SANS IDENTIFICATION
C           LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
            CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOSU, MNSOSU )
            IF( NTSOSU .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:10),'(I10)') NOSU
               KERR(1) = 'obvomimx: ERREUR SURFACE SANS XYZ ' //
     %                    KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 5
               GOTO 9999
            ENDIF
            NBSOCO = NBSOCO + MCN( MNSOSU + WNBSOM )

C           LE CADRE EXTREME DES COORDONNEES EST MIS A JOUR
            CALL CADEXT( MNSOSU, COIN )
            NBCOOR = MCN( MNSOSU + WBCOOR )
            DO K=1,NBCOOR
               HEXAVOLU(K,1) = MIN( HEXAVOLU(K,1), COIN(K,1) )
               HEXAVOLU(K,2) = MAX( HEXAVOLU(K,2), COIN(K,2) )
            ENDDO

 10      ENDDO
      ENDDO

C     LES POINTS SONT ILS COPLANAIRES DANS X=Cte ou Y=Cte ou Z=Cte ?
      DO I=1,3
C        L'AMPLITUDE DE L'OBJET DANS LA DIRECTION I
         H = HEXAVOLU(I,2) - HEXAVOLU(I,1)
         IF( H .LT. 1E-3 * (ABS(HEXAVOLU(I,1))+ABS(HEXAVOLU(I,2)))) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='obvomimx: TOUS LES POINTS SONT COPLANAIRES'
            ELSE
               KERR(1)='obvomimx: ALL POINTS ARE ON ONE PLANE'
            ENDIF
            CALL LEREUR
            IERR = 6
            GOTO 9999
         ENDIF
      ENDDO

C     AUGMENTATION DE 1% DE L'HEXAEDRE POUR ENGLOBER STRICTEMENT LE VOLUME
      DO I=1,3
C        L'AMPLITUDE DE L'OBJET DANS LA DIRECTION I EST AUGMENTE DE 1%
         H = ( HEXAVOLU(I,2) - HEXAVOLU(I,1) ) * 1.01 / 2
         R = ( HEXAVOLU(I,2) + HEXAVOLU(I,1) ) / 2
         HEXAVOLU(I,1) = R - H
         HEXAVOLU(I,2) = R + H
      ENDDO

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10000, KNMVOL, ((HEXAVOLU(I,J),I=1,3),J=1,2)
         PRINT*,'NOMBRE FACES DES SURFACES CONTOUR DU VOLUME NBFACO=',
     %           NBFACO
         PRINT*,'NOMBRE SOMMETS INITIAUX DU VOLUME NBSOCO=',NBSOCO
      ELSE
         PRINT 20000, KNMVOL, ((HEXAVOLU(I,J),I=1,3),J=1,2)
         PRINT*,'FACES NUMBER of VOLUME CONTOUR SURFACES NBFACO=',
     %           NBFACO
         PRINT*,'INITIAL VERTICES NUMBER of VOLUME NBSOCO=',NBSOCO
      ENDIF

10000 FORMAT(
     % ' COORDONNEES MIN-MAX des SOMMETS du VOLUME ',A/
     % ' XMIN=',G15.7,'  YMIN=',G15.7,'  ZMIN=',G15.7 /
     % ' XMAX=',G15.7,'  YMAX=',G15.7,'  ZMAX=',G15.7 )

20000 FORMAT(
     % ' MIN-MAX VERTICES COORDINATES of VOLUME ',A/
     % ' XMIN=',G15.7,'  YMIN=',G15.7,'  ZMIN=',G15.7 /
     % ' XMAX=',G15.7,'  YMAX=',G15.7,'  ZMAX=',G15.7 )

 9999 RETURN
      END
