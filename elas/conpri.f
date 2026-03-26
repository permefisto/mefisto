      SUBROUTINE CONPRI( NUELEM, NBELEM, NDIM, NCOMP, NPI, NDSM, NOFORE,
     %                   COOPTC, STRESS,
     %                   CONTMI, NOPMIN, NOCMIN,
     %                   CONTMX, NOPMAX, NOCMAX,
     %                   CONTPR, DIRPRI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES CONTRAINTES ET DIRECTIONS PRINCIPALES
C -----    EN LES POINTS D'INTEGRATION D'UN EF D'ELASTICITE
C
C ENTREES :
C ---------
C NUELEM : NUMERO DE L'ELEMENT FINI POUR CE TYPE D'ELEMENT
C NBELEM : NOMBRE TOTAL D'ELEMENTS FINIS DE CE TYPE
C NDIM   : DIMENSION DE L'ESPACE 2 OU 3
C NCOMP  : NOMBRE DE COMPOSANTES DU TENSEUR DES CONTRAINTES EN UN POINT
C          ( 2D => 3 ,  2D AXISYMETRIQUE => 4 ,  3D => 6 )
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DE L ELEMENT FINI
C NDSM   : NOMBRE DE CAS DE CHARGE OU SECONDS MEMBRES DU SYSTEME
C NOFORE : NUMERO DE LA FONCTION 'REGION' 0 SI ABSENTE
C
C COOPTC : LES COORDONNEES DES POINTS DE CALCUL DES CONTRAINTES
C STRESS : LES CONTRAINTES AUX POINTS ET POUR CHAQUE CAS DE CHARGE
C          STRESS ( NCOMP , NPI , NDSM )
C
C MODIFIES:
C ---------
C CONTMI : LA CONTRAINTE MINIMALE DE L'EF NUELEM
C NOPMIN : LE NUMERO DU POINT DE LA CONTRAINTE MINIMALE
C NOCMIN : LE NUMERO DU CAS   DE LA CONTRAINTE MINIMALE
C
C CONTMX : LA CONTRAINTE MAXIMALE DE L'EF NUELEM
C NOPMAX : LE NUMERO DU POINT DE LA CONTRAINTE MAXIMALE
C NOCMAX : LE NUMERO DU CAS   DE LA CONTRAINTE MAXIMALE
C
C CONTPR : LES CONTRAINTES PRINCIPALES
C DIRPRI : LES DIRECTIONS  PRINCIPALES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456...............................................................012
      PARAMETER        (MXEFAF=1)
C     NOMBRE MAXIMAL D'EF POUR L'AFFICHAGE DES CONTRAINTES PRINCIPALES
      include"./incl/ctemps.inc"
      include"./incl/langue.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  CONTPR(NBELEM,NPI,NDIM,NDSM), CONTMI, CONTMX
      REAL              DIRPRI(NBELEM,NPI,NDIM,NDIM,NDSM)
      REAL              COOPTC(NBELEM,NPI,NDIM)
      DOUBLE PRECISION  STRESS(NCOMP,NPI,NDSM)
      INTRINSIC         REAL
C
      DOUBLE PRECISION  VP(3),VP1,VP2,VP3,VMIN,VMAX,
     %                  DIRP(3,3),XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,
     %                  A(9)
      EQUIVALENCE       (VP(1),VP1),(VP(2),VP2),(VP(3),VP3)
      EQUIVALENCE       (DIRP(1,1),XV1),(DIRP(1,2),XV2),(DIRP(1,3),XV3),
     %                  (DIRP(2,1),YV1),(DIRP(2,2),YV2),(DIRP(2,3),YV3),
     %                  (DIRP(3,1),ZV1),(DIRP(3,2),ZV2),(DIRP(3,3),ZV3)
C
      CONTMI =  1D38
      CONTMX = -1D38
      DO 100 L=1,NPI
C
C        LES CONTRAINTES EN DIMENSION 3 OU EN DIMENSION 2 OU AXISYMETRIQUE
C        =================================================================
         DO 90 M=1,NDSM
C
            IF( NCOMP .EQ. 6 ) THEN
C
C              NDIM=3 STRESS= SX SY SZ SXY SYZ SZX
C              -----------------------------------
C              AXX AXY AXZ
C              AYX AYY AYZ
C              AZX AZY AZZ
C
               A(1) = STRESS(1,L,M)
               A(2) = STRESS(4,L,M)
               A(3) = STRESS(6,L,M)
C
               A(4) = STRESS(4,L,M)
               A(5) = STRESS(2,L,M)
               A(6) = STRESS(5,L,M)
C
               A(7) = STRESS(6,L,M)
               A(8) = STRESS(5,L,M)
               A(9) = STRESS(3,L,M)
C
C              CALCUL DES VALEURS ET VECTEURS PROPRES
               CALL VPM3DD( A, VP, DIRP )
C
C              LES CONTRAINTES PRINCIPALES
               CONTPR(NUELEM,L,1,M) = VP1
               CONTPR(NUELEM,L,2,M) = VP2
               CONTPR(NUELEM,L,3,M) = VP3
C
               VMIN = MIN( VP1, VP2, VP3 )
               IF( VMIN .LT. CONTMI ) THEN
                  CONTMI = VMIN
                  NOPMIN = L
                  NOCMIN = M
               ENDIF
C
               VMAX = MAX( VP1, VP2, VP3 )
               IF( VMAX .GT. CONTMX ) THEN
                  CONTMX = VMAX
                  NOPMAX = L
                  NOCMAX = M
               ENDIF
C
C              LES DIRECTIONS PRINCIPALES
C              LA PREMIERE
               DIRPRI(NUELEM,L,1,1,M) = REAL( XV1 )
               DIRPRI(NUELEM,L,2,1,M) = REAL( YV1 )
               DIRPRI(NUELEM,L,3,1,M) = REAL( ZV1 )
C              LA SECONDE
               DIRPRI(NUELEM,L,1,2,M) = REAL( XV2 )
               DIRPRI(NUELEM,L,2,2,M) = REAL( YV2 )
               DIRPRI(NUELEM,L,3,2,M) = REAL( ZV2 )
C              LA TROISIEME
               DIRPRI(NUELEM,L,1,3,M) = REAL( XV3 )
               DIRPRI(NUELEM,L,2,3,M) = REAL( YV3 )
               DIRPRI(NUELEM,L,3,3,M) = REAL( ZV3 )
C
               IF( NUELEM .LE. MXEFAF ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10030) NUELEM,L,
     %                    (COOPTC(NUELEM,L,J),J=1,3),M,VP,DIRP
                  ELSE
                     WRITE(IMPRIM,20030) NUELEM,L,
     %                    (COOPTC(NUELEM,L,J),J=1,3),M,VP,DIRP
                  ENDIF
               ENDIF
10030 FORMAT(' EF',I5,': PT',I3,' X=',G12.4,' Y=',G12.4,' Z=',G12.4/
     %' CAS',I4,
     %': CP1',G12.4,' CP2',G12.4,' CP3',G12.4,
     % ' DP1',3F7.3,' DP2',3F7.3,' DP3',3F7.3)
20030 FORMAT(' FE',I6,': PT',I3,' X=',G12.4,' Y=',G12.4,' Z=',G12.4/
     %' CASE',I4,
     %': MS1',G12.4,' MS2',G12.4,' MS3',G12.4,
     % ' MD1',3F7.3,' MD2',3F7.3,' MD3',3F7.3)
C
            ELSE
C
               IF ( NCOMP .EQ. 4 ) THEN
C
C                 CAS AXISYMETRIQUE
C                 STRESS= S(Z) S(R) S(TETA) S(RZ) => A= A(R) A(RZ) A(Z)
C                 -----------------------------------------------------
                  A(1) = STRESS(2,L,M)
                  A(2) = STRESS(4,L,M)
                  A(3) = STRESS(1,L,M)
C
               ELSE IF( NCOMP .EQ. 3 ) THEN
C
C                 NDIM=2 STRESS= SX SY SXY => A= A(X) A(XY) A(Y)
C                 ----------------------------------------------
                  A(1) = STRESS(1,L,M)
                  A(2) = STRESS(3,L,M)
                  A(3) = STRESS(2,L,M)
               ENDIF
C
C              CALCUL DES VALEURS ET VECTEURS PROPRES DE LA MATRICE A
               CALL VPM2DD( A, VP1,VP2, XV1,YV1 )
C
C              LES CONTRAINTES PRINCIPALES
               CONTPR(NUELEM,L,1,M) = VP1
               CONTPR(NUELEM,L,2,M) = VP2
C
               VMIN = MIN( VP1, VP2 )
               IF( VMIN .LT. CONTMI ) THEN
                  CONTMI = VMIN
                  NOPMIN = L
                  NOCMIN = M
               ENDIF
C
               VMAX = MAX( VP1, VP2 )
               IF( VMAX .GT. CONTMX ) THEN
                  CONTMX = VMAX
                  NOPMAX = L
                  NOCMAX = M
               ENDIF
C
C              LES DIRECTIONS PRINCIPALES
C              LA PREMIERE
               DIRPRI(NUELEM,L,1,1,M) = REAL( XV1 )
               DIRPRI(NUELEM,L,2,1,M) = REAL( YV1 )
C              LA SECONDE
               DIRPRI(NUELEM,L,1,2,M) = REAL( -YV1 )
               DIRPRI(NUELEM,L,2,2,M) = REAL(  XV1 )
C
C              IMPRESSION REDUITE DES CONTRAINTES
               J = 0
               IF( NOFORE .LE. 0 .AND. NUELEM .LE. MXEFAF ) J=1
               IF( NOFORE .GT. 0 ) THEN
C
C                 LES 4 PARAMETRES D'APPEL DE LA FONCTION 'CHOIX'
C                 LE TEMPS EN 1-ER PARAMETRE
                  A(1) = TEMPS
C                 PUIS LES 3 COORDONNEES X Y Z DU NOEUD
                  A(2) = COOPTC(NUELEM,L,1)
                  A(3) = COOPTC(NUELEM,L,2)
                  A(4) = COOPTC(NUELEM,L,3)
C                 FONCTION CHOIX(TEMPS,X,Y,Z)
                  CALL FONVAL( NOFORE, 4, A, NCODEV, VMAX )
                  IF( NINT(VMAX) .NE. 0 ) J=1
C
                ENDIF
                IF( J .NE. 0 ) THEN
                  IF( NCOMP .EQ. 4 ) THEN
C
C                    AXISYMETRIQUE
                     IF( LANGAG .EQ . 0 ) THEN
                        WRITE(IMPRIM,10025) NUELEM,L,
     %                       (COOPTC(NUELEM,L,J),J=1,2),M,
     %                       VP1,VP2,XV1,YV1,-YV1,XV1
                     ELSE
                        WRITE(IMPRIM,20025) NUELEM,L,
     %                       (COOPTC(NUELEM,L,J),J=1,2),M,
     %                       VP1,VP2,XV1,YV1,-YV1,XV1
                     ENDIF
10025 FORMAT(' EF',I5,' : PT',I3,' R=',G12.4,' Z=',G12.4,
     %' CAS',I4,': CP1',G12.4,' CP2',G12.4,
     %           '  DP1',2F7.3,' DP2',2F7.3)
20025 FORMAT(' FE',I5,' : PT',I3,' R=',G12.4,' Z=',G12.4,
     %' CASE',I4,': MS1',G12.4,' MS2',G12.4,
     %           '  MD1',2F7.3,' MD2',2F7.3)
C
                  ELSE
C
C                    2D
                     IF( LANGAG .EQ . 0 ) THEN
                        WRITE(IMPRIM,10020) NUELEM,L,
     %                       (COOPTC(NUELEM,L,J),J=1,2),M,
     %                       VP1,VP2,XV1,YV1,-YV1,XV1
                     ELSE
                        WRITE(IMPRIM,20020) NUELEM,L,
     %                       (COOPTC(NUELEM,L,J),J=1,2),M,
     %                       VP1,VP2,XV1,YV1,-YV1,XV1
                     ENDIF
10020 FORMAT(' EF',I5,' : PT',I3,' X=',G12.4,' Y=',G12.4,
     %' CAS',I4,': CP1',G12.4,' CP2',G12.4,
     %           '  DP1',2F7.3,' DP2',2F7.3)
20020 FORMAT(' FE',I5,' : PT',I3,' X=',G12.4,' Y=',G12.4,
     %' CASE',I4,': MS1',G12.4,' MS2',G12.4,
     %           '  MD1',2F7.3,' MD2',2F7.3)
C
                  ENDIF
               ENDIF
            ENDIF
  90     CONTINUE
 100  CONTINUE
C
      RETURN
      END
