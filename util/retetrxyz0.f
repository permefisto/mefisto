      SUBROUTINE RETETRXYZ0( XYZP0, MNNPEF, MNXYZP,
     %                       NEF, CB1, CB2, CB3, CB4, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    RECHERCHE EXHAUSTIVE DU TETRAEDRE NEF CONTENANT LE POINT XYZP0
C ----    PAR UNE BOUCLE SUR TOUS LES TETRAEDRES DU MAILLAGE
C         SI AUCUN TETRAEDRE NE CONTIENT XYZ0, ALORS XYZ0 DEVIENT
C         LE BARYCENTRE DU TETRAEDRE LE PLUS PROCHE DE LUI

C ENTREES:
C --------
C XYZP0  : 3 COORDONNEES DU POINT
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNXYZP : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET

C SORTIES:
C --------
C NEF    : NUMERO DE L'EF CONTENANT XYZP0 SI IERR=0
C CB1,CB2,CB3,CB4: COORDONNEES BARYCENTRIQUES DE XYZP0 DANS LE TETRAEDRE NEF
C IERR   : =0 PAS D'ERREUR, NEF>0 EST LE NUMERO DE L'EF CONTENANT XYZP0
C          =1 PAS DE TETRAEDRE CONTENANT LE POINT XYZP0
C          =2 EXISTENCE D'UN TETRAEDRE D'AIRE <=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY NOVEMBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      REAL              XYZP0(3), BARY(3), BARYMIN(3), DIST, DISTMIN
      DOUBLE PRECISION  VOLTER, CB1, CB2, CB3, CB4, D

C     L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
      NTYEF = 1
      MNELE = MCN( MNNPEF -1 + NTYEF )

C     LE NOMBRE D'ELEMENTS FINIS
      NBELEM = MCN( MNELE + WBELEM )

C     RECHERCHE DU POINT XYZP0 DANS LES EF
C     ------------------------------------
      MNP = MNXYZP + WYZPOI -3
      DO 100 NEF = 1, NBELEM

C        LE NO DES NOEUDS DU TETRAEDRE NEF
         MNN = MNELE  + WUNDEL -1 + NEF
C        LES 4 SOMMETS SONT EN TETE DES NOEUDS DE L'EF
         MNS1 = MNP + 3 * MCN(MNN)
         MNN  = MNN + NBELEM
         MNS2 = MNP + 3 * MCN(MNN)
         MNN  = MNN + NBELEM
         MNS3 = MNP + 3 * MCN(MNN)
         MNN  = MNN + NBELEM
         MNS4 = MNP + 3 * MCN(MNN)
         MNN  = MNN + NBELEM

C        LE VOLUME DU TETRAEDRE P1 P2 P3 P4
C        CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT P
         D = VOLTER( RMCN(MNS1),RMCN(MNS2),RMCN(MNS3),RMCN(MNS4) )

         IF( D .GT. 0 ) THEN

C           TETRAEDRE DE VOLUME NON DEGENERE
C           CALCUL DES 4 COORDONNEES BARYCENTRIQUES DU
C           POINT XYZP0 DANS LE TETRAEDRE NEF
            CB1 = VOLTER(XYZP0,RMCN(MNS2),RMCN(MNS3),RMCN(MNS4)) / D
            CB2 = VOLTER(XYZP0,RMCN(MNS3),RMCN(MNS1),RMCN(MNS4)) / D
            CB3 = VOLTER(XYZP0,RMCN(MNS1),RMCN(MNS2),RMCN(MNS4)) / D
            CB4 = VOLTER(XYZP0,RMCN(MNS1),RMCN(MNS3),RMCN(MNS2)) / D
            IF( CB1 .GE. -1D-8 .AND. CB1 .LE. 1D0 .AND.
     %          CB2 .GE. -1D-8 .AND. CB2 .LE. 1D0 .AND.
     %          CB3 .GE. -1D-8 .AND. CB3 .LE. 1D0 .AND.
     %          CB4 .GE. -1D-8 .AND. CB4 .LE. 1D0 ) THEN
C              LE POINT EST DANS LE TETRAEDRE NEF
C              OU SUR UNE DE SES 4 FACES
               IERR = 0
               GOTO 9000
            ENDIF

         ELSE

C           TETRAEDRE DEGENERE
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'retetrxyz0: TETRAEDRE',NEF,'de VOLUME',D,'<0 !'
            ELSE
              PRINT*,'retetrxyz0: TETRAHEDRON',NEF,' of VOLUME',D,'<0 !'
            ENDIF
            IERR = 2
            GOTO 100

         ENDIF

 100  ENDDO

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'retetrxyz0: AUCUN ELEMENT FINI CONTIENT le POINT X=',
     %           XYZP0(1),' Y=',XYZP0(2),' Z=',XYZP0(3)
      ELSE
         PRINT*,'retetrxyz0: NO FINITE ELEMENT CONTAINS THIS POINT X=',
     %           XYZP0(1),' Y=',XYZP0(2),' Z=',XYZP0(3)
      ENDIF
      IERR = 1

C     RECHERCHE DU TETRAEDRE DE BARYCENTRE LE PLUS PROCHE DU POINT
C     ------------------------------------------------------------
      DISTMIN = 1E28
      MNP = MNXYZP + WYZPOI -3
      DO 200 NEF = 1, NBELEM

C        LE NO DES NOEUDS DU TETRAEDRE NEF
         MNN = MNELE  + WUNDEL -1 + NEF
C        LES 4 SOMMETS SONT EN TETE DES NOEUDS DE L'EF
         MNS1 = MNP + 3 * MCN(MNN) -1
         MNN  = MNN + NBELEM
         MNS2 = MNP + 3 * MCN(MNN) -1
         MNN  = MNN + NBELEM
         MNS3 = MNP + 3 * MCN(MNN) -1
         MNN  = MNN + NBELEM
         MNS4 = MNP + 3 * MCN(MNN) -1
         MNN  = MNN + NBELEM

         DIST = 0
         DO K=1,3
            BARY( K ) = ( RMCN(MNS1+K) + RMCN(MNS2+K)
     %                  + RMCN(MNS3+K) + RMCN(MNS4+K) ) / 4
            DIST = DIST + ( XYZP0( K ) - BARY( K ) ) **2
         ENDDO

         IF( DIST .LT. DISTMIN ) THEN
            DISTMIN = DIST
            NEFMIN  = NEF
            BARYMIN(1) = BARY(1)
            BARYMIN(2) = BARY(2)
            BARYMIN(3) = BARY(3)
         ENDIF

 200  ENDDO

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'retetrxyz0: Le BARYCENTRE X=',BARYMIN(1),
     %          ' Y=',BARYMIN(2),' Z=',BARYMIN(3),
     %          ' du TETRAEDRE',NEFMIN
         PRINT*,'retetrxyz0: EST a la DISTANCE MINIMALE',
     %           SQRT(DISTMIN),' du POINT X=',
     %           XYZP0(1),' Y=',XYZP0(2),' Z=',XYZP0(3),
     %           ' et le REMPLACE'
      ELSE
         PRINT*,'retetrxyz0: the BARYCENTRE X=',BARYMIN(1),
     %          ' Y=',BARYMIN(2),' Z=',BARYMIN(3),
     %          ' of TETRAHEDRON',NEFMIN
         PRINT*,'retetrxyz0: is at DISTANCE MIN',SQRT(DISTMIN),
     %          ' of POINT X=',
     %           XYZP0(1),' Y=',XYZP0(2),' Z=',XYZP0(3),
     %          ' and REPLACE IT'
      ENDIF

      XYZP0(1) = BARYMIN(1)
      XYZP0(2) = BARYMIN(2)
      XYZP0(3) = BARYMIN(3)
      CB1      = 0.25D0
      CB2      = 0.25D0
      CB3      = 0.25D0
      CB4      = 0.25D0
      NEF      = NEFMIN
      IERR = 0

 9000 RETURN
      END
