      SUBROUTINE RETRIAXY0( XP0, YP0, MNNPEF, MNXYZP,
     %                      NEF, CB1, CB2, CB3, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:  RECHERCHE EXHAUSTIVE DU TRIANGLE NEF CONTENANT LE POINT XYP0
C ----  PAR UNE BOUCLE SUR TOUS LES TRIANGLES DU MAILLAGE
C
C ENTREES:
C --------
C XP0,YP0: COORDONNEES DU POINT
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNXYZP : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C
C SORTIES:
C --------
C NEF    : NUMERO DE L'EF CONTENANT XYP0 SI IERR=0
C CB1,CB2,CB3: COORDONNEES BARYCENTRIQUES DE XYP0 DANS LE TRIANGLE NEF
C IERR   : =0 PAS D'ERREUR, NEF>0 EST LE NUMERO DE L'EF CONTENANT XYP0
C          =1 PAS DE TRIANGLE CONTENANT LE POINT XYP0
C          =2 TRIANGLE D'AIRE <=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY NOVEMBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      DOUBLE PRECISION  XP0, YP0, X1,X2,X3, Y1,Y2,Y3, CB1,CB2,CB3, D
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
      NTYEF = 1
      MNELE = MCN( MNNPEF -1 + NTYEF )
C
C     LE NOMBRE D'ELEMENTS FINIS
      NBELEM = MCN( MNELE + WBELEM )
C
C     RECHERCHE DU POINT XYP0 DANS LES EF
      MNP = MNXYZP + WYZPOI -3
      DO 20 NEF = 1, NBELEM
C
C        LE NO DES NOEUDS DU TRIANGLE NEF
         MNN = MNELE  + WUNDEL -1 + NEF
C
C        LES 3 SOMMETS SONT EN TETE DES NOEUDS DE L'EF
         N1  = MCN(MNN)
         MNN = MNN + NBELEM
         MN  = MNP + 3*N1
         X1  = RMCN(MN)
         Y1  = RMCN(MN+1)
C
         N2  = MCN(MNN)
         MNN = MNN + NBELEM
         MN  = MNP + 3*N2
         X2  = RMCN(MN)
         Y2  = RMCN(MN+1)
C
         N3  = MCN(MNN)
         MNN = MNN + NBELEM
         MN  = MNP + 3*N3
         X3  = RMCN(MN)
         Y3  = RMCN(MN+1)
C
C        2 FOIS LA SURFACE DU TRIANGLE = DETERMINANT DE LA MATRICE
C        DE CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT P
         D = ( X2 - X1 ) * ( Y3 - Y1 ) - ( X3 - X1 ) * ( Y2 - Y1 )
C
         IF( D .GT. 0 ) THEN
C           TRIANGLE NON DEGENERE (CAS DEGENERE TRAITE dans ptdatr.f)
C           CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU
C           POINT XP0 YP0 DANS LE TRIANGLE NEF
            CB1 = ( (X2-XP0)*(Y3-YP0) - (X3-XP0)*(Y2-YP0) ) / D
            CB2 = ( (X3-XP0)*(Y1-YP0) - (X1-XP0)*(Y3-YP0) ) / D
            CB3 = ( (X1-XP0)*(Y2-YP0) - (X2-XP0)*(Y1-YP0) ) / D
            IF( CB1 .GE. 0D0 .AND. CB1 .LE. 1D0 .AND.
     %          CB2 .GE. 0D0 .AND. CB2 .LE. 1D0 .AND.
     %          CB3 .GE. 0D0 .AND. CB3 .LE. 1D0 ) THEN
C              LE POINT EST DANS LE TRIANGLE NEF OU SUR UNE DE SES 3 ARETES
               IERR = 0
               GOTO 9000
            ENDIF
         ELSE
C           TRIANGLE DEGENERE
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*)'TRIANGLE',NEF,'d''AIRE',D,'<0 !'
            ELSE
               WRITE(IMPRIM,*)'TRIANGLE',NEF,'with a SURFACE',D,'<0 !'
            ENDIF
            IERR = 2
            GOTO 9000
         ENDIF
 20   CONTINUE
C
      NBLGRC(NRERR) = 2
      WRITE(KERR(2)(1:35),'(A2,G15.6,A3,G15.6)') 'X=',XP0,' Y=',YP0
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='AUCUN ELEMENT FINI CONTIENT LE POINT'
      ELSE
         KERR(1) ='NO FINITE ELEMENT CONTAINS THE POINT'
      ENDIF
      CALL LEREUR
      IERR = 1
C
 9000 RETURN
      END
