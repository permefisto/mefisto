      SUBROUTINE RETRIAXY1( XP0, YP0, XP1, YP1, MNNPEF, MNXYZP,
     %                      MOARET, MXARET, MNLARE,
     %                      NEF, CB1, CB2, CB3, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:  RECHERCHE EXHAUSTIVE DU TRIANGLE NEF CONTENANT LE POINT XYP1
C ----  PAR UN CHEMIN XYP0-XYP1 A TRAVERS LES ARETES DU TRIANGLE NEF
C
C ENTREES:
C --------
C XP0,YP0: COORDONNEES DU POINT INITIAL DU SEGMENT DE DROITE
C XP1,YP1: COORDONNEES DU POINT FINAL   DU SEGMENT DE DROITE
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNXYZP : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MOARET : LE NOMBRE D'ENTIERS PAR ARETE
C MXARET : LA MAJORATION DU NOMBRE D'ARETES
C MNLARE : ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE
C NEF    : >0 NUMERO DE L'EF DE DEPART POUR TROVER L'EF CONTENANT XYP1
C
C SORTIES:
C --------
C NEF    : >0 NUMERO DE L'EF CONTENANT XYP1 SI IERR=0
C          =0 SI PARTICULE SORTANT DU MAILLAGE XYP1 EST EXTERIEUR
C CB1,CB2,CB3: COORDONNEES BARYCENTRIQUES DE XYP1 DANS LE TRIANGLE NEF
C IERR   : =0 PAS D'ERREUR, NEF>0 EST LE NUMERO DE L'EF CONTENANT XYP1
C                           NEF=0 SI PAS D'EF CONTENANT XYP1
C          =1 ARETE NON RETROUVEE POUR TRAVERSER LE TRIANGLE
C          =2 TRIANGLE D'AIRE <=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY NOVEMBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/nctyef.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      DOUBLE PRECISION  XP0, YP0, XP1, YP1,
     %                  X1,X2,X3, Y1,Y2,Y3, X0,Y0, CBMIN, CB1,CB2,CB3, D
      INTEGER           NOSOTR(3), NGS(2)
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
      NTYEF = 1
      MNELE = MCN( MNNPEF -1 + NTYEF )
C
C     LE NOMBRE D'ELEMENTS FINIS
      NBELEM = MCN( MNELE + WBELEM )
C
C     LE NO DES NOEUDS DU TRIANGLE NEF
 10   MNN = MNELE + WUNDEL -1 + NEF
      DO J=1,3
         NOSOTR(J) = MCN(MNN)
         MNN = MNN + NBELEM
      ENDDO
C
C     LES 3 SOMMETS SONT EN TETE DES NOEUDS DE L'EF
      MNP= MNXYZP + WYZPOI -3
      MN = MNP + 3*NOSOTR(1)
      X1 = RMCN(MN)
      Y1 = RMCN(MN+1)
C
      MN = MNP + 3*NOSOTR(2)
      X2 = RMCN(MN)
      Y2 = RMCN(MN+1)
C
      MN = MNP + 3*NOSOTR(3)
      X3 = RMCN(MN)
      Y3 = RMCN(MN+1)
C
C     2 FOIS LA SURFACE DU TRIANGLE = DETERMINANT DE LA MATRICE
      D = ( X2 - X1 ) * ( Y3 - Y1 ) - ( X3 - X1 ) * ( Y2 - Y1 )
      IF( D .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'TRIANGLE',NEF,'d''AIRE',D,'<=0 !'
         ELSE
            WRITE(IMPRIM,*) 'TRIANGLE',NEF,'with a SURFACE',D,'<=0 !'
         ENDIF
         IERR = 2
         GOTO 9000
      ENDIF
C
C     CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU
C     POINT XP1 YP1 DANS LE TRIANGLE NEF
      CB1 = ( (X2-XP1)*(Y3-YP1) - (X3-XP1)*(Y2-YP1) ) / D
      CB2 = ( (X3-XP1)*(Y1-YP1) - (X1-XP1)*(Y3-YP1) ) / D
      CB3 = ( (X1-XP1)*(Y2-YP1) - (X2-XP1)*(Y1-YP1) ) / D
C
      IF( CB1 .GE. 0D0 .AND. CB1 .LE. 1D0 .AND.
     %    CB2 .GE. 0D0 .AND. CB2 .LE. 1D0 .AND.
     %    CB3 .GE. 0D0 .AND. CB3 .LE. 1D0 ) THEN
C
C        LE POINT XYP1 EST DANS LE TRIANGLE NEF OU SUR UNE DE SES 3 ARETES
C        -----------------------------------------------------------------
         IERR = 0
         GOTO 9000
C
      ENDIF
C
C     LE POINT XP1 YP1 EST EXTERIEUR AU TRIANGLE NEF
C     PARCOURS PAR LES ARETES DE NEF POUR TROUVER LE TRIANGLE CONTENANT XY1
C     ---------------------------------------------------------------------
      CBMIN = MIN( CB1, CB2, CB3 )
C
      IF( CBMIN .EQ. CB3 ) THEN
C        INTERSECTION DU VECTEUR VXYP0 ISSU DE XYP0
C        AVEC L'ARETE 1 DU TRIANGLE
         CALL INTARSE( X1,Y1, X2,Y2, XP0,YP0, XP1,YP1,
     %                 LINTER, X0,Y0 )
C        LINTER : -1 SI (XP0,YP0)-(XP1,YP1) PARALLELE A (X1,Y1)-(X2,Y2)
C                  0 SI (XP0,YP0)-(XP1,YP1) N'INTERSECTE PAS (X1,Y1)-(X2,Y2)
C                                           ENTRE CES 2 SOMMETS
C                  1 SI (XP0,YP0)-(XP1,YP1)   INTERSECTE     (X1,Y1)-(X2,Y2)
C                                           ENTRE CES 2 SOMMETS
C        X0,Y0  :  2 COORDONNEES DU POINT D'INTERSECTION DU SEGMENT
C                    (XP0,YP0)-(XP1,YP1)
C        SUR L'ARETE (X1,Y1)-(X2,Y2)
         N1 = NOSOTR(1)
         N2 = NOSOTR(2)
         GOTO 60
      ENDIF
C
      IF( CBMIN .EQ. CB1 ) THEN
C        INTERSECTION DU VECTEUR VXYP0 ISSU DE XYP0
C        AVEC L'ARETE 2 DU TRIANGLE
         CALL INTARSE( X2,Y2, X3,Y3, XP0,YP0, XP1,YP1,
     %                 LINTER, X0,Y0 )
         N1 = NOSOTR(2)
         N2 = NOSOTR(3)
         GOTO 60
      ENDIF
C
      IF( CBMIN .EQ. CB2 ) THEN
C        INTERSECTION DU VECTEUR VXYP0 ISSU DE XYP0
C        AVEC L'ARETE 3 DU TRIANGLE
         CALL INTARSE( X3,Y3, X1,Y1, XP0,YP0, XP1,YP1,
     %                 LINTER, X0,Y0 )
         N1 = NOSOTR(3)
         N2 = NOSOTR(1)
         GOTO 60
      ENDIF
C
 60   IF( LINTER .LT. 0 ) THEN
C        DIRECTIONS PARALLELES => LA PARTICULE RESTE DANS L'EF NEF
         GOTO 10
      ENDIF
C
C     RECHERCHE DU TRIANGLE DE L'AUTRE COTE DE L'ARETE
      IF( N1 .LT. N2 ) THEN
         NGS(1) = N1
         NGS(2) = N2
      ELSE
         NGS(1) = N2
         NGS(2) = N1
      ENDIF
      LIBREF = MXARET
      CALL HACHAG( 2, NGS, MOARET, MXARET, MCN(MNLARE), 3,
     &             LIBREF, NOAR )
      IF( NOAR .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:24),'(2I12)') N1,N2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'RETRIAXY1: ARETE PERDUE'
            KERR(2) = 'SOMMETS ' // KERR(MXLGER)(1:24)
         ELSE
            KERR(1) = 'RETRIAXY1: LOST EDGE'
            KERR(2) = 'VERTICES ' // KERR(MXLGER)(1:24)
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF
C
C     ADRESSE DE L'ARETE NOAR DANS LARETE
      MNA = MNLARE + MOARET * (NOAR-1) - 1
C     LE CODAGE DU NUMERO DE TYPE D'EF POUR LE NO D'EF DANS LARETE
      NUCYEF = NCTYEF * NTYEF
      NEF1 = ABS( MCN( MNA + 4 ) )
      IF( NEF1 .GT. NUCYEF ) NEF1 = NEF1 - NUCYEF
      NEF2 = ABS( MCN( MNA + 5 ) )
      IF( NEF2 .GT. NUCYEF ) NEF2 = NEF2 - NUCYEF
C
C     L'EF DE L'AUTRE COTE DE L'ARETE
      IF( NEF1 .EQ. NEF ) THEN
         NEF = NEF2
      ELSE
         NEF = NEF1
      ENDIF
C
      IF( NEF .LE. 0 ) THEN
C
C        PAS DE TRIANGLE ADJACENT => LA PARTICULE SORT DU MAILLAGE
C        ---------------------------------------------------------
         NEF  = 0
         IERR = 0
         GOTO 9000
C
      ENDIF
C
C     PASSAGE A L'EF NEF
      GOTO 10
C
 9000 RETURN
      END
