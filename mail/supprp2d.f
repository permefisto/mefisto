       SUBROUTINE SUPPRP2D( NX, NY, MNXYZ1, MNNSEF1,  NUMST,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER UN SOMMET DU MAILLAGE D'UNE SURFACE 2D
C -----    SI LE SOMMET EST INTERNE TRIANGULATION SANS LE SOMMET
C          SI LE SOMMET EST FRONTALIER SUPPRESSION DES TRIANGLES DE CE SOMMET
C ENTREES:
C --------
C NX     : ABSCISSE DU POINT A SUPPRIMER
C NY     : ORDONNEE DU POINT A SUPPRIMER
C MNXYZ1 : ADRESSE DU TABLEAU XYZSOMMET
C MNNSEF1: ADRESSE DU TABLEAU NSEF
C
C SORTIES:
C --------
C NUMST  : LE NUMERO DU SOMMET SUPPRIME
C IERR   : 0 SI PAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: LEVI-OVSIOUK DEA ANALYSE NUMERIQUE UPMC PARIS    JANVIER 2000
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NUMST,N1,N2
      REAL          COORD(1:3)
      REAL          CN1(1:2)
      REAL          CN2(1:2)
      REAL          CNPT(1:2)
C
      IERR   = 0
      MNNTAB = 0
      MSUPTR = 0

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZ1+WNBSOM)

C     NBEF NOMBRE DES EF
      NBEF   = MCN(MNNSEF1+WBEFOB)
      MNSOEL = MNNSEF1 + WUSOEF

C     COORDONNEES OBJET DU POINT CLIQUE
C     =================================
      X=XOB2PX(NX)
      Y=YOB2PX(NY)
C
C      RECHERCHE DU SOMMET LE PLUS PROCHE NUMST
C      ========================================
       COORD(1)=RMCN(MNXYZ1+WYZSOM)
       COORD(2)=RMCN(MNXYZ1+WYZSOM+1)
       COORD(3)=RMCN(MNXYZ1+WYZSOM+2)
       DISTANCE=SQRT((X-COORD(1))**2+(Y-COORD(2))**2)
       NUMST=1
       DISTMP=DISTANCE
C      ON BOUCLE SUR LES POINTS DU MAILLAGE
       DO 10 I=1,NBSOM
         COORD(1)=RMCN(MNXYZ1+WYZSOM+3*(I-1))
         COORD(2)=RMCN(MNXYZ1+WYZSOM+3*(I-1)+1)
         COORD(3)=RMCN(MNXYZ1+WYZSOM+3*(I-1)+2)
         DISTMP=SQRT((X-COORD(1))**2+(Y-COORD(2))**2)
         IF (DISTMP .LE. DISTANCE) THEN
           DISTANCE=DISTMP
           NUMST=I
         ENDIF
10     CONTINUE
C
C      RECHERCHE DES QUADRANGLES OU TRIANGLES LE CONTENANT
C      ===================================================
       CALL TNMCDC( 'ENTIER', NBSOM, MNNTAB )
       IF (MNNTAB.EQ.0) THEN
          NBLGRC(NRERR) = 1
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'PAS ASSEZ DE PLACE EN MCN'
          ELSE
             KERR(1) = 'NOT ENOUGH MEMORY IN MCN'
          ENDIF
          CALL LEREUR
          IERR = 1
          RETURN
       ENDIF
       NBQTAJ=0
       DO 20 J=1,NBEF
         NTEMP=MNSOEL+4*(J-1)
         IF( (MCN(NTEMP)  .EQ.NUMST).OR.(MCN(NTEMP+1).EQ.NUMST)
     %   .OR.(MCN(NTEMP+2).EQ.NUMST).OR.(MCN(NTEMP+3).EQ.NUMST)) THEN
           MCN(MNNTAB+NBQTAJ)=J
           NBQTAJ=NBQTAJ+1
         ENDIF
20     CONTINUE
C      
C     ON INITIALISE LE TABLEAU DES TRIANGLES A SUPPRIMER
C     ==================================================
      LSUPTR = NBQTAJ
      IF( LSUPTR .LE. 0 ) LSUPTR=1
      CALL TNMCDC( 'ENTIER', LSUPTR, MSUPTR )
      NSUPTR=0
      DO 25 I=1,NBQTAJ
         MCN(MSUPTR+I-1)=0
 25   CONTINUE
C
C      ON DETERMINE SI NUMST EST AU BORD
C      =================================
       CALL PFRONT( MNNTAB,NBQTAJ,MNNSEF1,MNXYZ1,NUMST,NFR )
C      NFR=1 SI A LA FRONTIERE
       IF (NFR.EQ.1) THEN
C        ON EST A LA FRONTIERE
C        DONC TOUS LES TRIANGLES DOIVENT ETRE SUPPRIMES
         DO 27 I=1,NBQTAJ
           MCN(MSUPTR+I-1)=1
27       CONTINUE
         NSUPTR=NBQTAJ
       ELSE
C        ON EST PAS A LA FRONTIERE
         NSUPTR=2      
C
C        ON CHOISIT UN POINT NPT DIFFERENT DE NUMST
C        ==========================================
C        IL FAUDRAIT TROUVER NPT TEL QUE LA FIGURE
C        SOIT ETOILEE PAR RAPPORT A NPT
       CALL ETOILED( NUMST,MNNTAB,NBQTAJ,MNXYZ1,MNNSEF1,NPT )
       IF ( NPT .EQ. 0 ) GOTO 500
C
C        ON BOUCLE SUR LES TRIANGLES
C        ===========================
         NN=0
         DO 50 I=1,NBQTAJ
           NUMTR=MCN(MNNTAB+I-1)
           NTEMP=MNNSEF1+WUSOEF+4*(NUMTR-1)
           IF ((MCN(NTEMP).EQ.NPT).OR.(MCN(NTEMP+1).EQ.NPT)
     %       .OR.(MCN(NTEMP+2).EQ.NPT)
     %       .OR.(MCN(NTEMP+3).EQ.NPT)) THEN
              IF (NN.EQ.0) THEN
C                IL FAUDRA SUPPRIMER NN1
                 NN1=NUMTR
                 MCN(MSUPTR+I-1)=1
                 NN=1
              ELSE
C                IL FAUDRA SUPPRIMER NN2
                 NN2=NUMTR
                 MCN(MSUPTR+I-1)=1
              ENDIF
           ELSE
C             ON NE CONTIENT PAS LE POINT NPT
C             ON REMPLACE NUMST PAR NPT SI LE TRIANGLE
C             RESULTANT N'EST PAS PLAT SINON ON SUPPRIME
C             ON RECUPERE LES DEUX AUTRES POINTS
              N1=0
              N2=0
              IF ((MCN(NTEMP)  .NE.NUMST).AND.(N1.EQ.0)) N1=MCN(NTEMP)
              IF ((MCN(NTEMP+1).NE.NUMST).AND.(N1.EQ.0)) N1=MCN(NTEMP+1)
              IF ((MCN(NTEMP+2).NE.NUMST).AND.(N1.EQ.0)) N1=MCN(NTEMP+2)
C
              IF ((MCN(NTEMP).NE.NUMST).AND.(N1.NE.MCN(NTEMP)))
     %          N2=MCN(NTEMP)
              IF ((MCN(NTEMP+1).NE.NUMST).AND.(N1.NE.MCN(NTEMP+1)))
     %              N2=MCN(NTEMP+1)
              IF ((MCN(NTEMP+2).NE.NUMST).AND.(N1.NE.MCN(NTEMP+2)))
     %              N2=MCN(NTEMP+2)
C              LES DEUX AUTRES POINTS SONT N1 ET N2
C             ON RECUPERE LES COORDONNEES
              CN1(1)=RMCN(MNXYZ1+WYZSOM+3*(N1-1))
              CN1(2)=RMCN(MNXYZ1+WYZSOM+3*(N1-1)+1)
              CN2(1)=RMCN(MNXYZ1+WYZSOM+3*(N2-1))
              CN2(2)=RMCN(MNXYZ1+WYZSOM+3*(N2-1)+1)
              CNPT(1)=RMCN(MNXYZ1+WYZSOM+3*(NPT-1))
              CNPT(2)=RMCN(MNXYZ1+WYZSOM+3*(NPT-1)+1)
C             ON TESTE LE PRODUIT SCALAIRE
              SCAL=((CNPT(2)-CN1(2))*(CNPT(1)-CN2(1)))-
     %             ((CNPT(1)-CN1(1))*(CNPT(2)-CN2(2)))
              IF (SCAL.NE.0) THEN
                 IF (MCN(NTEMP)  .EQ.NUMST) MCN(NTEMP  )=NPT
                 IF (MCN(NTEMP+1).EQ.NUMST) MCN(NTEMP+1)=NPT
                 IF (MCN(NTEMP+2).EQ.NUMST) MCN(NTEMP+2)=NPT
                 IF (MCN(NTEMP+3).EQ.NUMST) MCN(NTEMP+3)=NPT
              ELSE
C                LE TRIANGLE EST A SUPPRIMER
                 MCN(MSUPTR+I-1)=1
                 NSUPTR=NSUPTR+1
              ENDIF
           ENDIF
 50     CONTINUE
C
      ENDIF
C     ON SORT DE LA BOUCLE QUI DETERMINE QUELS TRIANGLES SUPPRIMER
C
C     LES TRIANGLES A SUPPRIMER SONT DANS MSUPTR
C     IL Y EN A NSUPTR
C     ==========================================
C     ON RECOPIE LES ELEMENTS
      NUMEF=NBEF
C     ON BOUCLE DANS MSUPTR
      DO 61 J=1,NBQTAJ
         IF (MCN(MSUPTR+J-1).EQ.0) THEN
C           LE TRIANGLE N'EST PAS A SUPPRIMER
            GOTO 61
         ENDIF
C        LE TRIANGLE NUMTRS EST A SUPPRIMER
         NUMTRS=MCN(MNNTAB+J-1)
C        LE TRIANGLE NUMEF EST A SUPPRIMER??
C        EST-IL DANS MNNTAB ?
62       INDTRI=0
         DO 65 I=1,NBQTAJ
           NUMTR=MCN(MNNTAB+I-1)
           IF(NUMTR.EQ.NUMEF) THEN
              INDTRI=I
              GOTO 70
           ENDIF
65       CONTINUE
C        SI INDTRI=0 NUMEF EST PAS DANS MNNTAB
C        SINON POSITION INDTRI
C        SI DANS MNNTAB, A SUPPRIMER ??
70       IF (INDTRI.NE.0) THEN
           IF (MCN(MSUPTR+INDTRI-1).EQ.1) THEN
C            IL EST A SUPPRIMER
             NUMEF=NUMEF-1
             GOTO 62
           ENDIF
         ENDIF
C        ARRIVE ICI ON PEUT COPIER NUMEF DANS NUMTRS
C        SI NUMEF SUPERIEUR A NUMTRS
         IF (NUMEF.GT.NUMTRS) THEN
           NTEMP1=MNNSEF1+WUSOEF+4*(NUMEF-1)
           NTEMP2=MNNSEF1+WUSOEF+4*(NUMTRS-1)
           MCN(NTEMP2)=MCN(NTEMP1)
           MCN(NTEMP2+1)=MCN(NTEMP1+1)
           MCN(NTEMP2+2)=MCN(NTEMP1+2)
           MCN(NTEMP2+3)=MCN(NTEMP1+3)
         ENDIF
C        ET ON MET 0 DANS MSUPTR(J)
         MCN(MSUPTR+J-1)=0
         NUMEF=NUMEF-1
61     CONTINUE
C
C      RESTE A GERER LE TABLEAU NSEF
C      ON A ENLEVE NSUPTR ELEMENTS, CE SONT LES DERNIERS
       NBEF=MCN(MNNSEF1+WBEFOB)
C      ON MET LES ELEMENTS ENLEVES A ZERO
       DO 80 I=1,NSUPTR
         NTEMP=MNNSEF1+WUSOEF+4*(NBEF-I)
         MCN(NTEMP)=0
         MCN(NTEMP+1)=0
         MCN(NTEMP+2)=0
         MCN(NTEMP+3)=0
80     CONTINUE
       MCN(MNNSEF1+WBEFOB)=NBEF-NSUPTR
       MCN(MNNSEF1+WUTFMA)=0
       MCN(MNNSEF1+WBEFTG)=0
       MCN(MNNSEF1+WBEFAP)=0
       CALL ECDATE(MCN(MNNSEF1))
       MCN(MNNSEF1+MOTVAR(6))=NONMTD('~>>>NSEF')

C      RESTE A GERER LE TABLEAU XYZSOMMET
C      ON A SUPPRIME LE SOMMET NUMST
C      ON VA LE PERMUTER AVEC LE DERNIER
C      EST-CE LE DERNIER???
       IF( NBSOM .NE. NUMST ) THEN
C        PAS LE DERNIER
         NTEMP1=MNXYZ1+WYZSOM+3*(NUMST-1)
         NTEMP2=MNXYZ1+WYZSOM+3*(NBSOM-1)
         RMCN(NTEMP1  )=RMCN(NTEMP2)
         RMCN(NTEMP1+1)=RMCN(NTEMP2+1)
         RMCN(NTEMP1+2)=RMCN(NTEMP2+2)
C        ON BOUCLE DANS NSEF POUR REMPLACER NBSOM PAR NUMST
         DO 110 K=1,NBEF
           NTEMP=MNNSEF1+WUSOEF+4*(K-1)
           IF (MCN(NTEMP)  .EQ.NBSOM) MCN(NTEMP  )=NUMST
           IF (MCN(NTEMP+1).EQ.NBSOM) MCN(NTEMP+1)=NUMST
           IF (MCN(NTEMP+2).EQ.NBSOM) MCN(NTEMP+2)=NUMST
           IF (MCN(NTEMP+3).EQ.NBSOM) MCN(NTEMP+3)=NUMST
110         CONTINUE
       ENDIF
C
C      MAINTENANT ON SUPPRIME LE DERNIER
       NTEMP=MNXYZ1+WYZSOM+3*(NBSOM-1)
       RMCN(NTEMP  )=0
       RMCN(NTEMP+1)=0
       RMCN(NTEMP+2)=0
       MCN(MNXYZ1+WNBSOM)=NBSOM-1
       CALL ECDATE(MCN(MNXYZ1))
       MCN(MNXYZ1+MOTVAR(6))=NONMTD('~>>>XYZSOMMET')
C
 500   IF( MNNTAB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOM,  MNNTAB )
       IF( MSUPTR .GT. 0 ) CALL TNMCDS( 'ENTIER', LSUPTR, MSUPTR )
       RETURN
       END
