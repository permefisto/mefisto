       SUBROUTINE SUPPRP3D( MNXYZ1, MNNSEF1,  NUMST, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER UN SOMMET DU MAILLAGE D'UNE SURFACE 3D en
C -----    SUPPRIMANT TOUS LES TRIANGLES QUADRANGLES DONT IL EST SOMMET

C ENTREES:
C --------
C MNXYZ1 : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF1: ADRESSE DU TABLEAU NSEF DE LA SURFACE

C SORTIES:
C --------
C NUMST  : >0 NUMERO DU SOMMET SUPPRIME
C          =0 SI PAS DE SOMMET CLIQUE
C IERR   : =0 SI PAS D'ERREUR
C           1 PAS DE POINT CLIQUE RETROUVE
C           2 PRESENCE DE QUADRANGLE
C           9 SATURATION DU SUPER-TABLEAU MCN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C MODIFS : PERRONNET ALAIN Saint PIERRE du PERRAY               Mai 2020
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
      INTEGER       NUMST, N1, N2
      REAL          CN1(1:2)
      REAL          CN2(1:2)
      REAL          CNPT(1:2)

      IERR   = 0
      MNNTAB = 0
      MNSUTQ = 0

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZ1+WNBSOM)

C     NBEFOB NOMBRE DES EF
      NBEFOB = MCN(MNNSEF1+WBEFOB)

C     NBSOEF NOMBRE DE SOMMETS PAR EF
      NBSOEF = MCN(MNNSEF1+WBSOEF)

C     ADRESSE MCN DE DEBUT DES NO DE SOMMETS DES EF
      MNSOEL = MNNSEF1 + WUSOEF

C     SELECTIONNER A L'AIDE D'UN CLIC SOURIS UN SOMMET DU MAILLAGE
C     PRESENT DANS LA FENETRE GRAPHIQUE C-A-D
C     RETOURNER LE NUMERO NUMST0 DU SOMMET VISIBLE LE PLUS PROCHE
C     DU POINT CLIQUE EN NX, NY DANS UN MAILLAGE 3D
C     ------------------------------------------------------------
      CALL SESTCLIC( RMCN(MNXYZ1+WYZSOM), ITST, NUMST )
      IF( NUMST .LE. 0 ) GOTO 9990

C     RECHERCHE DES QUADRANGLES OU TRIANGLES CONTENANT LE SOMMET NUMST
C     ----------------------------------------------------------------
      CALL TNMCDC( 'ENTIER', NBSOM, MNNTAB )
      IF( MNNTAB .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'supprp3d: PAS ASSEZ de PLACE en MCN'
         ELSE
            KERR(1) = 'supprp3d: NOT ENOUGH MEMORY in MCN'
         ENDIF
         CALL LEREUR
         IERR = 9
         RETURN
      ENDIF

      NBTQAS=0
      DO J=1,NBEFOB
         NTEMP = MNSOEL+4*(J-1)
         IF( (MCN(NTEMP)  .EQ.NUMST).OR.(MCN(NTEMP+1).EQ.NUMST)
     %   .OR.(MCN(NTEMP+2).EQ.NUMST).OR.(MCN(NTEMP+3).EQ.NUMST)) THEN
C            NO J DU TQ CONTENANT LE SOMMET NUMST
            MCN( MNNTAB+NBTQAS ) = J
            NBTQAS = NBTQAS + 1
         ENDIF
      ENDDO
      IF( NBTQAS .LE. 0 ) GOTO 9900

C     ON INITIALISE LE TABLEAU DES TRIANGLES ou QUADRANGLES A SUPPRIMER
      LSUPTR = NBTQAS
      CALL TNMCDC( 'ENTIER', LSUPTR, MNSUTQ )
      NSUPTR=0
      DO I=1,NBTQAS
         MCN(MNSUTQ+I-1)=0
      ENDDO

C     ON DETERMINE SI NUMST EST AU BORD DE LA TRIANGULATION DE LA SURFACE
      CALL PFRONT( MNNTAB, NBTQAS, MNNSEF1, MNXYZ1, NUMST, NFR )
C     NFR=1 SI A LA FRONTIERE, 0 SI INTERNE, -1 SI QUADRANGLE
      IF( NFR .EQ. -1 ) THEN

C        IL EXISTE UN QUADRANGLE
         IERR = 2
         GOTO 9900

      ELSE IF( NFR .EQ. 1 ) THEN

C        LE SOMMET NUMST EST A LA FRONTIERE DE LA TRIANGULATION
C        DONC TOUS LES TRIANGLES DOIVENT ETRE SUPPRIMES
         DO I=1,NBTQAS
            MCN(MNSUTQ+I-1)=1
         ENDDO
         NSUPTR=NBTQAS

      ELSE

C        LE SOMMET NUMST N'EST PAS A LA FRONTIERE DE LA TRIANGULATION
         NSUPTR=2

C        ON CHOISIT UN SOMMET NPT DIFFERENT DE NUMST
C        IL FAUDRAIT TROUVER NPT TEL QUE LA FIGURE
C        SOIT ETOILEE PAR RAPPORT A NPT
         CALL ETOILED( NUMST, MNNTAB, NBTQAS, MNXYZ1, MNNSEF1, NPT )
         IF ( NPT .EQ. 0 ) GOTO 9900

C        BOUCLE SUR LES TRIANGLES
         NN=0
         DO I=1,NBTQAS

C          NUMERO DU TRIANGLE A SUPPRIMER
           NUMTR=MCN(MNNTAB+I-1)
           NTEMP=MNNSEF1+WUSOEF+4*(NUMTR-1)
           IF (  (MCN(NTEMP  ).EQ.NPT)
     %       .OR.(MCN(NTEMP+1).EQ.NPT)
     %       .OR.(MCN(NTEMP+2).EQ.NPT)
     %       .OR.(MCN(NTEMP+3).EQ.NPT)) THEN
              IF (NN.EQ.0) THEN
C                IL FAUDRA SUPPRIMER NN1
                 NN1=NUMTR
                 MCN(MNSUTQ+I-1)=1
                 NN=1
              ELSE
C                IL FAUDRA SUPPRIMER NN2
                 NN2=NUMTR
                 MCN(MNSUTQ+I-1)=1
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
     %          N2=MCN(NTEMP+1)
              IF ((MCN(NTEMP+2).NE.NUMST).AND.(N1.NE.MCN(NTEMP+2)))
     %          N2=MCN(NTEMP+2)
C              LES DEUX AUTRES POINTS SONT N1 ET N2
C             ON RECUPERE LES COORDONNEES
              MNX = MNXYZ1+WYZSOM-3
              CN1(1) =RMCN(MNX+3*N1  )
              CN1(2) =RMCN(MNX+3*N1+1)
              CN2(1) =RMCN(MNX+3*N2  )
              CN2(2) =RMCN(MNX+3*N2+1)
              CNPT(1)=RMCN(MNX+3*NPT  )
              CNPT(2)=RMCN(MNX+3*NPT+1)
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
                 MCN(MNSUTQ+I-1)=1
                 NSUPTR=NSUPTR+1
              ENDIF
           ENDIF
        ENDDO

      ENDIF
C     ON SORT DE LA BOUCLE QUI DETERMINE QUELS TRIANGLES SUPPRIMER

C     LES TRIANGLES A SUPPRIMER SONT DANS MNSUTQ: IL Y EN A NSUPTR
C     ============================================================
      CALL SUPEFMAR( NSUPTR, MCN(MNSUTQ), NBSOEF, NBEFOB, MCN(MNSOEL) )

 9900 IF( MNNTAB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOM,  MNNTAB )
      IF( MNSUTQ .GT. 0 ) CALL TNMCDS( 'ENTIER', LSUPTR, MNSUTQ )

      IF(LANGAG .EQ. 0 ) THEN
         PRINT*,'supprp3d: le SOMMET',NUMST,' est SUPPRIME'
      ELSE
         PRINT*,'supprp3d: the VERTEX',NUMST,' is DELETED'
      ENDIF
      GOTO 9999


C     PAS DE POINT CLIQUE et PAS DE SOMMET SUPPRIME
C     ---------------------------------------------
 9990 NUMST = 0
      IERR  = 1
      IF(LANGAG .EQ. 0 ) THEN
         PRINT*,'supprp3d: PAS de SOMMET SUPPRIME'
      ELSE
         PRINT*,'supprp3d: NO VERTEX DELETED'
      ENDIF

 9999 RETURN
      END
