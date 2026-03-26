       SUBROUTINE ETOILED( NUMPT, MNNTAB, INDICE, MNXYZ, MNNSEF, NPT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER UN POINT DANS LE MAILLAGE A L'INTERIEUR D'UN TRIANGLE
C -----
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C MNXYZ1 : ADRESSE DU TABLEAU XYZSOMMET
C MNNSEF1: ADRESSE DU TABLEAU NSEF
C NBEF   : NOMBRE D'EF
C MNSOEL : ADRESSE DU TABLEAU DE TOUS LES SOMMETS DU MAILLAGE
C NX     : ABSCISSE DU POINT A RAJOUTER
C NY     : ORDONNEE DU POINT A RAJOUTER
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C NBSOM  : LE NOMBRE DE SOMMETS
C NBEF   : LE NOMBRE D'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : KCEHA & NICO    LJLL UPMC                                2000
C MODIFS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE   (RMCN(1),MCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER NB
      INTEGER TSOM(1:4)
      INTEGER TR1,TR2
      INTEGER SOM1
      INTEGER SOM2
      INTEGER SOM
      INTEGER SOMTMP
      REAL    CSOM1(1:2)
      REAL    CSOM2(1:2)
      REAL    CSOM(1:2)
      REAL    CSOMTMP(1:2)
      INTEGER TROUVE
      INTEGER S1,S2,S3
      INTEGER L
      REAL    SCTMP1,SCTMP2
      REAL    SCAL1
C
C     ON REPERE LES SOMMETS MIS EN JEU
C     ================================
      NB=0
      NBSOM=MCN(MNXYZ+WNBSOM)
      CALL TNMCDC( 'ENTIER', NBSOM, MTEMP )
      DO 10 I=1,NBSOM
         MCN(MTEMP+I-1)=0
 10   CONTINUE
      DO 20 I=1,INDICE
         NUMTR=MCN(MNNTAB+I-1)
         TSOM(1)=MCN(MNNSEF+WUSOEF+4*(NUMTR-1))
         TSOM(2)=MCN(MNNSEF+WUSOEF+4*(NUMTR-1)+1)
         TSOM(3)=MCN(MNNSEF+WUSOEF+4*(NUMTR-1)+2)
         DO 30 K=1,3
            IF (TSOM(K).NE.NUMPT) THEN
               IF (MCN(MTEMP+TSOM(K)-1).EQ.0) NB=NB+1
               MCN(MTEMP+TSOM(K)-1)=1
            ENDIF
 30      CONTINUE
 20   CONTINUE
C
C     ON VERIFIE MTEMP
C     ================
      DO 35 L=1,NBSOM
         PRINT *,'SOMMET ',L,' : ',MCN(MTEMP+L-1)
 35   CONTINUE
C
C     ON GENERE LE TABLEAU CONTENANT LES NUMEROS DES SOMMETS
C     ======================================================
      LMNSOM = NB
      IF( NB .LE. 0 ) LMNSOM = 1
      CALL TNMCDC( 'ENTIER', LMNSOM, MNSOM )
      INDTMP=0
      DO 40 I=1,NBSOM
         IF (MCN(MTEMP+I-1).EQ.1) THEN
            MCN(MNSOM+INDTMP)=I
            INDTMP=INDTMP+1
         ENDIF
 40   CONTINUE
C
C     ON DETRUIT LE TABLEAU TEMPORAIRE
C     ================================
      CALL TNMCDS('ENTIER',NBSOM,MTEMP)
C
C     VERIFICATION DU TABLEAU
C     =======================
      PRINT *,'NB : ',NB
      DO 45 L=1,NB
         PRINT *,'SOMMET ',L,' NUMERO : ',MCN(MNSOM+L-1)
 45   CONTINUE
C
C     ON COMMENCE L'ALGORITHME
C     ========================
C     ON BOUCLE SUR LES SOMMETS
      NPT=0
      DO 50 I=1,NB
         TR1=0
         TR2=0
         SOM=MCN(MNSOM+I-1)
C        BOUCLE SUR LES TRIANGLES POUR TROUVER 2 TRIANGLES DE SOMMET SOM
         DO 60 J=1,INDICE
            NUMTR=MCN(MNNTAB+J-1)
            TSOM(1)=MCN(MNNSEF+WUSOEF+4*(NUMTR-1))
            TSOM(2)=MCN(MNNSEF+WUSOEF+4*(NUMTR-1)+1)
            TSOM(3)=MCN(MNNSEF+WUSOEF+4*(NUMTR-1)+2)
            S1=TSOM(1)
            S2=TSOM(2)
            S3=TSOM(3)
            IF ((S1.EQ.SOM).OR.(S2.EQ.SOM).OR.(S3.EQ.SOM)) THEN
               IF (TR1.EQ.0) THEN
                  TR1=NUMTR
               ELSE
                  TR2=NUMTR
               ENDIF
            ENDIF
 60      CONTINUE
C
C        ON A REPERE LES DEUX TRIANGLES TR1 ET TR2 CONTENANT SOM
C        QUELS SONT LES DEUX AUTRES POINTS ?
         TSOM(1)=MCN(MNNSEF+WUSOEF+4*(TR1-1))
         TSOM(2)=MCN(MNNSEF+WUSOEF+4*(TR1-1)+1)
         TSOM(3)=MCN(MNNSEF+WUSOEF+4*(TR1-1)+2)
         IF ((TSOM(1).NE.NUMPT).AND.(TSOM(1).NE.SOM))SOM1=TSOM(1)
         IF ((TSOM(2).NE.NUMPT).AND.(TSOM(2).NE.SOM)) SOM1=TSOM(2)
         IF ((TSOM(3).NE.NUMPT).AND.(TSOM(3).NE.SOM)) SOM1=TSOM(3)
C        ON A RECUPERE LE PREMIER
         TSOM(1)=MCN(MNNSEF+WUSOEF+4*(TR2-1))
         TSOM(2)=MCN(MNNSEF+WUSOEF+4*(TR2-1)+1)
         TSOM(3)=MCN(MNNSEF+WUSOEF+4*(TR2-1)+2)
         IF ((TSOM(1).NE.NUMPT).AND.(TSOM(1).NE.SOM)) SOM2=TSOM(1)
         IF ((TSOM(2).NE.NUMPT).AND.(TSOM(2).NE.SOM)) SOM2=TSOM(2)
         IF ((TSOM(3).NE.NUMPT).AND.(TSOM(3).NE.SOM)) SOM2=TSOM(3)
C        ON A RECUPERE LE DEUXIEME
         PRINT *,'SOM : ',SOM
         PRINT *,'TR1 : ',TR1
         PRINT *,'TR2 : ',TR2
         PRINT *,'SOM1 : ',SOM1
         PRINT *,'SOM2 : ',SOM2
C
C        ON INTIALISE LES COORDONNEES
         CSOM1(1)=RMCN(MNXYZ+WYZSOM+3*(SOM1-1))
         CSOM1(2)=RMCN(MNXYZ+WYZSOM+3*(SOM1-1)+1)
         CSOM2(1)=RMCN(MNXYZ+WYZSOM+3*(SOM2-1))
         CSOM2(2)=RMCN(MNXYZ+WYZSOM+3*(SOM2-1)+1)
         CSOM(1)=RMCN(MNXYZ+WYZSOM+3*(SOM-1))
         CSOM(2)=RMCN(MNXYZ+WYZSOM+3*(SOM-1)+1)
         PRINT *,'CSOM1 : ',CSOM1(1),' ',CSOM1(2)
         PRINT *,'CSOM2 : ',CSOM2(1),' ',CSOM2(2)
         PRINT *,'CSOM : ',CSOM(1),' ',CSOM(2)
C      
C      
C        ON EST ETOILE PAR RAPPORT A SOM ?
         SCAL1 = (CSOM1(2)-CSOM(2)) * (CSOM2(1)-CSOM(1))
     %         - (CSOM1(1)-CSOM(1)) * (CSOM2(2)-CSOM(2))
         TROUVE=1
         DO 70 K=1,NB
            SOMTMP=MCN(MNSOM+K-1)
            PRINT *,'SOMTMP : ',SOMTMP
            CSOMTMP(1)=RMCN(MNXYZ+WYZSOM+3*(SOMTMP-1))
            CSOMTMP(2)=RMCN(MNXYZ+WYZSOM+3*(SOMTMP-1)+1)
            PRINT *,'SOMTMP : ',CSOMTMP(1),' ',CSOMTMP(2)
            SCTMP1=(CSOM1(2)-CSOM(2))*(CSOMTMP(1)-CSOM(1))-
     %             (CSOM1(1)-CSOM(1))*(CSOMTMP(2)-CSOM(2))
            SCTMP2=(CSOM2(2)-CSOM(2))*(CSOMTMP(1)-CSOM(1))-
     %             (CSOM2(1)-CSOM(1))*(CSOMTMP(2)-CSOM(2))
            PRINT *,'SCAL1 : ',SCAL1
            PRINT *,'SCTMP1 : ',SCTMP1
            PRINT *,'SCTMP2 : ',SCTMP2
            IF( (SCTMP1*SCTMP2.GE.0) .AND. (SOMTMP.NE.SOM) .AND.
     %          (SOMTMP.NE.SOM1) .AND. (SOMTMP.NE.SOM2) ) THEN
               TROUVE=0
               GOTO 80
            ENDIF
 70      CONTINUE
 80      PRINT *,'TROUVE : ',TROUVE
         IF (TROUVE.EQ.1) THEN
            NPT=SOM
            GOTO 90
         ENDIF
 50   CONTINUE
C     ON SORT DE LA BOUCLE SUR LES SOMMETS
C
90    CALL TNMCDS( 'ENTIER', LMNSOM, MNSOM )
      RETURN
      END
