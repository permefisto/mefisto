      SUBROUTINE VOEX24( NTLXVO , LADEFI , RADEFI ,
     &                   NTCUVO , MNCUVO , NTSOCU , MNSOCU , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE MAILLAGE D'UN CONE OU D'UN DEMI CONE
C ----
C
C ENTREES :
C --------
C NTLXVO  : NUMERO DU TABLEAU TS DU LEXIQUE DU CONE
C LADEFI  : TABLEAU ENTIER DE DEFINITION DU VOLUME PARTITIONNEE
C RADEFI  : TABLEAU REEL   DE DEFINITION DU VOLUME PARTITIONNEE
C           CF '~TD/D/A_VOLUME__DEFINITION'
C
C SORTIES :
C ---------
C NTCUVO  : NUMERO      DU TMS 'NSEF' DES NUMEROS DES EF
C MNCUVO  : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES EF
C           CF '~TD/D/A___NSEF'
C NTSOCU  : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOCU  : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C           CF '~TD/D/A___XYZSOMMET'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR  : MICHEL BULIK DEA D'A.N. MAI 1989
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
      IERR = 0
C
C     PARAMETRES DU MAILLAGE
C     ======================
      NBINRC = LADEFI(WBINRC)
      NBINHC = LADEFI(WBINHC)
      NBINSC = LADEFI(WBINSC)
      RAYONC = RADEFI(WAYONC)
      NUPOSC = LADEFI(WUPOSC)
      IF( LADEFI(WUTYVO) .EQ. 23 ) THEN
C        CONE
         NCHOIX = 2
      ELSE
C        DEMI CONE
         NCHOIX = 1
      ENDIF
      RAGERA = RADEFI(WAGERA)
      RAGEHA = RADEFI(WAGEHA)
      RAGEAN = RADEFI(WAGEAN)
C
C     VERIFICATION DES PARAMETRES
C     ===========================
      IF( NBINRC .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBINRC
         KERR(1) =  'NOMBRE D''ARETES EN R INCORRECT ='
     %             //KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBINHC .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBINHC
         KERR(1) =  'NOMBRE D''ARETES EN H INCORRECT ='
     %             //KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBINSC .LT. 1+NCHOIX ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBINSC
         KERR(1) =  'NOMBRE DE SECTEURS INCORRECT ='
     %             //KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( RAYONC .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAYONC
         KERR(1) =  'RAYON INCORRECT ='
     %             //KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBINHC .LT. NBINRC ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE D''ARETES EN R > NOMBRE D''ARETES EN H !'
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGERA.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGERA
         KERR(1) = 'RAGERA INCORRECT ='
     %             //KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGEHA.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGEHA
         KERR(1) =  'RAGEHA INCORRECT ='
     %             //KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGEAN.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGEAN
         KERR(1) =  'RAGEAN INCORRECT ='
     %             //KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     RECUPERATION DES 3 COORDONNEES DU SOMMET DU CONE
C     ================================================
      CALL LXNLOU( NTPOIN , NUPOSC , NTLXSC , MN )
      CALL LXTSOU( NTLXSC , 'XYZSOMMET' , NTSOSC , MNSOSC )
      MN = MNSOSC + WYZSOM
      X = RMCN( MN     )
      Y = RMCN( MN + 1 )
      Z = RMCN( MN + 2 )
C
      IF ( Z.EQ.0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'SOMMET DU CONE INCORRECT'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS
C     ======================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      NH = NBINHC + 1
      NR = NBINRC + 1
      NC = NBINSC
      NPAS = INTDIV(NH-1,NR-2)
      NBSOM = NH + (NC+2-NCHOIX) * (1 + 2 * ( NH-2-NPAS*(NR-3) ) +
     &                      NPAS * (NR+2) * (NR-3) / 2 )
      CALL LXTNDC( NTLXVO , 'XYZSOMMET' , 'MOTS' , WYZSOM+3*NBSOM )
      CALL LXTSOU( NTLXVO , 'XYZSOMMET' , NTSOCU , MNSOCU )
C
C     GENERATION DES NSEF
C     ==========================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
      NBCUVO = NC * ( NPAS * (NR-3) * (NR-2) / 2 + 2 * NH - 3)
      CALL LXTNDC(NTLXVO,'NSEF','ENTIER',WUSOEF+8*NBCUVO)
      CALL LXTSOU(NTLXVO,'NSEF',NTCUVO,MNCUVO)
C
C     GENERATION DES SOMMETS ET ELEMENTS DU CONE
C     ==========================================
      CALL CONE3( NH , NR , NC , X , Y , Z , RAYONC ,
     &            NCHOIX , RAGERA , RAGEHA , RAGEAN ,
     &            RMCN(MNSOCU+WYZSOM) , MCN(MNCUVO+WUSOEF) , IERR )
C
C     MISE A JOUR DES TABLEAUX
C     ========================
C     MISE A JOUR DU TABLEAU 'XYZSOMMET' DE CE VOLUME
C     NBSOM 'nombre de sommets'
      MCN( MNSOCU + WNBSOM ) = NBSOM
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOCU) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOCU + WNBTGS ) = 0
      MCN( MNSOCU + WBCOOR ) = 3
      MCN( MNSOCU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE VOLUME
C     TYPE DE L'OBJET : VOLUME
      MCN( MNCUVO + WUTYOB ) = 4
C
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNCUVO + WUTFMA ) = -1
C
C     PAS DE TANGENTES STOCKEES
      MCN( MNCUVO + WBTGEF ) = 0
      MCN( MNCUVO + WBEFAP ) = 0
      MCN( MNCUVO + WBEFTG ) = 0
C
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN( MNCUVO + WUTYMA ) = 0
C
C     NBSOEF 'nombre de sommets par sous-objets'
      MCN( MNCUVO + WBSOEF ) = 8
C
C     NBEFOB 'nombre de sous-objets de l''objet'
      MCN( MNCUVO + WBEFOB ) = NBCUVO
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNCUVO) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCUVO + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
      END
C234567--------------------------------------------------------------012
      SUBROUTINE CONE3(NH , NR , NC , X , Y , Z , R ,
     &                 NCHOIX , RAGERA ,RAGEHA , RAGEAN ,
     &                 XYZSOM , NOEUDS , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT    : GENERER LE MAILLAGE D'UN CONE
C
C ENTREES :
C --------
C NH     : NOMBRE DE POINTS SUR LA HAUTEUR
C NR     : NOMBRE DE POINTS SUR LE RAYON DE LA BASE DU CONE
C NC     : NOMRE DE SECTIONS ANGULAIRES
C X      : COORDONNEE X DU SOMMET DU CONE
C Y      : COORDONNEE Y DU SOMMET DU CONE
C Z      : COORDONNEE Z DU SOMMET DU CONE
C R      : RAYON DE LA BASE DU CONE
C NCHOIX : CHOIX ENTRE CONE ET DEMICONE (1-DEMICONE,2-CONE)
C RAGERA : RAISON GEOMETRIQUE DE CONCENTRATION SUR LE RAYON
C RAGEHA : RAISON GEOMETRIQUE DE CONCENTRATION SUR LA HAUTEUR
C RAGEAN : RAISON GEOMETRIQUE DE CONCENTRATION SUR L'ANGLE
C
C CONTRAINTES SUR LES ENTREES :
C ----------------------------
C NH >= 3
C NR >= 3
C NH >= NR
C NC >= 3 POUR LE CONE ET >= 2 POUR LE DEMICONE
C (X,Y,Z) DIFFERENT DE (0,0,0) ET Z DIFFERENT DE 0
C R > 0
C NCHOIX = 1  OU  NCHOIX = 2
C RAGERA > 0 ( 1 - MAILLAGE HOMOGENE )
C RAGEHA > 0 ( 1 - MAILLAGE HOMOGENE )
C RAGEAN > 0 ( 1 - MAILLAGE HOMOGENE ) CA MARCHE QUAND NC = PAIR
C
C SORTIES :
C --------
C XYZSOM : TABLEAU DES COORDONNEES DES NOEUDS
C NOEUDS : TABLEAU D'AFFECTATION DES SOMMETS AUX ELEMENTS
C IERR   : 0 SI PAS D'ERREUR
C          > 0 SINON
C
C VARIABLES UTILISEES :
C --------------------
C NPAS   : VOIR VOEX24.DOC POUR EXPLICATION
C NEFACT : NOMBRE DES ELEMENTS SUR LE RAYON AU NIVEAU I
C NUMEF  : NUMERO D'ELEMENT
C NBEF   : NOMBRE TOTAL DES ELEMENTS
C NBS    : NOMBRE TOTAL DES SOMMETS
C I      : NIVEAU ( VARIE DE 1 A NH )
C J      : DISTANCE DE L'AXE ( VARIE DE 1=AXE A NRACT )
C K      : ANGLE ( VARIE DE 1 A NC+1 )
C NUMSOM : NUMERO DU SOMMET
C NRACT  : NOMBRE DE SOMMETS SUR LE RAYON AU NIVEAU I
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : MICHEL BULIK DEA D'A.N. MAI 1989
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/gsmenu.inc"
      INTEGER NPAS,NR,NH,NC,NEFACT,NUMEF,NBEF,NBS,I,J,K,NUMSOM,NRACT
      INTEGER NOEUDS(8,*),NCHOIX
      REAL XYZSOM(3,*)
      REAL R,X,Y,Z,RI,XI,YI,ZI,ALFA,PI,RAGERA,RAGEHA,RAGEAN
C
C     CALCUL DU NOMBRE DES SOMMETS ET DES ELEMENTS
C
      NPAS = INTDIV(NH-1,NR-2)
      NBEF = NC*(NPAS*(NR-3)*(NR-2)/2 + 2*NH - 3)
      NBS = NH + (NC+2-NCHOIX)*(1 + 2*(NH-2-NPAS*(NR-3)) +
     &                   NPAS*(NR+2)*(NR-3)/2 )
C
C     INITIALISATION DES VARIABLES
C     'NUMSOM' EST LE NUMERO DU SOMMET
C     'NRACT' EST LE NOMBRE DE POINTS SUR LE RAYON AU NIVEAU I
C
      NUMSOM = 0
      NRACT = NR+1
      PI    = ATAN(1.0) * 4
      IERR  = 0
C
C     BOUCLE SUR LA PARTIE REGULIERE DU MAILLAGE
      DO I=1,NPAS*(NR-3)
         RI = DIVISP(R,0.,RAGEHA,1,NH,I)
         XI = DIVISP(0.,X,RAGEHA,1,NH,I)
         YI = DIVISP(0.,Y,RAGEHA,1,NH,I)
         ZI = DIVISP(0.,Z,RAGEHA,1,NH,I)
         NUMSOM = NUMSOM+1
         XYZSOM(1,NUMSOM) = XI
         XYZSOM(2,NUMSOM) = YI
         XYZSOM(3,NUMSOM) = ZI
         IF ((NPAS .EQ. 1) .OR. (INTMOD(I,NPAS) .EQ. 1)) THEN
            NRACT = NRACT-1
         ENDIF
         DO J=2,NRACT
            DO K=1,NC+2-NCHOIX
               IF ((INTMOD(NC,2).EQ.0).AND.(RAGEAN.NE.1)) THEN
                  IF (K.LE.(NC/2+1)) THEN
                     ALFA = DIVISP(0.,NCHOIX*PI/2.0,RAGEAN,1,
     %                             INTDIV(NC,2)+1,K)
                     GOTO 90
                  ENDIF
                  ALFA = DIVISP(NCHOIX*PI/2.0,NCHOIX*PI,1/RAGEAN,
     %                          INTDIV(NC,2)+1,NC+1,K)
                  GOTO 90
               ENDIF
               ALFA = DIVISP(0.,NCHOIX*PI,1.0,1,NC+1,K)

 90            NUMSOM = NUMSOM+1
               XYZSOM(1,NUMSOM) = XI+DIVISP(0.,RI,RAGERA,1,NRACT,J)
     %                              *COS(ALFA)
               XYZSOM(2,NUMSOM) = YI+DIVISP(0.,RI,RAGERA,1,NRACT,J)
     %                              *SIN(ALFA)
               XYZSOM(3,NUMSOM) = ZI
            ENDDO
         ENDDO
      ENDDO
C
C     BOUCLE SUR LA PARTIE OU NRACT = 3
      NRACT = NRACT-1
      DO I=1+NPAS*(NR-3),NH-2
         RI = DIVISP(R,0.,RAGEHA,1,NH,I)
         XI = DIVISP(0.,X,RAGEHA,1,NH,I)
         YI = DIVISP(0.,Y,RAGEHA,1,NH,I)
         ZI = DIVISP(0.,Z,RAGEHA,1,NH,I)
         NUMSOM = NUMSOM+1
         XYZSOM(1,NUMSOM) = XI
         XYZSOM(2,NUMSOM) = YI
         XYZSOM(3,NUMSOM) = ZI
         DO J=2,NRACT
            DO K=1,NC+2-NCHOIX
               IF ((INTMOD(NC,2).EQ.0).AND.(RAGEAN.NE.1)) THEN
                  IF (K.LE.(NC/2+1)) THEN
                     ALFA = DIVISP( 0.,NCHOIX*PI/2.0,RAGEAN,1,
     %                              INTDIV(NC,2)+1,K )
                     GOTO 190
                  ENDIF
                  ALFA = DIVISP( NCHOIX*PI/2.0,NCHOIX*PI,1/RAGEAN,
     &                           INTDIV(NC,2)+1,NC+1,K )
                  GOTO 190
               ENDIF
               ALFA = DIVISP(0.,NCHOIX*PI,1.0,1,NC+1,K)
 190           NUMSOM = NUMSOM+1
               XYZSOM(1,NUMSOM) = XI+DIVISP(0.,RI,RAGERA,1,NRACT,J)
     %                              *COS(ALFA)
               XYZSOM(2,NUMSOM) = YI+DIVISP(0.,RI,RAGERA,1,NRACT,J)
     %                              *SIN(ALFA)
               XYZSOM(3,NUMSOM) = ZI
            ENDDO
         ENDDO
      ENDDO
C
C     CALCUL AU NIVEAU  I = NH-1 ET NRACT = 2
      I = NH-1
      RI = DIVISP(R,0.,RAGEHA,1,NH,I)
      XI = DIVISP(0.,X,RAGEHA,1,NH,I)
      YI = DIVISP(0.,Y,RAGEHA,1,NH,I)
      ZI = DIVISP(0.,Z,RAGEHA,1,NH,I)
      NUMSOM = NUMSOM+1
      XYZSOM(1,NUMSOM) = XI
      XYZSOM(2,NUMSOM) = YI
      XYZSOM(3,NUMSOM) = ZI
      DO K=1,NC+2-NCHOIX
         IF ((INTMOD(NC,2).EQ.0).AND.(RAGEAN.NE.1)) THEN
            IF (K.LE.(NC/2+1)) THEN
               ALFA = DIVISP(0.,NCHOIX*PI/2.0,RAGEAN,1,INTDIV(NC,2)+1,K)
               GOTO 290
            ENDIF
            ALFA = DIVISP( NCHOIX*PI/2.0,NCHOIX*PI,1/RAGEAN,
     &                     INTDIV(NC,2)+1,NC+1,K )
            GOTO 290
         ENDIF
         ALFA = DIVISP(0.,NCHOIX*PI,1.0,1,NC+1,K)
 290     NUMSOM = NUMSOM+1
         XYZSOM(1,NUMSOM) = XI+RI*COS(ALFA)
         XYZSOM(2,NUMSOM) = YI+RI*SIN(ALFA)
         XYZSOM(3,NUMSOM) = ZI
      ENDDO
C
C     CALCUL AU NIVEAU I = NH I.E. LE SOMMET DU CONE
C
      NUMSOM = NUMSOM+1
      XYZSOM(1,NUMSOM) = X
      XYZSOM(2,NUMSOM) = Y
      XYZSOM(3,NUMSOM) = Z
C
      IF ( NUMSOM.NE.NBS ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VOEX24:NUMSOM<>NBS'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     FIN DE CALCUL DES COORDONNEES DES NOEUDS
C     DEBUT D'AFFECTATION DES NOEUDS AUX ELEMENTS
C     INITIALISATION DES VARIABLES
C     'NEFACT' EST LE NOMBRE D'ELEMENTS SUR LE RAYON AU NIVEAU I
C     'NUMEF' EST LE NUMERO D'ELEMENT
C
      NEFACT = NR
      NUMEF = 0
C
C     BOUCLE SUR LA PARTIE REGULIERE DU CONE
C
      DO 105 I=1,(NR-3)*NPAS
      J=1
      DO K=1,NC
C        PENTAEDRES AUTOUR DE L'AXE
         NUMEF = NUMEF+1
         NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,1)
         NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
         NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
         NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,1)
         NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K)
         NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K+1)
         NOEUDS(7,NUMEF)=0
         NOEUDS(8,NUMEF)=0
      ENDDO
      IF ((NPAS .EQ. 1) .OR. (INTMOD(I,NPAS) .EQ. 1)) THEN
      NEFACT = NEFACT-1
      ENDIF
      DO J=2,NEFACT-1
         DO K=1,NC
C           HEXAEDRES A L'INTERIEUR DU CONE
            NUMEF = NUMEF+1
            NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K)
            NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
            NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
            NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K+1)
            NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K)
            NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K)
            NOEUDS(7,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K+1)
            NOEUDS(8,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K+1)
         ENDDO
      ENDDO
      J = NEFACT
      DO 130 K=1,NC
         IF (INTMOD(I,NPAS).EQ.0) THEN
C           A LA FIN IL Y A OU BIEN UN PENTAEDRE HORIZONTAL
            NUMEF = NUMEF+1
            NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K)
            NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
            NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K)
            NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K+1)
            NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
            NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K+1)
            NOEUDS(7,NUMEF)=0
            NOEUDS(8,NUMEF)=0
            GOTO 130
         ENDIF

C        OU BIEN UN HEXAEDRE
         NUMEF = NUMEF+1
         NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K)
         NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
         NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
         NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K+1)
         NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K)
         NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K)
         NOEUDS(7,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K+1)
         NOEUDS(8,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K+1)
 130  ENDDO
 105  ENDDO
C
C     BOUCLE SUR LA PARTIE OU NEFACT=2
      DO I=(NR-3)*NPAS+1,NH-3
         J=1
         DO K=1,NC
C     PENTAEDRES AUTOUR DE L'AXE
            NUMEF = NUMEF+1
            NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,1)
            NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
            NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
            NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,1)
            NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K)
            NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K+1)
            NOEUDS(7,NUMEF)=0
            NOEUDS(8,NUMEF)=0
         ENDDO
         J=2
         DO K=1,NC
C           HEXAEDRES
            NUMEF = NUMEF+1
            NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K)
            NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
            NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
            NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K+1)
            NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K)
            NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K)
            NOEUDS(7,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K+1)
            NOEUDS(8,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K+1)
         ENDDO
      ENDDO
C
C     NIVEAU AVANT DERNIER I.E. I=NH-2
C
      I=NH-2
      J=1
      DO K=1,NC
C        PENTAEDRES AUTOUR DE L'AXE
         NUMEF = NUMEF+1
         NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,1)
         NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
         NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
         NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,1)
         NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K)
         NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J+1,K+1)
         NOEUDS(7,NUMEF)=0
         NOEUDS(8,NUMEF)=0
      ENDDO
      J=2
      DO K=1,NC
C        PENTAEDRES HORIZONTAUX
         NUMEF = NUMEF+1
         NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K)
         NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
         NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K)
         NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K+1)
         NOEUDS(6,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
         NOEUDS(5,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K+1)
         NOEUDS(7,NUMEF)=0
         NOEUDS(8,NUMEF)=0
      ENDDO
C
C     DERNIER NIVEAU I.E. I=NH-1
C
      I=NH-1
      J=1
      DO K=1,NC

C        TETRAEDRES AU SOMMET DU CONE
         NUMEF = NUMEF+1
         NOEUDS(1,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J,K)
         NOEUDS(2,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K)
         NOEUDS(3,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I,J+1,K+1)
         NOEUDS(4,NUMEF)=NUMS(NCHOIX,NH,NR,NC,I+1,J,K)
         NOEUDS(5,NUMEF)=0
         NOEUDS(6,NUMEF)=0
         NOEUDS(7,NUMEF)=0
         NOEUDS(8,NUMEF)=0

      ENDDO
C
      IF ( NUMEF.NE.NBEF ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VOEX24: NUMEF <> NBEF !'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

      RETURN
      END
C234567--------------------------------------------------------------012
      INTEGER FUNCTION INTDIV(I,J)
      INTEGER K
      K = 0
 100  IF (K*J .LE. I) THEN
         K = K+1
         GOTO 100
      ENDIF
      INTDIV = K-1
      RETURN
      END
C234567--------------------------------------------------------------012
      INTEGER FUNCTION INTMOD(I,J)
      INTEGER K
      K = 0
 100  IF (K*J .LE. I) THEN
      K = K+1
      GOTO 100
      ENDIF
      INTMOD = I - (K-1)*J
      RETURN
      END
C234567--------------------------------------------------------------012
      REAL FUNCTION DIVISP(DEBUT,FIN,GRAD,NDEB,NFIN,N)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA DIVISION D'UN INTERVALLE
C
C ENTREES :
C --------
C DEBUT  : DEBUT D'INTERVALLE A DIVISER
C FIN    : FIN D'INTERVALLE A DIVISER
C GRAD   : RAISON GEOMETRIQUE DE CONCENTRATION
C          0 - 1       DEBUT -> FIN
C          1 - INFINI  FIN -> DEBUT
C NDEB   : NUMERO DU DEBUT D'INTERVALLE
C NFIN   : NUMERO DE LA FIN D'INTERVALLE
C N      : NUMERO DE LA DIVISION DESIREE
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : MICHEL BULIK DEA D'A.N. MAI 1989
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (NDEB.GT.NFIN) THEN
      STOP 'ERREUR DANS DIVISP : NDEB > NFIN !'
      ELSE IF ((N .GT. NFIN) .OR. (N .LT. NDEB)) THEN
      STOP 'ERREUR DANS DIVISP : N INCORRECT !'
      ELSE IF (GRAD.LE.0) THEN
      STOP 'ERREUR DANS DIVISP : RAISON GEOMETRIQUE INCORRECT !'
      ENDIF
      IF ( GRAD .EQ. 1.0 ) THEN
      DIVISP = DEBUT + (FIN-DEBUT)*REAL(N-NDEB)/REAL(NFIN-NDEB)
      GOTO 100
      ENDIF
      DIVISP = DEBUT+(FIN-DEBUT)*(1-GRAD**(N-NDEB))/
     &                                (1-GRAD**(NFIN-NDEB))
 100  RETURN
      END
C234567--------------------------------------------------------------012
      INTEGER FUNCTION NUMS(NCHOIX,NH,NR,NC,I,J,K)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DONNER LE NUMERO DU NOEUD QUI SE TROUVE
C        AU NIVEAU I
C        A LA DISTANCE J DE L'AXE DU CONE ( 1 = AXE )
C        A L'ANGLE K
C
C ENTREES :
C --------
C NCHOIX : CHOIX ENTRE CONE ET DEMI CONE (1-DEMICONE,2-CONE)
C NH     : NOMBRE DE POINTS SUR LA HAUTEUR DU CONE
C NR     : NOMBRE DE POINTS SUR LE RAYON DE LA BASE DU CONE
C NC     : NOMBRE DE SECTIONS ANGULAIRES
C I      : NIVEAU DU SOMMET ( 1=BASE DU CONE )
C J      : DISTANCE DU SOMMET DE L'AXE ( 1=AXE )
C K      : ANGLE DU SOMMET
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : MICHEL BULIK DEA D'A.N. MAI 1989
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER NCHOIX,NH,NR,NC,I,J,K,NPAS,NRTEMP,ITEMP
C
C     VERIFICATION DES COORDONNEES
C
      IF ((I.LT.1).OR.(I.GT.NH)) THEN
      STOP 'ERREUR DANS NUMS : NIVEAU INCORRECT !'
      ELSE IF ((J.LT.1).OR.(J.GT.NR)) THEN
      STOP 'ERREUR DANS NUMS : RAYON INCORRECT !'
      ELSE IF ((K.LT.1).OR.(K.GT.NC+1)) THEN
      STOP 'ERREUR DANS NUMS : ANGLE INCORRECT !'
      ENDIF
C
C     INITIALISATION DES VARIABLES
C
      NPAS = INTDIV(NH-1,NR-2)
      NUMS = 0
      NRTEMP = NR+1
C
C     DEBUT DE CALCUL
C     PARTIE REGULIERE, AUDESSOUS DE I
C
      DO ITEMP=1,INTMIN(NPAS*(NR-3),I-1)
         IF ((NPAS.EQ.1).OR.(INTMOD(ITEMP,NPAS).EQ.1)) THEN
            NRTEMP = NRTEMP-1
         ENDIF
         NUMS = NUMS+1+(NC+2-NCHOIX)*(NRTEMP-1)
      ENDDO
      IF (J.GT.NRTEMP) THEN
      STOP 'ERREUR # 1 DANS NUMS : J > NRTEMP !'
      ENDIF
C
C     SI I EST DANS LA PARTIE REGULIERE
C
      IF (I.LE.NPAS*(NR-3)) THEN
      NUMS = NUMS+1
      IF (J.GE.2) THEN
         NUMS = NUMS+(1+INTMOD(K-1,NC))*(NCHOIX-1)+
     &         (NC+2-NCHOIX)*(J-2)+K*(2-NCHOIX)
      ENDIF
      GOTO 1000
      ENDIF
C
C     PARTIE OU NRTEMP = 3, AUDESSOUS DE I
C
      NRTEMP = NRTEMP-1
      IF (J.GT.NRTEMP) THEN
      STOP 'ERREUR # 2 DANS NUMS : J > NRTEMP !'
      ENDIF
      DO ITEMP=(NR-3)*NPAS+1,INTMIN(I-1,NH-2)
         NUMS = NUMS+1+(NC+2-NCHOIX)*(NRTEMP-1)
      ENDDO
C
C     SI I SE TROUVE DANS LA PARTIE OU NRTEMP=3
C
      IF (I.LE.NH-2) THEN
      NUMS = NUMS+1
      IF (J.GE.2) THEN
      NUMS = NUMS+(1+INTMOD(K-1,NC))*(NCHOIX-1)+
     & (NC+2-NCHOIX)*(J-2)+K*(2-NCHOIX)
      ENDIF
      GOTO 1000
      ENDIF
C
C     PARTIE OU NRTEMP=2 I.E. NIVEAU NH-1
C
      NRTEMP = NRTEMP-1
      IF (J.GT.NRTEMP) THEN
      STOP 'ERREUR # 3 DANS NUMS : J > NRTEMP !'
      ENDIF
      IF (I.EQ.NH-1) THEN
      NUMS = NUMS+1
      IF (J.EQ.2) THEN
      NUMS = NUMS+(1+INTMOD(K-1,NC))*(NCHOIX-1)+K*(2-NCHOIX)
      ENDIF
      GOTO 1000
      ENDIF
C
C     SI I=NH I.E. LE SOMMET DU CONE
      NRTEMP = NRTEMP-1
      IF (J.GT.NRTEMP) THEN
      STOP 'ERREUR # 4 DANS NUMS : J > NRTEMP !'
      ENDIF
      NUMS = NUMS+2+(NC+2-NCHOIX)
C
C     TOUS LE GOTO SONT DIRIGES ICI
C
 1000 RETURN
      END
C234567--------------------------------------------------------------012
      INTEGER FUNCTION INTMIN(I,J)
      IF (I.LT.J) THEN
      INTMIN = I
      GOTO 10
      ENDIF
      INTMIN = J
 10   RETURN
      END
