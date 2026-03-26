      SUBROUTINE VOEX22( NTLXVO , LADEFI , RADEFI ,
     &                   NTCUVO , MNCUVO , NTSOCU , MNSOCU , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE MAILLAGE D'UN DEMI CYLINDRE
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
C           CF '~TD/D/A___SOMMETS'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR  : ALI MARS DEA D'A.N. MAI 1989
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
      IERR = 0
C
C     PARAMETRES DU MAILLAGE
C     ======================
      NBARCY = LADEFI(WBARCY) + 1
      NBAHCY = LADEFI(WBAHCY) + 1
      NBSECY = LADEFI(WBSECY)
      RAYOCS = RADEFI(WAYOCS)
      RAYOCI = RADEFI(WAYOCI)
      NUPOCY = LADEFI(WUPOCY)
      RAGRCY = RADEFI(WAGRCY)
      RAGHCY = RADEFI(WAGHCY)
      RAGACY = RADEFI(WAGACY)
C
C     VERIFICATION DES PARAMETRES
C     ===========================
      IF( NBARCY .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:5),'(I5)') NBARCY
         KERR(1) = 'NOMBRE DE POINTS SUR R INCORRECT ='
     %              // KERR(MXLGER)(1:5)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBAHCY .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:5),'(I5)') NBAHCY
         KERR(1) =  'NOMBRE DE POINTS SUR H INCORRECT ='
     %              // KERR(MXLGER)(1:5)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBSECY .LT. 3 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:5),'(I5)') NBSECY
         KERR(1) =  'NOMBRE DE SECTEURS INCORRECT ='
     %              // KERR(MXLGER)(1:5)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( RAYOCS .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAYOCS
         KERR(1) =  'RAYON SUPERIEUR INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( RAYOCI .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAYOCI
         KERR(1) =  'RAYON INFERIEUR INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGRCY.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGRCY
         KERR(1) =  'RAGRCY INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGHCY.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGHCY
         KERR(1) =  'RAGHCY INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGACY.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGACY
         KERR(1) =  'RAGACY INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     RECUPERATION DES 3 COORDONNEES DU CENTRE DU DEMI-CERCLE SUP
C     ===========================================================
      CALL LXNLOU( NTPOIN , NUPOCY , NTLXSC , MN )
      CALL LXTSOU( NTLXSC , 'XYZSOMMET' , NTSOSC , MNSOSC )
      MN = MNSOSC + WYZSOM
      X = RMCN( MN     )
      Y = RMCN( MN + 1 )
      Z = RMCN( MN + 2 )
C
      IF ( Z.EQ.0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SOMMET DU DEMI-CYLINDRE INCORRECT'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS
C     ======================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      NBH = NBAHCY
      NBR = NBARCY
      NBSEC = NBSECY
      NBSOM = NBH * ((NBR-1)*NBSEC + 1)
      CALL LXTNDC( NTLXVO , 'XYZSOMMET' , 'MOTS' , WYZSOM+3*NBSOM )
      CALL LXTSOU( NTLXVO , 'XYZSOMMET' , NTSOCU , MNSOCU )
C
C     GENERATION DES NSEF
C     ==========================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
      NBCUVO = (NBR-1)*(NBSEC-1)*(NBH-1)
      CALL LXTNDC(NTLXVO,'NSEF','ENTIER',WUSOEF+8*NBCUVO)
      CALL LXTSOU(NTLXVO,'NSEF',NTCUVO,MNCUVO)
C
C     GENERATION DES SOMMETS ET ELEMENTS DU DEMI-CYLINDRE
C     ==================================================
      CALL DEMICY( NBSEC , NBR , NBH , X , Y , Z , RAYOCI ,
     &             RAYOCS , RAGHCY , RAGRCY , RAGACY ,
     &             RMCN(MNSOCU+WYZSOM) , MCN(MNCUVO+WUSOEF) )
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
      SUBROUTINE DEMICY( NBSEC , NBR , NBH , X , Y , Z , RI , RS ,
     &                   QH , QR ,QS ,
     &                   XYZSOM , NOEUDS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT    : GENERER LE MAILLAGE D'UN DEMI-CYLINDRE
C
C ENTREES :
C --------
C NBH    : NOMBRE DE POINTS SUR LA HAUTEUR
C NBR    : NOMBRE DE POINTS SUR LE RAYON
C NBSEC  : NOMBRE DE SECTIONS ANGULAIRES
C X      : COORDONNEE X DU SOMMET DU CENTRE DU DEMICERCLE SUPERIEUR
C Y      : COORDONNEE Y DU SOMMET DU CENTRE DU DEMICERCLE SUPERIEUR
C Z      : COORDONNEE Z DU SOMMET DU CENTRE DU DEMICERCLE SUPERIEUR
C RI     : RAYON DU DEMICERCLE INFERIEUR
C RS     : RAYON DU DEMICERCLE SUPERIEUR
C QR     : RAISON GEOMETRIQUE DE CONCENTRATION SUR LE RAYON
C QH     : RAISON GEOMETRIQUE DE CONCENTRATION SUR LA HAUTEUR
C QS     : RAISON GEOMETRIQUE DE CONCENTRATION SUR L'ANGLE
C
C CONTRAINTES SUR LES ENTREES :
C ----------------------------
C NBH >= 2
C NBR >= 2
C NBSEC >= 3
C (X,Y,Z) DIFFERENT DE (0,0,0) ET Z DIFFERENT DE 0
C RS > 0
C RI > 0
C QR > 0 ( 1 - MAILLAGE HOMOGENE )
C QH > 0 ( 1 - MAILLAGE HOMOGENE )
C QS > 0 ( 1 - MAILLAGE HOMOGENE )
C
C SORTIES :
C --------
C XYZSOM : TABLEAU DES COORDONNEES DES NOEUDS
C NOEUDS : TABLEAU D'AFFECTATION DES SOMMETS AUX ELEMENTS
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALI MARS DEA D'A.N. MAI 1989
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION XYZSOM(3,*)
      DIMENSION NOEUDS(8,*)
      REAL X,Y,Z
      REAL RI,RS
C
CCC      NBS = NBH * ((NBR-1)*NBSEC + 1)
CCC      NBELF = (NBSEC-1)*(NBR-1)*(NBH-1)
      NBSN = 1+NBSEC*(NBR-1)
      NBELFN = (NBR-1)*(NBSEC-1)
C
      PI = 3.1415926535897932385
C
      DO K=1,NBH
         NUMS = (K-1)*NBSN+1
         XYZSOM(1,NUMS) = DIVISN( QH , K , NBH ) * X
         XYZSOM(2,NUMS) = DIVISN( QH , K , NBH ) * Y
         XYZSOM(3,NUMS) = DIVISN( QH , K , NBH ) * Z
      ENDDO

      DO K=1,NBH
         DO J=2,NBR
            DO I=1,NBSEC
               NUMS=(K-1)*NBSN+(J-2)*NBSEC+1+I
               RK = RI+DIVISN(QH,K,NBH)*(RS-RI)
               VJ = DIVISN(QR,J,NBR)
               WI = DIVISN(QS,I,NBSEC)
               XYZSOM(1,NUMS)=DIVISN(QH,K,NBH)*X+RK*VJ*COS(PI*WI)
               XYZSOM(2,NUMS)=DIVISN(QH,K,NBH)*Y+RK*VJ*SIN(PI*WI)
               XYZSOM(3,NUMS)=DIVISN(QH,K,NBH)*Z
            ENDDO
         ENDDO
      ENDDO

      DO K=1,NBH-1
         DO I=1,NBSEC-1
            NUMEF = (K-1)*NBELFN+I
            NOEUDS(1,NUMEF) = (K-1)*NBSN+1
            NOEUDS(2,NUMEF) = (K-1)*NBSN+I+1
            NOEUDS(3,NUMEF) = (K-1)*NBSN+2+I
            NOEUDS(4,NUMEF) = K*NBSN+1
            NOEUDS(5,NUMEF) = K*NBSN+1+I
            NOEUDS(6,NUMEF) = K*NBSN+2+I
            NOEUDS(7,NUMEF) = 0
            NOEUDS(8,NUMEF) = 0
            DO J=2,NBR-1
               NUMEF = (K-1)*NBELFN+(J-1)*(NBSEC-1)+I
               NOEUDS(1,NUMEF) = (K-1)*NBSN+1+I+(J-2)*NBSEC
               NOEUDS(2,NUMEF) = (K-1)*NBSN+I+1+(J-2)*NBSEC+NBSEC
               NOEUDS(3,NUMEF) = (K-1)*NBSN+2+I+NBSEC+(J-2)*NBSEC
               NOEUDS(4,NUMEF) = (K-1)*NBSN+2+I+(J-2)*NBSEC
               NOEUDS(5,NUMEF) = K*NBSN+1+I+(J-2)*NBSEC
               NOEUDS(6,NUMEF) = K*NBSN+1+I+NBSEC+(J-2)*NBSEC
               NOEUDS(7,NUMEF) = K*NBSN+2+I+NBSEC+(J-2)*NBSEC
               NOEUDS(8,NUMEF) = K*NBSN+2+I+(J-2)*NBSEC
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
