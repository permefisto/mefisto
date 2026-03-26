        SUBROUTINE VOEX61( NTLXVO, LADEFI, RADEFI,
     %                     NTEF6C, MNEF6C, NTST6C, MNST6C, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN 6CUBES CENTRE DE COTES
C -----    DE NARETE DE TAILLE D'ARETE CONSTANTE
C
C ENTREES:
C --------
C NTLXVO : NUMERO DU TABLEAU TS DU LEXIQUE DU 6-CUBE
C LADEFI : TABLEAU DE DEFINITION DES VOLUMES
C          CF '~td/d/a_volume__definition'
C
C SORTIES:
C --------
C NTEF6C : NUMERO      DU TMS 'NSEF' DES NUMEROS DES CUBES
C MNEF6C : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES CUBES
C          CF '~td/d/a___nsef'
C NTST6C : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME 6-CUBE
C MNST6C : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME 6-CUBE
C          CF '~td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  TEXAS A & M UNIVERSITY           JULY  2005
C2345X7..............................................................012
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
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      REAL              COIN(6,2)
C
      IERR = 0
C
C     DIMENSION DE L'ESPACE DU 6-CUBES
C     ================================
      NBCOOR = 6
C
C     MIN A L'ORIGINE OU CUBE6 CENTRE?
C     NC6CUM =0 : origine au minimum des coordonnees du cube6
C            =1 : origine centree au centre du cube6
C            =2 : origine centree au centre du sous cube6 milieu
C     ==========================================================
      NC6CUM = LADEFI(WUBE6M)
C
C     NOMBRE D'ARETES DANS LES 6 DIMENSIONS
C     =====================================
      NA6CUB = LADEFI(WUBE6N)
      WRITE(IMPRIM,*)'NOMBRE D''ARETES DANS UNE DIRECTION=',NA6CUB
      IF( NA6CUB .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NA6CUB
         KERR(1)='NOMBRE INCORRECT (<2) D''ARETES ='//KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     TAILLE CONSTANTE DES ARETES DANS LES 6 DIMENSIONS
C     =================================================
      AR6CUB = RADEFI(WUBE6A)
      IF( AR6CUB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:14),'(E14.6)') AR6CUB
         KERR(1)='TAILLE INCORRECTE D''ARETE ='//KERR(MXLGER)(1:14)
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      WRITE(IMPRIM,*)'LONGUEUR D''UNE ARETE DES 6-CUBES =',AR6CUB
      WRITE(IMPRIM,*)'LARGEUR DU 6-CUBE =',AR6CUB*NA6CUB
C
C     GENERATION DES SOMMETS DE CE 6-CUBES CENTRE
C     ===========================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      NBSOM = ( NA6CUB + 1 ) ** NBCOOR
      CALL LXTNDC( NTLXVO, 'XYZSOMMET', 'MOTS' , WYZSOM+NBCOOR*NBSOM )
      CALL LXTSOU( NTLXVO, 'XYZSOMMET',  NTST6C, MNST6C )
C     LA DIMENSION DE L'ESPACE = NOMBRE DES COORDONNEES D'UN SOMMET
      MCN( MNST6C + WBCOOR ) = NBCOOR
C     LE NOMBRE DE SOMMETS
      MCN( MNST6C + WNBSOM ) = NBSOM
C
C     CALCUL DES NBCOOR COORDONNEES DES NBSOM SOMMETS
C
      IF( NC6CUM .EQ. 0 ) THEN
C
C        MIN DES COORDONNEES DU CUBE 6 A l'ORIGINE
         X0 = 0.0
C
      ELSE
C
C        CUBE6 CENTRE AU MIEUX
C        L'ORIGINE DU REPERE NE DOIT PAS ETRE SUR UNE FACE DU 6-CUBES
C        CECI EST NECESSAIRE POUR L'ATOME D'HELIUM POUR LAQUELLE
C        LES POIDS ET COORDONNEES DE LA FORMULE D INTEGRATION EMPLOIENT
C        NPI=64  POINTS SUR LE 6CUBES XYZUVW OBTENUE PAR PRODUIT TENSORIEL DE
C        CELLE SUR L HEXAEDRE UNITE XYZ EXACTE POUR Q1 A 8 POINTS=SOMMETS ET
C        CELLE SUR L HEXAEDRE UNITE UVW EXACTE POUR Q3 A 8 POINTS DE GAUSS
         X0 = - ( AR6CUB * NA6CUB ) / 2
         IF( NC6CUM .EQ. 2 .AND. MOD(NA6CUB, 2) .EQ. 0 ) THEN
            X0 = X0 - AR6CUB / 2
         ENDIF
      ENDIF
C
      MN = MNST6C + WYZSOM
      DO 60 N=0,NA6CUB
         DO 50 M=0,NA6CUB
            DO 40 L=0,NA6CUB
               DO 30 K=0,NA6CUB
                  DO 20 J=0,NA6CUB
                     DO 10 I=0,NA6CUB
                        RMCN( MN     ) = X0 + I * AR6CUB
                        RMCN( MN + 1 ) = X0 + J * AR6CUB
                        RMCN( MN + 2 ) = X0 + K * AR6CUB
                        RMCN( MN + 3 ) = X0 + L * AR6CUB
                        RMCN( MN + 4 ) = X0 + M * AR6CUB
                        RMCN( MN + 5 ) = X0 + N * AR6CUB
                        MN  = MN + NBCOOR
 10                  CONTINUE
 20               CONTINUE
 30            CONTINUE
 40         CONTINUE
 50      CONTINUE
 60   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNST6C) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNST6C + WNBTGS ) = 0
      MCN( MNST6C + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     GENERATION DES NSEF
C     ===================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE 6-CUBE
      CALL LXTNDC( NTLXVO, 'NSEF', 'ENTIER', WDNCUB+1)
      CALL LXTSOU( NTLXVO, 'NSEF',  NTEF6C,  MNEF6C )
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE VOLUME
C     TYPE DE L'OBJET : VOLUME BIEN QU'EN FAIT CE SOIT UN 6-CUBE!
      MCN( MNEF6C + WUTYOB ) = 4
C     LE TYPE I6CONNU DE FERMETURE DU MAILLAGE
      MCN( MNEF6C + WUTFMA ) = -1
C     NUMERO DU TYPE DU MAILLAGE : 6CUBES STRUCTURE
      MCN( MNEF6C + WUTYMA ) = 8
C     NA6CUB NOMBRE D'ARETES DANS UNE DIRECTION
      MCN( MNEF6C + WANCUB ) = NA6CUB
C     ND6CUB DIMENSION N DU 6-CUBE
      MCN( MNEF6C + WDNCUB ) = NBCOOR
C     LE NOMBRE DE SOMMETS PAR 6-CUBE
      MCN( MNEF6C + WBSOEF ) = 2 ** NBCOOR
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : VOLUME C0
      MCN( MNEF6C + WBTGEF ) = 0
      MCN( MNEF6C + WBEFAP ) = 0
      MCN( MNEF6C + WBEFTG ) = 0
C     LE NOMBRE D'EF DU VOLUME
      MCN( MNEF6C + WBEFOB ) = NA6CUB ** NBCOOR
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNEF6C) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNEF6C + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     MIN ET MAX DES COORDONNEES
      CALL CADEXT( MNST6C, COIN )
C
C     RENUMEROTATION DES SOMMETS POUR REDUIRE LE PROFIL DE LA MATRICE
ccc      CALL RE6CUB( NTLXVO, NTEF6C, MNEF6C, MNST6C )
C     AUGMENTE LEGEREMENT LE PROFIL DE LA MATRICE! => SUPPRESSION
C
      RETURN
      END
