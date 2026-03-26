      SUBROUTINE SUEX41( NTLXSU, LADEFI, RADEFI,
     %                   NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXTRAIRE UN TRIANGLE OU QUADRANGLE D'UNE SURFACE SI EN
C -----    TOUS CES SOMMETS (NORMALE,VECTEUR)>=0
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA SURFACE
C          CF ~/td/d/a_surface__definition
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DE LA SURFACE
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DE LA SURFACE
C          CF ~/td/d/a___nsef
C NTSOSU : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          1 SI SURFACE INITIALE INCONNUE
C          2 SI SURFACE INITIALE SANS NSEF
C          3 SI SURFACE INITIALE SANS XYZSOMMET
C          4 SI AUCUNE FACE EXTRAITE
C          5 SI VECTEUR DIRECTION EST NUL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Septembre 2010
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*), VECDIR(1:3)
      INTEGER           NOSOEL(12)
      CHARACTER*24      NMSUIN
      REAL              COORNO(3)
C
      IERR   = 0
      MNFALV = 0
C
C     LA SURFACE INITIALE
C     ===================
C     LE NOM DE CETTE SURFACE
      NUSUIN = LADEFI( WUSUIN )
C     LE TABLEAU LEXIQUE DE CETTE SURFACE
      CALL LXNLOU( NTSURF, NUSUIN, NTLXSF, MN )
      IF( NTLXSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE INITIALE INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN INITIAL SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE NOM DE LA SURFACE INITIALE
      CALL NMOBNU( 'SURFACE', NUSUIN, NMSUIN )
C     ASSURER L'ORIENTATION DU MAILLAGE DE LA SURFACE
C     PAR PARCOURS DES EF A TRAVERS LES ARETES COMMUNES
C     ET PERMUTATION DES SOMMETS 2-NBSTEF DU SECOND EF
C     D'UNE ORIENTATION DIFFERENTE DU PREMIER EF
      CALL SUORIENT( NMSUIN, NTFASF, MNFASF,  NTSOSF, MNSOSF, IERR )
C
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'NSEF', NTFASF, MNFASF )
      IF( NTFASF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE SANS TMS NSEF'
         ELSE
            KERR(1) = 'SURFACE WITHOUT NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTSOSF, MNSOSF )
      IF( NTSOSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE SANS TMS XYZSOMMET'
         ELSE
            KERR(1) = 'SURFACE WITHOUT XYZSOMMET TMS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
      NBSOM = MCN( MNSOSF + WNBSOM )
      NBTGS = MCN( MNSOSF + WNBTGS )
C
C     LES 3 COMPOSANTES CARTESIENNES DU VECTEUR
      VECDIR(1) = RADEFI(WECDIR)
      VECDIR(2) = RADEFI(WECDIR+1)
      VECDIR(3) = RADEFI(WECDIR+2)
      IF( VECDIR(1)**2 + VECDIR(2)**2 + VECDIR(3)**2 .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VECTEUR NUL'
         ELSE
            KERR(1) = 'NULL VECTOR'
         ENDIF
         CALL LEREUR
         IERR = 5
         RETURN
      ENDIF
C
C     LE NOMBRE DE FACES DE LA SURFACE EXTRAITE
      NBFASU = 0
      NBSOSU = 0
      MNFALV = 0
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNFASF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF,
     %             NBEFOB, NX, NY, NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU TEMPORAIRE
C     DE STOCKAGE DES NBSOEF SOMMETS DE CHAQUE EF EXTRAIT
C                 DES NBTGEF TANGENTES ET DU CODE GEOMETRIQUE
      MXFALV = NBEFOB * ( NBSOEF + NBTGEF + 1 )
      CALL TNMCDC( 'ENTIER', MXFALV, MNFALV )
      MNAVS = MNFALV - 1
      MNAVT = MNAVS + NBEFOB * NBSOEF
      MNAVG = MNAVT + NBEFOB * NBTGEF
C
C     LA BOUCLE SUR LES EF DU MAILLAGE
C     --------------------------------
      MNS = MNSOSF + WYZSOM - 3
      DO 100 N=1,NBEFOB
C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNFASF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C        LE NOMBRE DE SOMMETS DE CET ELEMENT = NCOGEL
C
         IF( NCOGEL .EQ. 3 ) THEN
C
C           CALCUL DU VECTEUR NORMAL AU TRIANGLE
            CALL NORFA3( RMCN(MNS+3*NOSOEL(1)),
     %                   RMCN(MNS+3*NOSOEL(2)),
     %                   RMCN(MNS+3*NOSOEL(3)), COORNO, IERR )
C
C           LE PRODUIT SCALAIRE (VECTEUR NORMAL, VECTEUR DIRECTION)
            PRSCNV = PROSCR( COORNO, VECDIR, 3 )
            IF( PRSCNV .LT. 0 ) GOTO 100
C
         ELSE
C
C           CALCUL DU VECTEUR NORMAL EN CHAQUE SOMMET DU QUADRANGLE
            I0 = 4
            DO 20 I=1,NCOGEL
C
C              CALCUL DU VECTEUR NORMAL AU SOMMET I DU QUADRANGLE
               IF( I .EQ. 4 ) THEN
                  I1 = 1
               ELSE
                  I1 = I+1
               ENDIF
               CALL NORFA3( RMCN(MNS+3*NOSOEL(I )),
     %                      RMCN(MNS+3*NOSOEL(I1)),
     %                      RMCN(MNS+3*NOSOEL(I0)), COORNO, IERR )
C
C              LE PRODUIT SCALAIRE (VECTEUR NORMAL, VECTEUR DIRECTION)
               PRSCNV = PROSCR( COORNO, VECDIR, 3 )
               IF( PRSCNV .LT. 0 ) GOTO 100
C
C              PASSAGE AU SOMMET SUIVANT
               I0 = I
C
 20         CONTINUE
C
         ENDIF
C
C        L'EF VERIFIE LE CRITERE
         NBFASU = NBFASU + 1
C        LES NUMEROS DES SOMMETS SONT STOCKES
         DO 30 I=1,NCOGEL
            MCN( MNAVS + I ) = NOSOEL( I )
 30      CONTINUE
C        COMPLETION EVENTUELLE AVEC DES ZEROS
         DO 40 I=NCOGEL+1,NBSOEF
            MCN( MNAVS + I ) = 0
 40      CONTINUE
         MNAVS = MNAVS + NBSOEF
C
C        LES NUMEROS DES NBTGEF TANGENTES OU 0 TANGENTE
         IF( NBTGEF .EQ. 0 .OR. NUEFTG .EQ. 0 ) THEN
C           LES NBTGEF TANGENTES SONT NULLES
            DO 50 I=1,NBTGEF
               MCN( MNAVT + I ) = 0
 50         CONTINUE
C           CODE GEOMETRIQUE
            MCN( MNAVG + 1 ) = 0
         ELSE
            DO 60 I=1,NBTGEF
               MCN( MNAVT + I ) = NOSOEL( NBSOEF + I )
 60         CONTINUE
C           CODE GEOMETRIQUE
            MCN( MNAVG + 1 ) = NUGEEF
         ENDIF
         MNAVT = MNAVT + NBTGEF
         MNAVG = MNAVG + 1
C
 100  CONTINUE
C
C     RENUMEROTATION DES SOMMETS ET TANGENTES DES EF EXTRAITS
      CALL REEFEX( 3,      NTLXSU, NBSOEF, NBTGEF,
     %             NBSOM,  NBTGS,  MNSOSF,
     %             NBFASU, MNFALV, MNFALV+NBEFOB*NBSOEF,
     %             MNFALV+NBEFOB*(NBSOEF+NBTGEF),
     %             NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
      CALL TNMCDS( 'ENTIER', MXFALV, MNFALV )
C
 9999 RETURN
      END
