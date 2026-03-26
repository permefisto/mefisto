      SUBROUTINE SUEX36( NTLXSU, LADEFI,
     %                   NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EXTRAIRE UNE SURFACE D'UNE SURFACE A PARTIR D'UN CRITERE LOGIQUE
C -----
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF ~/TD/D/A_SURFACE__DEFINITION
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DE LA SURFACE
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DE LA SURFACE
C          CF ~/TD/D/A___NSEF
C NTSOSU : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~/TD/D/A___XYZSOMMET
C IERR   : 0 SI PAS D'ERREUR
C          1 SI SURFACE INITIALE SANS NSEF
C          2 SI SURFACE INITIALE SANS SOMMETS
C          3 SI FONCTION INCONNUE
C          4 SI AUCUNE FACE EXTRAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1996
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      INTEGER           NOSOEL(12)
      DOUBLE PRECISION  DXYZ(3),DOUI
      LOGICAL           OUI
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
         IERR = 1
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
         IERR = 2
         RETURN
      ENDIF
      NBSOM = MCN( MNSOSF + WNBSOM )
      NBTGS = MCN( MNSOSF + WNBTGS )
C
C     LE NUMERO DE LA FONCTION
      NUFOCS = LADEFI( WUFOCS )
      IF( NUFOCS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FONCTION CRITERE INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN CRITERION FUNCTION'
         ENDIF
         CALL LEREUR
         IERR = 3
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
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU CRITERE
      CALL TNMCDC( 'ENTIER', NBSOM, MNOUIS )
      MNO = MNOUIS - 1
C
C     CALCUL DU CRITERE EN CHACUN DES SOMMETS DU MAILLAGE
      MNS = MNSOSF + WYZSOM - 3
      DO 10 N=1,NBSOM
C
C        LES 3 COORDONNEES EN DOUBLE PRECISION DU SOMMET N
         MN = MNS + 3 * N
         DXYZ(1) = RMCN( MN   )
         DXYZ(2) = RMCN( MN+1 )
         DXYZ(3) = RMCN( MN+2 )
C
C        LE SOMMET N VERIFIE T IL LE CRITERE ?
         CALL FONVAL( NUFOCS, 3, DXYZ, NCODEV, DOUI )
         IF( NCODEV .EQ. 0 ) DOUI = 0
C        1 SI CRITERE VERIFIE, 0 SINON
         MCN( MNO+N ) = NINT( DOUI )
 10   CONTINUE
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
C     LA BOUCLE SUR LES NO SOMMET DU MAILLAGE
C     -----------------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      DO 100 N=1,NBEFOB
C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNFASF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C        LE NOMBRE DE SOMMETS DE CET ELEMENT = NCOGEL
C
         OUI = .TRUE.
         DO 20 I=1,NCOGEL
C           LE SOMMET VERIFIE T IL LE CRITERE ?
            OUI = OUI .AND. (MCN(MNO+NOSOEL(I)).EQ. 1)
 20      CONTINUE
C
         IF( OUI ) THEN
C           L'EF VERIFIE LE CRITERE
            NBFASU = NBFASU + 1
C           LES NUMEROS DES SOMMETS SONT STOCKES
            DO 30 I=1,NCOGEL
               MCN( MNAVS + I ) = NOSOEL( I )
 30         CONTINUE
C           COMPLETION EVENTUELLE AVEC DES ZEROS
            DO 40 I=NCOGEL+1,NBSOEF
               MCN( MNAVS + I ) = 0
 40         CONTINUE
            MNAVS = MNAVS + NBSOEF
C
C           LES NUMEROS DES NBTGEF TANGENTES OU 0 TANGENTE
            IF( NBTGEF .EQ. 0 .OR. NUEFTG .EQ. 0 ) THEN
C              LES NBTGEF TANGENTES SONT NULLES
               DO 50 I=1,NBTGEF
                  MCN( MNAVT + I ) = 0
 50            CONTINUE
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = 0
            ELSE
               DO 60 I=1,NBTGEF
                  MCN( MNAVT + I ) = NOSOEL( NBSOEF + I )
 60            CONTINUE
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = NUGEEF
            ENDIF
            MNAVT = MNAVT + NBTGEF
            MNAVG = MNAVG + 1
         ENDIF
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
      CALL TNMCDS( 'ENTIER', NBSOM,  MNOUIS )
 9999 RETURN
      END
