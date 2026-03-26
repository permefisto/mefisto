      SUBROUTINE SUEX38( NTLXSU, LADEFI,
     %                   NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXTRAIRE UNE SURFACE D'UNE SURFACE A PARTIR DE LA DEFINITION
C -----    AVEC LA SOURIS D'UN RECTANGLE CONTENANT TOTALEMENT LES FACES
C       ATTENTION: EN 3D, LA NORMALE A LA FACE DOIT POINTER VERS L'OEIL
C                  POUR ETRE EXTRAITE
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Septembre 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      INTEGER           NOSOEL(12)
      LOGICAL           OUI
      CHARACTER*24      NMSUIN
      REAL              R, PT1(3), PT2(3), AXYZ(3), COORNO(3)
C
C     SI NON INTERACTIF => RETOUR
      IF( INTERA .LT. 3 ) THEN
         IERR = 1
         RETURN
      ENDIF
      MNFALV = 0
      MNOUIS = 0
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
      IF( NDIMLI .GE. 3 ) THEN
C        ASSURER L'ORIENTATION DU MAILLAGE DE LA SURFACE
C        PAR PARCOURS DES EF A TRAVERS LES ARETES COMMUNES
C        ET PERMUTATION DES SOMMETS 2-NBSTEF DU SECOND EF
C        D'UNE ORIENTATION DIFFERENTE DU PREMIER EF
         CALL SUORIENT( NMSUIN, NTFASF, MNFASF,  NTSOSF, MNSOSF, IERR )
      ENDIF
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
C     NON ou OUI PRISE EN COMPTE DE LA NORMALE POINTANT VERS L'OEIL EN 3D
      NORMON = LADEFI( WORMON )
C
C     LE NOMBRE DE FACES DE LA SURFACE EXTRAITE
      NBFASU = 0
      NBSOSU = 0
      MNFALV = 0
C
C     TRACE DE LA SURFACE ET DEFINITION DE LA DIRECTION DE VISEE ADEQUATE
      LORBITE = 1
      CALL TRAFAC( NMSUIN, NUSUIN, MNFASF, MNSOSF )
C
C     CLIQUER 2 FOIS POUR SELECTIONNER UN RECTANGLE SUR L'ECRAN
      CALL CLICRECT( NOTYEV, PT1, PT2 )
      LORBITE = 0
      IF( NOTYEV .LE. 0 ) GOTO 9999
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNFASF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF,
     %             NBEFOB, NX, NY, NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     CHAQUE SOMMET EST IL DANS LE RECTANGLE DE SELECTION?
      CALL TNMCDC( 'ENTIER', NBSOM, MNOUIS )
      MNO = MNOUIS - 1
      MNS = MNSOSF + WYZSOM - 3
      NBSTOK = 0
      DO N=1,NBSOM
C
C        LES 3 COORDONNEES DU SOMMET N
         MN = MNS + 3 * N
         IF( NDIMLI .EQ. 3 ) THEN
C           LES COORDONNEES AXONOMETRIQUES
            CALL XYZAXO( RMCN(MN), AXYZ )
         ELSE
            AXYZ(1) = RMCN(MN)
            AXYZ(2) = RMCN(MN+1)
            AXYZ(3) = 0
         ENDIF
C
C        LE POINT AXYZ EST IL DANS LE RECTANGLE?
         IF( PT1(1) .LE. AXYZ(1) .AND. AXYZ(1) .LE. PT2(1)  .AND.
     %       PT1(2) .LE. AXYZ(2) .AND. AXYZ(2) .LE. PT2(2) ) THEN
C           OUI:
            NBSTOK = NBSTOK + 1
            MCN( MNO+N ) = 1
         ELSE
C           NON:
            MCN( MNO+N ) = 0
         ENDIF
C
      ENDDO
C
      IF( NBSTOK .EQ. 0 ) THEN
C        AUCUN SOMMET N'EST SELECTIONNE => PAS DE FACES EXTRAITES
         IERR = 4
         GOTO 9990
      ENDIF
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
C     LE VECTEUR DIRECTION DE VISEE PTV -> OEIL
      IF( NDIMLI .GE. 3 ) THEN
         DIREVI(1) = AXOEIL(1) - AXOPTV(1)
         DIREVI(2) = AXOEIL(2) - AXOPTV(2)
         DIREVI(3) = AXOEIL(3) - AXOPTV(3)
      ENDIF
C
C     LA BOUCLE SUR LES EF DU MAILLAGE
C     --------------------------------
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
C           LE SOMMET EST IL DANS LE RECTANGLE?
            OUI = OUI .AND. (MCN(MNO+NOSOEL(I)).EQ. 1)
 20      CONTINUE
C
C        EN 3D: CALCUL DE LA NORMALE A LA FACE
         IF( NDIMLI .EQ. 3 .AND. NORMON .NE. 0 ) THEN
            CALL NORFA3( RMCN(MNS+3*NOSOEL(1)),
     %                   RMCN(MNS+3*NOSOEL(2)),
     %                   RMCN(MNS+3*NOSOEL(3)), COORNO, IERR )
C
C           LA NORMALE VA T ELLE VERS L'OEIL?
            R = PROSCR( COORNO, DIREVI, 3 )
            IF( R .LE. 0 ) GOTO 100
         ENDIF
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
 9990 CALL TNMCDS( 'ENTIER', NBSOM,  MNOUIS )
C
 9999 RETURN
      END
