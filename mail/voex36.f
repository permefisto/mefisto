      SUBROUTINE VOEX36( NTLXVO , LADEFI ,
     %                   NTCUVO , MNCUVO , NTSOVO , MNSOVO , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EXTRAIRE UN VOLUME D'UN VOLUME A PARTIR D'UN CRITERE LOGIQUE
C -----
C
C ENTREES:
C --------
C NTLXVO : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME
C LADEFI : TABLEAU ENTIER DE DEFINITION DU VOLUME
C          CF ~/TD/D/A_VOLUME__DEFINITION
C
C SORTIES:
C --------
C NTCUVO : NUMERO      DU TMS 'NSEF' DU VOLUME
C MNCUVO : ADRESSE MCN DU TMS 'NSEF' DU VOLUME
C          CF ~/TD/D/A___NSEF
C NTSOVO : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOVO : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF ~/TD/D/A___XYZSOMMET
C IERR   : 0 SI PAS D'ERREUR
C          1 SI VOLUME INITIAL SANS NSEF
C          2 SI VOLUME INITIAL SANS XYZSOMMET
C          3 SI FONCTION INCONNUE
C          4 SI AUCUN CUBE EXTRAIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1996
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_volume__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      INTEGER           NOSOEL(64)
      DOUBLE PRECISION  DXYZ(3),DOUI
      LOGICAL           OUI
C
C     LE VOLUME INITIAL
C     ===================
C     LE NOM DE CE VOLUME
      NUVOIN = LADEFI( WUVOIN )
C     LE TABLEAU LEXIQUE DE CE VOLUME
      CALL LXNLOU( NTVOLU , NUVOIN , NTLXVL , MN )
      IF( NTLXVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME INITIAL INCONNU'
         ELSE
            KERR(1) = 'UNKNOWN INITIAL VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTSOU( NTLXVL , 'NSEF' , NTCUVL , MNCUVL )
      IF( NTCUVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS TMS NSEF'
         ELSE
            KERR(1) = 'VOLUME WITHOUT TMS NSEF'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CE VOLUME
      CALL LXTSOU( NTLXVL , 'XYZSOMMET' , NTSOVL , MNSOVL )
      IF( NTSOVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS TMS XYZSOMMET'
         ELSE
            KERR(1) = 'VOLUME WITHOUT TMS XYZSOMMET'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      NBSOM = MCN( MNSOVL + WNBSOM )
      NBTGS = MCN( MNSOVL + WNBTGS )
C
C     LE NUMERO DE LA FONCTION
      NUFOCV = LADEFI( WUFOCV )
      IF( NUFOCV .LE. 0 ) THEN
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
C     LE NOMBRE DE CUBES DU VOLUME EXTRAIT
      NBCUVO = 0
      NBSOVO = 0
      MNFALV = 0
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNCUVL) ,
     %             NUTYMA , NBSOEL , NBSOEF , NBTGEF,
     %             LDAPEF , LDNGEF , LDTGEF , NBEFOB ,
     %             NX     , NY     , NZ     ,
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU CRITERE
      CALL TNMCDC( 'ENTIER' , NBSOM , MNOUIS )
      MNO = MNOUIS - 1
C
C     CALCUL DU CRITERE EN CHACUN DES SOMMETS DU MAILLAGE
      MNS = MNSOVL + WYZSOM - 3
      DO 10 N=1,NBSOM
C
C        LES 3 COORDONNEES EN DOUBLE PRECISION DU SOMMET N
         MN = MNS + 3 * N
         DXYZ(1) = RMCN( MN   )
         DXYZ(2) = RMCN( MN+1 )
         DXYZ(3) = RMCN( MN+2 )
C
C        LE SOMMET N VERIFIE T IL LE CRITERE ?
         CALL FONVAL( NUFOCV , 3 , DXYZ , NCODEV , DOUI )
         IF( NCODEV .EQ. 0 ) DOUI = 0
C        1 SI CRITERE VERIFIE , 0 SINON
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
C     ----------------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      DO 100 N=1,NBEFOB
C        LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
         CALL NSEFNS( N      , NUTYMA , NBSOEF , NBTGEF,
     %                LDAPEF , LDNGEF , LDTGEF ,
     %                MNCUVL , NX , NY , NZ ,
     %                NCOGEL , NUGEEF , NUEFTG, NOSOEL , IERR )
C        LE NOMBRE DE SOMMETS DE CET ELEMENT
         NBSO = NBSOME( NCOGEL )
C
         OUI  = .TRUE.
         DO 20 I=1,NBSO
C           LE SOMMET VERIFIE T IL LE CRITERE ?
            OUI = OUI .AND. (MCN(MNO+NOSOEL(I)).EQ. 1)
 20      CONTINUE
C
         IF( OUI ) THEN
C           L'EF VERIFIE LE CRITERE
            NBCUVO = NBCUVO + 1
C           LES NUMEROS DES SOMMETS SONT STOCKES
            DO 30 I=1,NBSO
               MCN( MNAVS + I ) = NOSOEL( I )
 30         CONTINUE
C           COMPLETION EVENTUELLE AVEC DES ZEROS
            DO 40 I=NBSO+1,NBSOEF
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
      CALL REEFEX( 4,      NTLXVO, NBSOEF, NBTGEF,
     %             NBSOM,  NBTGS,  MNSOVL,
     %             NBCUVO, MNFALV, MNFALV+NBEFOB*NBSOEF,
     %             MNFALV+NBEFOB*(NBSOEF+NBTGEF),
     %             NTCUVO, MNCUVO, NTSOVO, MNSOVO, IERR )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
      CALL TNMCDS( 'ENTIER', MXFALV, MNFALV )
      CALL TNMCDS( 'ENTIER', NBSOM,  MNOUIS )
 9999 RETURN
      END
