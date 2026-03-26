        SUBROUTINE SUEX35( NTLXSU , LADEFI ,
     %                     NTFASU , MNFASU , NTSOSU , MNSOSU , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE DES FACES DESIGNEES D'UNE QUADRANGULATION
C -----    STRUCTUREE EXTRAITE D'UNE QUADRANGULATION STRUCTUREE
C          OU SON COMPLEMENTAIRE C-A-D TOUTES LES FACES QUI NE SONT PAS
C          DESIGNEES PAR LES INDICES MIN ET MAX (=>UN TROU DANS LE MAILLAGE)
C          ET LE MAILLAGE EST ALORS NON STRUCTURE
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE A CREER
C LADEFI : TABLEAU DE DEFINITION DE LA SURFACE PARTITIONNEE
C          CF '~/TD/D/A_SURFACE__DEFINITION'
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF '~/TD/D/A___NSEF'
C NTSOSU : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~/TD/D/A___XYZSOMMET'
C IERR   : =0 SI PAS D'ERREUR
C          >0 EN CAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1996
C AJOUTS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Octobre 2011
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*), NOSOEL(12), QUADES
C
      IERR = 0
C
C     QUADRANGLES DESIGNES OU NON (LE COMPLEMENTAIRE)?
      QUADES = LADEFI( WUADES )
C
C     RESTAURATION DU MAILLAGE DE LA SURFACE INITIALE
C     ===============================================
      CALL LXNLOU( NTSURF , LADEFI(WUSUST) , NTLXSF , MN )
      IF( NTLXSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SUEX35: SURFACE INCONNUE'
         ELSE
            KERR(1) = 'SUEX35: UNKNOWN SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     RESTAURATION DES TABLEAUX SOMMETS ET NSEF
      CALL LXTSOU( NTLXSF , 'XYZSOMMET' , NTSOSF , MNSOSF )
      IF( NTSOSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SUEX35: SURFACE SANS SOMMETS'
         ELSE
            KERR(1) = 'SUEX35: SURFACE WITHOUT VERTICES'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      CALL LXTSOU( NTLXSF , 'NSEF' , NTFASF , MNFASF )
      IF( NTFASF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SUEX35: SURFACE SANS NSEF'
         ELSE
            KERR(1) = 'SUEX35: SURFACE WITHOUT NSEF'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LE NOMBRE DE SOMMETS
      NBSOM = MCN( MNSOSF + WNBSOM )
C     LE NOMBRE DE TANGENTES
      NBTGS = MCN( MNSOSF + WNBTGS )
C
C     L'ADRESSE DU DEBUT DES COORDONNEES DES SOMMETS
C     ==============================================
C     TYPE DE L'OBJET : SURFACE
      NTYPVO = MCN(MNFASF + WUTYOB )
      IF( NTYPVO .NE. 3 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SUEX35: TYPE OBJET DIFFERENT DE SURFACE'
         ELSE
            KERR(1) = 'SUEX35: NOT A SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     NUMERO DU TYPE DU MAILLAGE: 4=QUADRANGLE STRUCTURE
      NTYPMA = MCN ( MNFASF + WUTYMA )
      IF ( NTYPMA .NE. 4 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SUEX35: SURFACE NON NON UN QUADRANGLE STRUCTURE'
         ELSE
            KERR(1) = 'SUEX35: SURFACE is NOT a STRUCTURED QUADRANGLE'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     VARIABLE NBARXH : NOMBRE DE SEGMENTS SUIVANT X
      NBXM1 = MCN ( MNFASF + WBARXQ )
C     VARIABLE NBARYH : NOMBRE DE SEGMENTS SUIVANT Y
      NBYM1 = MCN ( MNFASF + WBARYQ )
C
C     CALCUL ET VERIFICATION DES BORNES DE LA SURFACE EXTRAITE
C     ========================================================
      IMIN = LADEFI(WI2MIN)
      IMAX = LADEFI(WI2MAX)
      JMIN = LADEFI(WJ2MIN)
      JMAX = LADEFI(WJ2MAX)
C
      IF( (IMIN.LE.0)     .OR. (JMIN.LE.0)     ) GOTO 9990
      IF( (IMAX.GT.NBXM1) .OR. (JMAX.GT.NBYM1) ) GOTO 9990
      IF( (IMIN.GT.IMAX)  .OR. (JMIN.GT.JMAX)  ) GOTO 9990
C
C     LE NOMBRE DE FACES DE LA SURFACE EXTRAITE
      NBFASU = 0
      NBSOSU = 0
      MNFALV = 0
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNFASF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX   , NY   , NZ   ,
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9900
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
      IF( QUADES .EQ. 0 ) GOTO 105
C
C     TRAITEMENT DES QUADRANGLES DESIGNES PAR LES MIN MAX DES INDICES
C     ---------------------------------------------------------------
C     LA BOUCLE SUR LES QUADRANGLES DESIGNES DU MAILLAGE
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      DO 100 NJ=JMIN,JMAX
         DO 90 NI=IMIN,IMAX
C
C           LE NUMERO DE l'EF
            N = NI + NX * (NJ-1)
C           LE NUMERO DES NBSOEF SOMMETS DE L'EF N
            CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNFASF, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C           LE NOMBRE DE SOMMETS DE CET ELEMENT = NCOGEL
C
            NBFASU = NBFASU + 1
C           LES NUMEROS DES SOMMETS SONT STOCKES
            DO I=1,NCOGEL
               MCN( MNAVS + I ) = NOSOEL( I )
            ENDDO
C           COMPLETION EVENTUELLE AVEC DES ZEROS
            DO I=NCOGEL+1,NBSOEF
               MCN( MNAVS + I ) = 0
            ENDDO
            MNAVS = MNAVS + NBSOEF
C
C           LES NUMEROS DES NBTGEF TANGENTES OU 0 TANGENTE
            IF( NBTGEF .LE. 0 .OR. NUEFTG .LE. 0 ) THEN
C              LES NBTGEF TANGENTES SONT NULLES
               DO I=1,NBTGEF
                  MCN( MNAVT + I ) = 0
               ENDDO
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = 0
            ELSE
               DO I=1,NBTGEF
                  MCN( MNAVT + I ) = NOSOEL( NBSOEF + I )
               ENDDO
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = NUGEEF
            ENDIF
            MNAVT = MNAVT + NBTGEF
            MNAVG = MNAVG + 1
 90      CONTINUE
 100  CONTINUE
      GOTO 300
C
C     TRAITEMENT DES QUADRANGLES NON DESIGNES (=> UN TROU DANS LE MAILLAGE)
C     ---------------------------------------------------------------------
C     LA BOUCLE SUR LES QUADRANGLES DU MAILLAGE
 105  DO 200 NJ=1,NY
         DO 190 NI=1,NX
C
C           SELECTION DU QUADRANGLE OU NON
            IF( JMIN .LE. NJ .AND. NJ .LE. JMAX   .AND.
     %          IMIN .LE. NI .AND. NI .LE. IMAX ) GOTO 190
C
C           LE NUMERO DE L'EF
            N = NI + NX * (NJ-1)
C           LE NUMERO DES NBSOEF SOMMETS DE L'EF N
            CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNFASF, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C           LE NOMBRE DE SOMMETS DE CET ELEMENT = NCOGEL
C
            NBFASU = NBFASU + 1
C           LES NUMEROS DES SOMMETS SONT STOCKES
            DO I=1,NCOGEL
               MCN( MNAVS + I ) = NOSOEL( I )
            ENDDO
C           COMPLETION EVENTUELLE AVEC DES ZEROS
            DO I=NCOGEL+1,NBSOEF
               MCN( MNAVS + I ) = 0
            ENDDO
            MNAVS = MNAVS + NBSOEF
C
C           LES NUMEROS DES NBTGEF TANGENTES OU 0 TANGENTE
            IF( NBTGEF .LE. 0 .OR. NUEFTG .LE. 0 ) THEN
C              LES NBTGEF TANGENTES SONT NULLES
               DO I=1,NBTGEF
                  MCN( MNAVT + I ) = 0
               ENDDO
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = 0
            ELSE
               DO I=1,NBTGEF
                  MCN( MNAVT + I ) = NOSOEL( NBSOEF + I )
               ENDDO
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = NUGEEF
            ENDIF
            MNAVT = MNAVT + NBTGEF
            MNAVG = MNAVG + 1
 190     CONTINUE
 200  CONTINUE
C
C     RENUMEROTATION DES SOMMETS ET TANGENTES DES EF EXTRAITS
C     -------------------------------------------------------
 300  CALL REEFEX( 3,      NTLXSU, NBSOEF, NBTGEF,
     %             NBSOM,  NBTGS,  MNSOSF,
     %             NBFASU, MNFALV, MNFALV+NBEFOB*NBSOEF,
     %             MNFALV+NBEFOB*(NBSOEF+NBTGEF),
     %             NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
      CALL TNMCDS( 'ENTIER', MXFALV, MNFALV )
C
      IF( QUADES .NE. 0 ) THEN
C
C        RESTRUCTURATION DU MAILLAGE: QUADRANGLE STRUCTURE
C        -------------------------------------------------
         NBEFOB = MCN( MNFASU + WBEFOB )
         NBSOEF = MCN( MNFASU + WBSOEF )
         NBEFAP = MCN( MNFASU + WBEFAP )
         NBEFTG = MCN( MNFASU + WBEFTG )
         NBTGEF = MCN( MNFASU + WBTGEF )
C
C        TRANSLATION DE 4*NBEFOB AU DELA DU NUMERO DES SOMMETS
         NI = WUSOEF + NBSOEF * NBEFOB - WBARYQ - 1
         N  = MNFASU + WUSOEF + NBSOEF * NBEFOB
         DO I = N, N+NBEFAP+NBEFTG*(1+NBTGEF)-1
            MCN(I-NI) = MCN(I)
         ENDDO
C
C        NUMERO DU TYPE DU MAILLAGE : QUADRANGLE STRUCTURE
         MCN( MNFASU + WUTYMA ) = 4
C        VARIABLE NBARXQ : NOMBRE DE SEGMENTS SUIVANT X
         MCN( MNFASU + WBARXQ ) = IMAX-IMIN+1
C        VARIABLE NBARYQ : NOMBRE DE SEGMENTS SUIVANT Y
         MCN( MNFASU + WBARYQ ) = JMAX-JMIN+1
C
C        REDUCTION DU TMS 'NSEF'
         CALL TAMSRA( NTFASU, WBARYQ+1+NBEFAP+NBEFTG*(1+NBTGEF) )
C
      ENDIF
C
 9900 RETURN
C
C     ERREUR DANS LES INDICES
 9990 NBLGRC(NRERR) = 3
      WRITE(KERR(MXLGER-1)(1:10),'(I10)') NBXM1
      WRITE(KERR(MXLGER  )(1:10),'(I10)') NBYM1
      IF( LANGAG .EQ. 0 ) THEN
      KERR(1)='SUEX35:MAUVAISE DEFINITION DES BORNES MIN MAX DES ARETES'
      KERR(2)='I COMPRIS ENTRE  1 ET' // KERR(MXLGER-1)(1:10)
      KERR(3)='J COMPRIS ENTRE  1 ET' // KERR(MXLGER  )(1:10)
      ELSE
      KERR(1)='SUEX35: BAD DEFINITION of BOUNDS MIN MAX of EDGES'
      KERR(2)='I INCLUDED BETWEEN 1 and' // KERR(MXLGER-1)(1:10)
      KERR(3)='J INCLUDED BETWEEN 1 and' // KERR(MXLGER  )(1:10)
      ENDIF
      CALL LEREUR
      IERR = 1
      RETURN
      END
