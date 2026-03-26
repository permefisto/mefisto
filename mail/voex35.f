      SUBROUTINE VOEX35( NTLXSU , LADEFI ,
     %                   NTCUVO , MNCUVO , NTSOCU , MNSOCU , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE DES HEXAEDRES DESIGNES D'UNE HEXAEDRISATION
C -----    STRUCTUREE EXTRAITE D'UNE HEXAEDRISATION STRUCTUREE
C          OU SON COMPLEMENTAIRE C-A-D TOUTES LES HEXAEDRES QUI NE SONT PAS
C          DESIGNES PAR LES INDICES MIN ET MAX (=>UN TROU DANS LE MAILLAGE)
C          ET LE MAILLAGE EST ALORS NON STRUCTURE
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME A CREER
C LADEFI : TABLEAU DE DEFINITION DU VOLUME PARTITIONNE
C          CF '~/TD/D/A_VOLUME__DEFINITION'
C
C SORTIES:
C --------
C NTCUVO : NUMERO      DU TMS 'NSEF' DES NUMEROS DES CUBES
C MNCUVO : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES CUBES
C          CF '~/TD/D/A___NSEF'
C NTSOCU : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOCU : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
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
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))
      INTEGER        LADEFI(0:*), NOSOEL(64), HEXDES
C
      IERR = 0
C
C     HEXAEDRES DESIGNES OU NON (LE COMPLEMENTAIRE)?
      HEXDES = LADEFI( WEXDES )
C
C     RESTAURATION DU MAILLAGE DU VOLUME INITIAL
C     ==========================================
      CALL LXNLOU( NTVOLU, LADEFI(WUVOIN), NTLXVO, MN )
      IF( NTLXVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOEX35: VOLUME INCONNU'
         ELSE
            KERR(1) = 'VOEX35: UNKNOWN VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     RESTAURATION DES TABLEAUX SOMMETS ET NSEF
      CALL LXTSOU( NTLXVO , 'XYZSOMMET' , NTSOVL , MNSOVL )
      IF( NTSOVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOEX35: VOLUME SANS XYZ SOMMETS'
         ELSE
            KERR(1) = 'VOEX35: VOLUME WITHOUT XYZ VERTICES'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      CALL LXTSOU( NTLXVO , 'NSEF' , NTCUVL , MNCUVL )
      IF( NTCUVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOEX35: VOLUME SANS NSEF'
         ELSE
            KERR(1) = 'VOEX35: VOLUME WITHOUT NSEF'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS
      NBSOM = MCN( MNSOVL + WNBSOM )
C     LE NOMBRE DE TANGENTES
      NBTGS = MCN( MNSOVL + WNBTGS )
C
C     L'ADRESSE DU DEBUT DES COORDONNEES DES SOMMETS
C     ==============================================
C     TYPE DE L'OBJET : VOLUME
      NTYPVO = MCN(MNCUVL + WUTYOB )
      IF( NTYPVO .NE. 4 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOEX35: N''EST PAS UN VOLUME'
         ELSE
            KERR(1) = 'VOEX35: NOT A VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     NUMERO DU TYPE DU MAILLAGE : HEXAEDRE STRUCTURE
      NTYPMA = MCN ( MNCUVL + WUTYMA )
      IF (NTYPMA .NE. 7) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOEX35: VOLUME NON STRUCTURE'
         ELSE
            KERR(1) = 'VOEX35: NOT A STRUCTURED HEXAHEDRON'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     VARIABLE NBARXH : NOMBRE DE SEGMENTS SUIVANT X
      NBAX = MCN ( MNCUVL + WBARXH )
C     VARIABLE NBARYH : NOMBRE DE SEGMENTS SUIVANT Y
      NBAY = MCN ( MNCUVL + WBARYH )
C     VARIABLE NBARZH : NOMBRE DE SEGMENTS SUIVANT Y
      NBAZ = MCN ( MNCUVL + WBARZH )
C
C     CALCUL ET VERIFICATION DES BORNES DU VOLUME EXTRAIT
C     ===================================================
      IMIN = LADEFI(WI3MIN)
      IMAX = LADEFI(WI3MAX)
      JMIN = LADEFI(WJ3MIN)
      JMAX = LADEFI(WJ3MAX)
      KMIN = LADEFI(WK3MIN)
      KMAX = LADEFI(WK3MAX)
C
      IF( (IMIN.LE.0)    .OR. (JMIN.LE.0)    .OR.
     %    (KMIN.LE.0)    ) GOTO 9990
      IF( (IMAX.GT.NBAX) .OR. (JMAX.GT.NBAY) .OR.
     %    (KMAX.GT.NBAZ) ) GOTO 9990
      IF( (IMIN.GT.IMAX) .OR. (JMIN.GT.JMAX) .OR.
     %    (KMIN.GT.KMAX) ) GOTO 9990
C
C     LE NOMBRE DE CUBES DU VOLUME EXTRAITE
      NBCUVO = 0
      NBSOSU = 0
      MNFALV = 0
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNCUVL),
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
      IF( HEXDES .EQ. 0 ) GOTO 105
C
C     TRAITEMENT DES HEXAEDRES DESIGNES PAR LES MIN MAX DES INDICES
C     -------------------------------------------------------------
C     LA BOUCLE SUR LES HEXAEDRES DESIGNES DU MAILLAGE
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DES EF DU MAILLAGE
      DO 100 NK=KMIN,KMAX
         DO  90 NJ=JMIN,JMAX
            DO 80 NI=IMIN,IMAX
C
C              LE NUMERO DE l'EF
               N = NI + NX * ( (NJ-1) + NY * (NK-1) )
C              LE NUMERO DES NBSOEF SOMMETS DE L'HEXAEDRE N
               CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNCUVL, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
               NBCUVO = NBCUVO + 1
C              LES NUMEROS DES 8 SOMMETS SONT STOCKES
               DO I=1,8
                  MCN( MNAVS + I ) = NOSOEL( I )
               ENDDO
               MNAVS = MNAVS + 8
C
C              LES NUMEROS DES NBTGEF TANGENTES OU 0 TANGENTE
               IF( NBTGEF .EQ. 0 .OR. NUEFTG .EQ. 0 ) THEN
C                 LES NBTGEF TANGENTES SONT NULLES
                  DO I=1,NBTGEF
                     MCN( MNAVT + I ) = 0
                  ENDDO
C                 CODE GEOMETRIQUE
                  MCN( MNAVG + 1 ) = 0
               ELSE
                  DO I=1,NBTGEF
                     MCN( MNAVT + I ) = NOSOEL( NBSOEF + I )
                  ENDDO
C                 CODE GEOMETRIQUE
                  MCN( MNAVG + 1 ) = NUGEEF
               ENDIF
               MNAVT = MNAVT + NBTGEF
               MNAVG = MNAVG + 1
  80        CONTINUE
  90     CONTINUE
 100  CONTINUE
      GOTO 300
C
C     TRAITEMENT DES HEXAEDRES NON DESIGNES (=> UN TROU DANS LE MAILLAGE)
C     -------------------------------------------------------------------
C     LA BOUCLE SUR LES HEXAEDRES DESIGNES DU MAILLAGE
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DES EF DU MAILLAGE
 105  DO 200 NK=1,NZ
         DO 190 NJ=1,NY
            DO 180 NI=1,NX
C
C           SELECTION DE L'HEXAEDRE OU NON
            IF( KMIN .LE. NK .AND. NK .LE. KMAX   .AND.
     %          JMIN .LE. NJ .AND. NJ .LE. JMAX   .AND.
     %          IMIN .LE. NI .AND. NI .LE. IMAX ) GOTO 180
C
C              LE NUMERO DE l'EF
               N = NI + NX * ( (NJ-1) + NY * (NK-1) )
C              LE NUMERO DES NBSOEF SOMMETS DE L'HEXAEDRE N
               CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNCUVL, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
               NBCUVO = NBCUVO + 1
C              LES NUMEROS DES 8 SOMMETS SONT STOCKES
               DO I=1,8
                  MCN( MNAVS + I ) = NOSOEL( I )
               ENDDO
               MNAVS = MNAVS + 8
C
C              LES NUMEROS DES NBTGEF TANGENTES OU 0 TANGENTE
               IF( NBTGEF .EQ. 0 .OR. NUEFTG .EQ. 0 ) THEN
C                 LES NBTGEF TANGENTES SONT NULLES
                  DO I=1,NBTGEF
                     MCN( MNAVT + I ) = 0
                  ENDDO
C                 CODE GEOMETRIQUE
                  MCN( MNAVG + 1 ) = 0
               ELSE
                  DO I=1,NBTGEF
                     MCN( MNAVT + I ) = NOSOEL( NBSOEF + I )
                  ENDDO
C                 CODE GEOMETRIQUE
                  MCN( MNAVG + 1 ) = NUGEEF
               ENDIF
               MNAVT = MNAVT + NBTGEF
               MNAVG = MNAVG + 1
 180        CONTINUE
 190     CONTINUE
 200  CONTINUE
C
C     RENUMEROTATION DES SOMMETS ET TANGENTES DES EF EXTRAITS
C     -------------------------------------------------------
 300  CALL REEFEX( 4,      NTLXSU, NBSOEF, NBTGEF,
     %             NBSOM,  NBTGS,  MNSOVL,
     %             NBCUVO, MNFALV, MNFALV+NBEFOB*NBSOEF,
     %             MNFALV+NBEFOB*(NBSOEF+NBTGEF),
     %             NTCUVO, MNCUVO, NTSOCU, MNSOCU, IERR )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
      CALL TNMCDS( 'ENTIER', MXFALV, MNFALV )
C
      IF( HEXDES .NE. 0 ) THEN
C
C        RESTRUCTURATION DU MAILLAGE: HEXAEDRE STRUCTURE
C        -----------------------------------------------
         NBEFOB = MCN( MNCUVO + WBEFOB )
         NBSOEF = MCN( MNCUVO + WBSOEF )
         NBEFAP = MCN( MNCUVO + WBEFAP )
         NBEFTG = MCN( MNCUVO + WBEFTG )
         NBTGEF = MCN( MNCUVO + WBTGEF )
C
C        TRANSLATION DE NBSOEF*NBEFOB AU DELA DU NUMERO DES SOMMETS
         J = WUSOEF + NBSOEF * NBEFOB - WBARZH - 1
         N = MNCUVO + WUSOEF + NBSOEF * NBEFOB
         DO I = N, N+NBEFAP+NBEFTG*(1+NBTGEF)-1
            MCN(I-J) = MCN(I)
         ENDDO
C
C        NUMERO DU TYPE DU MAILLAGE : HEXAEDRE STRUCTURE
         MCN( MNCUVO + WUTYMA ) = 7
C        VARIABLE NBARXQ : NOMBRE DE SEGMENTS SUIVANT X OU I
         MCN( MNCUVO + WBARXH ) = IMAX-IMIN+1
C        VARIABLE NBARYQ : NOMBRE DE SEGMENTS SUIVANT Y OU J
         MCN( MNCUVO + WBARYH ) = JMAX-JMIN+1
C        VARIABLE NBARZQ : NOMBRE DE SEGMENTS SUIVANT Z OU K
         MCN( MNCUVO + WBARZH ) = KMAX-KMIN+1
C
C        REDUCTION DU TMS 'NSEF'
         CALL TAMSRA( NTCUVO, WBARZH+1+NBEFAP+NBEFTG*(1+NBTGEF) )
C
      ENDIF
C
 9900 RETURN
C
C     ERREUR DANS LES INDICES
 9990 NBLGRC(NRERR) = 3
      WRITE(KERR(MXLGER-2)(1:10),'(I10)') NBAX
      WRITE(KERR(MXLGER-1)(1:10),'(I10)') NBAY
      WRITE(KERR(MXLGER  )(1:10),'(I10)') NBAZ
         IF( LANGAG .EQ. 0 ) THEN
      KERR(1)='VOEX35:MAUVAISE DONNEE DES BORNES MIN MAX DES ARETES'
      KERR(2)='I A COMPRENDRE ENTRE 1 ET' // KERR(MXLGER-2)(1:10)
      KERR(3)='J A COMPRENDRE ENTRE 1 ET' // KERR(MXLGER-1)(1:10)
      KERR(3)='K A COMPRENDRE ENTRE 1 ET' // KERR(MXLGER  )(1:10)
         ELSE
      KERR(1)='VOEX35: BAD BOUNDS of MIN MAX of EDGES'
      KERR(2)='I MUST BE INCLUDED BETWEEN 1 and' // KERR(MXLGER-2)(1:10)
      KERR(3)='J MUST BE INCLUDED BETWEEN 1 and' // KERR(MXLGER-1)(1:10)
      KERR(3)='K MUST BE INCLUDED BETWEEN 1 and' // KERR(MXLGER  )(1:10)
         ENDIF
      CALL LEREUR
      IERR = 1
      RETURN
      END
