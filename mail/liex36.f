      SUBROUTINE LIEX36( NTLXLI , LADEFI ,
     %                   NTARLI , MNARLI , NTSOLI , MNSOLI , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXTRAIRE UNE LIGNE D'UNE LIGNE A PARTIR D'UN CRITERE LOGIQUE
C -----
C
C ENTREES:
C --------
C NTLXLI : NUMERO DU TABLEAU TS DU LEXIQUE DE LA LIGNE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA LIGNE
C          CF ~/TD/D/A_LIGNE__DEFINITION
C
C SORTIES:
C --------
C NTARLI : NUMERO      DU TMS 'NSEF' DE LA LIGNE
C MNARLI : ADRESSE MCN DU TMS 'NSEF' DE LA LIGNE
C          CF ~/TD/D/A___NSEF
C NTSOLI : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C IERR   : 0 SI PAS D'ERREUR
C          1 SI LIGNE INITIALE SANS NSEF
C          2 SI LIGNE INITIALE SANS SOMMETS
C          3 SI FONCTION INCONNUE
C          4 SI AUCUNE ARETE EXTRAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      IMPLICIT           INTEGER(W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      INTEGER           NOSOEL(4)
      DOUBLE PRECISION  DXYZ(3),DOUI
      LOGICAL           OUI
C
C     LA LIGNE INITIALE
C     =================
C     LE NOM DE CETTE LIGNE
      NULIIN = LADEFI( WULIIN )
C     LE TABLEAU LEXIQUE DE CETTE LIGNE
      CALL LXNLOU( NTLIGN , NULIIN , NTLXLG , MN )
      IF( NTLXLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE INITIALE INCONNUE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'NSEF' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG , 'NSEF' , NTARLG , MNARLG )
      IF( NTARLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS NSEF'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG , 'XYZSOMMET' , NTSOLG , MNSOLG )
      IF( NTSOLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS SOMMETS'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      NBSOM = MCN( MNSOLG + WNBSOM )
      NBTGS = MCN( MNSOLG + WNBTGS )
C
C
C     LE NUMERO DE LA FONCTION
      NUFOCL = LADEFI( WUFOCL )
      IF( NUFOCL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FONCTION CRITERE INCONNUE'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LE NOMBRE D'ARETES DE LA LIGNE EXTRAITE
      NBARLI = 0
      NBSOLI = 0
      MNARLV = 0
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNARLG) ,
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
      MNS = MNSOLG + WYZSOM - 3
      DO 10 N=1,NBSOM
C
C        LES 3 COORDONNEES EN DOUBLE PRECISION DU SOMMET N
         MN = MNS + 3 * N
         DXYZ(1) = RMCN( MN   )
         DXYZ(2) = RMCN( MN+1 )
         DXYZ(3) = RMCN( MN+2 )
C
C        LE SOMMET N VERIFIE T IL LE CRITERE ?
         CALL FONVAL( NUFOCL , 3 , DXYZ , NCODEV , DOUI )
         IF( NCODEV .EQ. 0 ) DOUI = 0
C        1 SI CRITERE VERIFIE , 0 SINON
         MCN( MNO+N ) = NINT( DOUI )
 10   CONTINUE
C
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU TEMPORAIRE
C     DE STOCKAGE DES 2 SOMMETS DE CHAQUE ARETE EXTRAITE
C                 DES 2 TANGENTES ET DU CODE GEOMETRIQUE
      MXARLV = NBEFOB * ( NBSOEF + NBTGEF + 1 )
      CALL TNMCDC( 'ENTIER' , MXARLV , MNARLV )
      MNAVS = MNARLV - 1
      MNAVT = MNAVS + NBEFOB * NBSOEF
      MNAVG = MNAVT + NBEFOB * NBTGEF
C
C     LA BOUCLE SUR LES NO SOMMET DU MAILLAGE
C     -----------------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      DO 100 N=1,NBEFOB
C        LE NUMERO DES NBSOEF SOMMETS DU EF N
         CALL NSEFNS( N      , NUTYMA , NBSOEF , NBTGEF,
     %                LDAPEF , LDNGEF , LDTGEF,
     %                MNARLG , NX , NY , NZ ,
     %                NCOGEL , NUGEEF , NUEFTG, NOSOEL , IERR )
C
         OUI  = .TRUE.
         DO 20 I=1,2
C           LE SOMMET VERIFIE T IL LE CRITERE ?
            OUI = OUI .AND. (MCN(MNO+NOSOEL(I)).EQ. 1)
 20      CONTINUE
C
         IF( OUI ) THEN
C
C           L'EF VERIFIE LE CRITERE
            NBARLI = NBARLI + 1
C
C           LES NUMEROS DES SOMMETS SONT STOCKES
            DO 30 I=1,2
               MCN( MNAVS + I ) = NOSOEL( I )
 30         CONTINUE
            MNAVS = MNAVS + 2
C
C           LES NUMEROS DES 2 TANGENTES OU 0
            IF( NBTGEF .EQ. 0 .OR. NUEFTG .EQ. 0 ) THEN
C              LES 2 TANGENTES SONT NULLES
               DO 40 I=1,2
                  MCN( MNAVT + I ) = 0
 40            CONTINUE
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = 0
            ELSE
               DO 42 I=1,2
                  MCN( MNAVT + I ) = NOSOEL( 2 + I )
 42            CONTINUE
C              CODE GEOMETRIQUE
               MCN( MNAVG + 1 ) = NUGEEF
            ENDIF
            MNAVT = MNAVT + 2
            MNAVG = MNAVG + 1
         ENDIF
 100  CONTINUE
C
C     RENUMEROTATION DES SOMMETS ET TANGENTES DES EF EXTRAITS
      CALL REEFEX( 2,      NTLXLI, NBSOEF, NBTGEF,
     %             NBSOM,  NBTGS,  MNSOLG,
     %             NBARLI, MNARLV, MNARLV+NBEFOB*NBSOEF,
     %             MNARLV+NBEFOB*(NBSOEF+NBTGEF),
     %             NTARLI, MNARLI, NTSOLI, MNSOLI, IERR )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
      CALL TNMCDS( 'ENTIER', MXARLV, MNARLV )
      CALL TNMCDS( 'ENTIER', NBSOM,  MNOUIS )
C
C     TENTATIVE DE STRUCTURATION DE LA LIGNE EXTRAITE
      IF( IERR .EQ. 0 ) THEN
         CALL LIGSTR( NTLXLI, NTARLI, MNARLI, NTSOLI, MNSOLI, IER )
      ENDIF
C
 9999 RETURN
      END
