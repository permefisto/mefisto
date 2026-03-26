      SUBROUTINE LIEX30( NTLXLI , LADEFI ,
     %                   NTARLI , MNARLI , NTSOLI , MNSOLI , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AMELIORATION DE LA QUALITE DES ARETES D'UNE LIGNE
C -----    A PARTIR DE LA FONCTION UTILISATEUR TAILLE_IDEALE(x,y,z)
C          ou EDGE_LENGTH(x,y,z) ou sinon LA TAILLE D'ARETE PAR DEFAUT
C          avec PERTE DES EVENTUELLES TANGENTES
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
C          CF '~/TD/D/A___XYZSOMMET'
C IERR   : 0 SI PAS D'ERREUR
C          1 SI LIGNE INITIALE SANS NSEF
C          2 SI LIGNE INITIALE SANS SOMMETS
c          3 SI NUMERO D'ARETE INCORRECT
C          4 SI PAS DE FONCTION TAILLE_ARETE et DARETE=0D0
C          5 SI TAILLE ET DARETE DONNENT 0D0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC $ St PIERRE du PERRAY  Octobre 2015
C2345X7..............................................................012
      IMPLICIT          INTEGER(W)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/darete.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      INTEGER           NOSOEL(4)
      DOUBLE PRECISION  XYZ1(3), XYZ2(3), D1, D2, DM, D
      INTRINSIC         SQRT

C     DARETE EST LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE dans incl/darete.inc
C     DARETE=0D0 INDIQUE UNE NON INITIALISATION DE CETTE VALEUR
C     SA VALEUR EST REMPLACEE DANS LES CALCULS PAR LA VALEUR
C     DE LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' ou 'EDGE_LENGTH(X,Y,Z)'
C     AU POINT (X,Y,Z) FOURNIE PAR L'UTILISATEUR SELON LE LANGAGE LU

C     EXISTENCE OU NON DE LA FONCTION TAILLE_IDEALE(x,y,z) ou EDGE_LENGTH(x,y,z)
      NOFOTI = NOFOTIEL()

      IF( NOFOTI .LE. 0 .AND. DARETE .EQ. 0D0 ) THEN
         IERR = 4
         RETURN
      ENDIF
C
C     LA LIGNE INITIALE
C     =================
C     LE NOM DE CETTE LIGNE
      NULIIN = LADEFI( WULIIN )
C     LE TABLEAU LEXIQUE DE CETTE LIGNE
      CALL LXNLOU( NTLIGN, NULIIN, NTLXLG, MN )
      IF( NTLXLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE INITIALE INCONNUE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE TABLEAU 'NSEF' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG, 'NSEF', NTARLG, MNARLG )
      IF( NTARLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS NSEF'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG, 'XYZSOMMET', NTSOLG, MNSOLG )
      IF( NTSOLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS SOMMETS'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      NBSOM  = MCN( MNSOLG + WNBSOM )
      NBTGS  = MCN( MNSOLG + WNBTGS )
      NBCOOR = MCN( MNSOLG + WBCOOR )
C
C     LE NOMBRE D'ARETES DE LA LIGNE FINALE
      NBARLI = 0
      NBSOLI = 0
      MNARLV = 0
      MNXYZS = 0
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE INITIAL
      CALL NSEFPA( MCN(MNARLG),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU TEMPORAIRE
C     DE STOCKAGE DES 2 SOMMETS DE CHAQUE ARETE EXTRAITE
      MXARLV = 10 * NBEFOB * NBSOEF
      CALL TNMCDC( 'ENTIER', MXARLV, MNARLV )
      MNAVS = MNARLV - 1

      MNXYZS = 0
      MXXYZS = 10 * NBSOM * NBCOOR
      CALL TNMCDC( 'REEL', MXXYZS, MNXYZS )
      MNS = MNXYZS - 1
C
C     LA BOUCLE SUR LES ARETES DU MAILLAGE DE LA LIGNE
C     ------------------------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DES ARETES DE LA LIGNE
      DO 100 N = 1, NBEFOB

C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNARLG, NX, NY, NZ ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )

C        L'ARETE N EST ELLE TROP COURTE OU TROP LONGUE?
         MNS1 = MNSOLG + WYZSOM + NBCOOR * ( NOSOEL(1) - 1 ) - 1
         MNS2 = MNSOLG + WYZSOM + NBCOOR * ( NOSOEL(2) - 1 ) - 1
         D    = 0D0
         DO K=1,3
            XYZ1(K) = RMCN( MNS1 + K )
            XYZ2(K) = RMCN( MNS2 + K )
            D       = D + ( XYZ2(K) - XYZ1(K) ) ** 2 
         ENDDO

C        LA TAILLE DE L'ARETE DROITE
         D = SQRT( D ) 

C        TAILLE DE L'ARETE IDEALE AU SOMMET 1 ET 2 DE L'ARETE
         D1 = 0D0
         D2 = 0D0

         IF( NOFOTI .GT. 0 ) THEN

            CALL FONVAL( NOFOTI, 3, XYZ1, NCODEV, D1 )
            IF( NCODEV .LE. 0 ) THEN
               D1   = 0D0
            ENDIF
C
            CALL FONVAL( NOFOTI, 3, XYZ2, NCODEV, D2 )
            IF( NCODEV .LE. 0 ) THEN
               D2   = 0D0
            ENDIF

         ELSE

            D1 = DARETE
            D2 = DARETE

         ENDIF

         IF( D1 .LE. 0D0 .OR. D2 .LE. 0D0 ) THEN
            IERR = 5
            GOTO 9000
         ENDIF

C        LONGUEUR MOYENNE
         DM = ( D1 + D2 ) / 2D0

C        NOMBRE D'ARETES SOUHAITEES
         NBA = NINT( D / DM )
         IF( NBA .GE. 1 ) THEN

C           PREMIER SOMMET DE L'ARETE
            NBSOLI = NBSOLI + 1
            DO K = 1, 3
               RMCN( MNS + K ) = REAL( XYZ1(K) )
            ENDDO
            MNS = MNS + 3

            DO NA = 1, NBA

C              UNE ARETE DE PLUS
               NBARLI = NBARLI + 1
               MCN( MNAVS + 1 ) = NBSOLI

C              LE SECOND SOMMET DE L'ARETE
               NBSOLI = NBSOLI + 1
               MCN( MNAVS + 2 ) = NBSOLI

               IF( NA .NE. NBA ) THEN
                  DO K = 1, 3
                     RMCN( MNS + K ) = REAL( XYZ1(K)
     %                               + NA * ( XYZ2(K) - XYZ1(K) ) / NBA)
                  ENDDO
               ELSE
                  DO K = 1, 3
                     RMCN( MNS + K ) = REAL( XYZ2(K) )
                  ENDDO
               ENDIF

               MNS = MNS + NBCOOR
               IF( MNS+NBCOOR .GT. MNXYZS+MXXYZS ) THEN
C                 AUGMENTATION DE LA TAILLE DU TABLEAU
                  K = NBCOOR * NBSOLI
                  CALL TNMCAU( 'REEL', MXXYZS, 2*MXXYZS, K, MNXYZS )
                  MXXYZS = 2 * MXXYZS
                  MNS    = MNXYZS - 1 + K
               ENDIF

               MNAVS = MNAVS + NBSOEF
               IF( MNAVS+NBSOEF .GT. MNARLV+MXARLV ) THEN
C                 AUGMENTATION DE LA TAILLE DU TABLEAU
                  K = NBSOEF * NBARLI
                  CALL TNMCAU( 'ENTIER', MXARLV, 2*MXARLV, K, MNARLV )
                  MXARLV = 2 * MXARLV
                  MNAVS  = MNARLV -1 + K
               ENDIF

            ENDDO

         ENDIF

 100  CONTINUE

C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     -----------------------------------
      CALL LXTNDC( NTLXLI, 'XYZSOMMET', 'MOTS',
     %             WYZSOM+NBCOOR*NBSOLI )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )

C     LE NOMBRE DE SOMMETS DE CE MAILLAGE
      MCN( MNSOLI + WNBSOM ) = NBSOLI
C     LE NOMBRE DE TANGENTES DE CE MAILLAGE
      MCN( MNSOLI + WNBTGS ) = 0
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOLI + WBCOOR ) = NBCOOR

C     COPIE DES 3 COORDONNEES DES NBSOLI SOMMETS
      MN1 = MNXYZS - 1
      MN2 = MNSOLI + WYZSOM - 1
      DO K=1,3*NBSOLI
         RMCN( MN2 + K ) = RMCN( MN1 + K )
      ENDDO

C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOLI) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOLI + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C     CONSTRUCTION DU TABLEAU 'NSEF'
C     ------------------------------
      CALL LXTNDC( NTLXLI, 'NSEF', 'MOTS', WUSOEF + NBARLI * NBSOEF )
      CALL LXTSOU( NTLXLI, 'NSEF',  NTARLI, MNARLI )

C     LE TYPE DE L'OBJET
      MCN( MNARLI + WUTYOB ) = 2
C     LIGNE UNKNOWN SUR SA FERMETURE
      MCN( MNARLI + WUTFMA ) = -1
C     LE NOMBRE DE SOMMETS DES EF
      MCN( MNARLI + WBSOEF ) = 2
C     LE NOMBRE DE TANGENTES DES EF   ICI PAS DE TG STOCKEES
      MCN( MNARLI + WBTGEF ) = 0
C     LE NOMBRE D'EF DU MAILLAGE
      MCN( MNARLI + WBEFOB ) = NBARLI
C     LE NOMBRE D EF POINT'ES
      MCN( MNARLI + WBEFAP ) = 0
C     LE NOMBRE D EF A TG
      MCN( MNARLI + WBEFTG ) = 0

C     LE TYPE DU MAILLAGE: ICI NON STRUCTURE
      MCN( MNARLI + WUTYMA ) = 0

C     COPIE DES 2 NUMEROS DES SOMMETS DES NBARLI ARETES
      MN1 = MNARLV - 1
      MN2 = MNARLI + WUSOEF - 1
      DO K=1,2*NBARLI
         RMCN( MN2 + K ) = RMCN( MN1 + K )
      ENDDO

C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARLI) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARLI + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C     DESTRUCTION DES TABLEAUX AUXILIAIRES
 9000 CALL TNMCDS( 'ENTIER', MXARLV, MNARLV )
      CALL TNMCDS( 'REEL',   MXXYZS, MNXYZS )

C     TENTATIVE DE STRUCTURATION DE LA LIGNE EXTRAITE
      IF( IERR .EQ. 0 ) THEN
         CALL LIGSTR( NTLXLI, NTARLI, MNARLI, NTSOLI, MNSOLI, IER )
      ENDIF
C
      RETURN
      END
