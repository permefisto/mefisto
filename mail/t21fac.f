      SUBROUTINE T21FAC( NMOBJT, NUOBJT, MNNSEF, MNXYZS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LE MAILLAGE DE LA SURFACE 2D DE NOM NMOBJT
C -----    DECRIT PAR LES 2 TABLEAUX 'XYZSOMMET ET 'NSEF'

C ENTREE :
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNXYZS : ADRESSE MCN DU TABLEAU 'XYZSOMMET'    A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1996
C ......................................................................
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/gsmenu.inc"

C     DECLARATION DU SUPER-TABLEAU NUMERIQUE MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)

      CHARACTER*(*)     NMOBJT

C     LES VARIABLES LOCALES
      CHARACTER*10      NMSOMM
      CHARACTER*11      KSYMBO
      INTEGER           NOSOEL(1:12)
      REAL              X(12), Y(12), XYB(2)

C     LE TYPE DE L'OBJET
      IF( MCN( MNNSEF + WUTYOB ) .NE. 3 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = NMOBJT // ' NON UNE SURFACE'
         CALL LEREUR
         RETURN
      ENDIF

C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN

C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET ET TG DU TMS 'XYZSOMMET'
      NBSOM = MCN( MNXYZS + WNBSOM )
      MNS   = MNXYZS + WYZSOM - 3
      MNTG  = MNS + 3 * NBSOM

C     REDUCTION DES FACES
      REDUCF = PREDUF * 0.01
      REDUC1 = 1.0 - REDUCF

C     LA BOUCLE SUR LES FACES DE LA SURFACE
C     =====================================
      MNSTS  = MNXYZS + WYZSOM
      QEFMIN = 2
      NEFMIN = 0

      DO 200 NEF = 1, NBEFOB

C        LE NUMERO DES NBSOEF SOMMETS DE L'EF NEF
         CALL NSEFNS( NEF   , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C        LE NOMBRE DE SOMMETS DE CET EF = NCOGEL 3 OU 4

C        LES 3 COORDONNEES DES NCOGEL SOMMETS DE LA FACE
         DO J=1,NCOGEL
            MN   = MNS + 3 * NOSOEL(J)
            X(J) = RMCN( MN     )
            Y(J) = RMCN( MN + 1 )
         ENDDO

         IF( NUEFTG .EQ. 0 .AND.
     %      (IAVNEF .NE. 0 .OR. REDUCF .GT. 0.0) ) THEN

C           REDUCTION DES FACES DE CET EF SANS TG
C           CALCUL DES COORDONNEES XY DU BARYCENTRE
            XX = 0
            YY = 0
            DO J=1,NCOGEL
               XX = XX + X(J)
               YY = YY + Y(J)
            ENDDO
            XYB(1) = XX / NCOGEL
            XYB(2) = YY / NCOGEL

            DO J=1,NCOGEL
               X(J) = X(J) * REDUC1 + XYB(1) * REDUCF
               Y(J) = Y(J) * REDUC1 + XYB(2) * REDUCF
            ENDDO
         ENDIF

C        TRACE DE LA FACE DE NCOGEL ARETES ET DES ARETES
         IF( IAVFAC .EQ. 0 ) THEN
            NCF = -1
         ELSE
C           LA COULEUR DE LA FACE
            NCF = NCOUFA
            IF( LCRITR .GT. 0 ) THEN
C              TRACE DE LA QUALITE DE L'ELEMENT FINI
               CALL QUALEF( NCOGEL,   NOSOEL, NBSOM, RMCN(MNS+3),
     %                      SURFVOLU, QUALIT, IERR )
               IF( IERR .NE. 0 ) RETURN
               IF( QUALIT .LT. QEFMIN ) THEN
C                 REPERAGE DE L'EF DE PLUS MAUVAISE QUALITE
                  QEFMIN = QUALIT
                  NEFMIN = NEF
               ENDIF
C              LA COULEUR DE LA FACE VISUALISE LA QUALITE
               NCF = N1COUL + 9 - INT( 10 * ( 1 - QUALIT ) )
               NCF = MAX( NCF, N1COUL )
            ENDIF
         ENDIF

         IF( IAVARE .EQ. 0 ) THEN
            NCA = -1
         ELSE
            NCA = NCOUAF
         ENDIF

         IF( NBTGEF .GT. 0 ) THEN

C           EF A TG?
            NUEFTG = MCN( MNNSEF + LDAPEF -1 + NEF )
            IF( NUEFTG .GT. 0 ) THEN
C              OUI : LES NUMEROS DES 8 TGS DU QUADRANGLE
C                    SONT RANGEES DANS NOSOEL(5:12)
               NT = 4
               NJ = NCOGEL
               DO 150 J=1,NCOGEL

C                 LE NUMERO DE LA PREMIERE TANGENTE DU SOMMET J
                  NT  = NT + 1
                  NTG = NOSOEL( NT )
                  NJ  = NJ + 1
                  IF( NTG .NE. 0 ) THEN
C                    IL EXISTE UNE TANGENTE
                     MN = MNTG + 3 * ABS(NTG)
                     IF( NTG .GT. 0 ) THEN
                        X(NJ) = RMCN( MN     )
                        Y(NJ) = RMCN( MN + 1 )
                     ELSE
                        X(NJ) = -RMCN( MN     )
                        Y(NJ) = -RMCN( MN + 1 )
                     ENDIF
                  ELSE
C                    PAS DE TANGENTE: COTE DROIT
                     IF( J .LT. NCOGEL ) THEN
                        J1 = J + 1
                     ELSE
                        J1 = 1
                     ENDIF
                     X(NJ) = X(J1) - X(J)
                     Y(NJ) = Y(J1) - Y(J)
                  ENDIF

C                 LE NUMERO DE LA SECONDE TANGENTE DU SOMMET J
                  NT  = NT + 1
                  NTG = NOSOEL( NT )
                  NJ  = NJ + 1
                  IF( NTG .NE. 0 ) THEN
C                    IL EXISTE UNE TANGENTE
                     MN = MNTG + 3 * ABS(NTG)
                     IF( NTG .GT. 0 ) THEN
                        X(NJ) = RMCN( MN     )
                        Y(NJ) = RMCN( MN + 1 )
                     ELSE
                        X(NJ) = -RMCN( MN     )
                        Y(NJ) = -RMCN( MN + 1 )
                     ENDIF
                  ELSE
C                    PAS DE TANGENTE: COTE DROIT
                     IF( J .NE. 1 ) THEN
                        J1 = J - 1
                     ELSE
                        J1 = NCOGEL
                     ENDIF
                     X(NJ) = X(J1) - X(J)
                     Y(NJ) = Y(J1) - Y(J)
                  ENDIF
 150           ENDDO

C              LE TRACE DE L'EF A TG
               CALL FAP32D( NCF, NCA, PREDUF, NCOGEL, X, Y, NOSOEL )
               GOTO 190
            ENDIF
         ENDIF

C        EF DROIT P1 SANS TG
C        NOMBRE D'EPAISSEURS DU TRACE DES ARETES DES FACES
         CALL XVEPAISSEUR( NEPARF )
         CALL FACE2D( NCF, NCA, NCOGEL, X, Y )

C        TRACE EVENTUEL DU NO DE L'EF
 190     IF( IAVNEF .NE. 0 ) THEN
            WRITE( NMSOMM, '(I7)' ) NEF
            KSYMBO = '.' // NMSOMM
            CALL SANSBL( KSYMBO, L )
            CALL SYMBOLE2D( NCONEF, XYB(1), XYB(2), KSYMBO(1:L) )
         ENDIF
 200  ENDDO

C     TRACE EVENTUEL DU NO DES SOMMETS
      IF( IAVNSO .NE. 0 ) THEN
ccc         CALL XVCOULEUR( NCONSO )
ccc         MN = MNXYZS + WYZSOM - 3
         DO J=1,NBSOM
            CALL TRST2D( NCONSO, J, RMCN(MNSTS) )
ccc            MN = MN + 3
ccc            WRITE( NMSOMM, '(I7)' ) J
ccc            KSYMBO = '.' // NMSOMM
ccc            CALL SANSBL( KSYMBO, L )
ccc            CALL SYMBOLE2D( NCONSO, RMCN(MN), RMCN(MN+1), KSYMBO(1:L) )
         ENDDO
      ENDIF

C     TRACE DE LA POIGNEE ET DU NOM DE LA SURFACE
C     LES COORDONNEES DE LA POIGNEE
      IF( IAVNEF .EQ. 0 ) THEN
         DO J=2,NCOGEL
            X(1) = X(1) + X(J)
            Y(1) = Y(1) + Y(J)
         ENDDO
         XYB(1) = X(1) / NCOGEL
         XYB(2) = Y(1) / NCOGEL
      ENDIF
      CALL ITEMS2( XYB, NMOBJT, NUOBJT )

C     TRACE DES ARETES EN ROUGE DE L'EF DE PLUS MAUVAISE QUALITE SI <0.5
C     ------------------------------------------------------------------
      IF( LCRITR .GT. 0 .AND. QEFMIN .LT. 0.5 ) THEN

C        LA COULEUR DE LA FACE VISUALISE LA QUALITE QEFMIN
         NCF = N1COUL + 9 - INT( 10 * ( 1 - QEFMIN ) )
         NCF = MAX( NCF, N1COUL )
         CALL XVEPAISSEUR( NEPARF+2 )

C        LE NUMERO DES NBSOEF SOMMETS DE L'EF NEFMIN
         CALL NSEFNS( NEFMIN, NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )

C        LES 3 COORDONNEES DES NCOGEL SOMMETS DE LA FACE
         MNS = MNXYZS + WYZSOM - 3
         DO J=1,NCOGEL
            MN   = MNS + 3 * NOSOEL(J)
            X(J) = RMCN( MN     )
            Y(J) = RMCN( MN + 1 )
         ENDDO

         IF( NBTGEF .GT. 0 ) THEN
C           EF A TG?
            NUEFTG = MCN( MNNSEF + LDAPEF -1 + NEFMIN )
            IF( NUEFTG .GT. 0 ) THEN

C              OUI : LES NUMEROS DES 8 TGS DU QUADRANGLE
C                    SONT RANGEES DANS NOSOEL(5:12)
               NT = 4
               NJ = NCOGEL
               DO 550 J=1,NCOGEL

C                 LE NUMERO DE LA PREMIERE TANGENTE DU SOMMET J
                  NT  = NT + 1
                  NTG = NOSOEL( NT )
                  NJ  = NJ + 1
                  IF( NTG .NE. 0 ) THEN
C                    IL EXISTE UNE TANGENTE
                     MN = MNTG + 3 * ABS(NTG)
                     IF( NTG .GT. 0 ) THEN
                        X(NJ) = RMCN( MN     )
                        Y(NJ) = RMCN( MN + 1 )
                     ELSE
                        X(NJ) = -RMCN( MN     )
                        Y(NJ) = -RMCN( MN + 1 )
                     ENDIF
                  ELSE
C                    PAS DE TANGENTE: COTE DROIT
                     IF( J .LT. NCOGEL ) THEN
                        J1 = J + 1
                     ELSE
                        J1 = 1
                     ENDIF
                     X(NJ) = X(J1) - X(J)
                     Y(NJ) = Y(J1) - Y(J)
                  ENDIF

C                 LE NUMERO DE LA SECONDE TANGENTE DU SOMMET J
                  NT  = NT + 1
                  NTG = NOSOEL( NT )
                  NJ  = NJ + 1
                  IF( NTG .NE. 0 ) THEN
C                    IL EXISTE UNE TANGENTE
                     MN = MNTG + 3 * ABS(NTG)
                     IF( NTG .GT. 0 ) THEN
                        X(NJ) = RMCN( MN     )
                        Y(NJ) = RMCN( MN + 1 )
                     ELSE
                        X(NJ) = -RMCN( MN     )
                        Y(NJ) = -RMCN( MN + 1 )
                     ENDIF
                  ELSE
C                    PAS DE TANGENTE: COTE DROIT
                     IF( J .NE. 1 ) THEN
                        J1 = J - 1
                     ELSE
                        J1 = NCOGEL
                     ENDIF
                     X(NJ) = X(J1) - X(J)
                     Y(NJ) = Y(J1) - Y(J)
                  ENDIF
 550           ENDDO

C              LE TRACE DE L'EF A TG AVEC DES ARETES ROUGE
               CALL FAP32D( NCF, NCROUG, PREDUF, NCOGEL, X, Y, NOSOEL )
               GOTO 600
            ENDIF
         ENDIF

C        EF SANS TG
         CALL FACE2D( NCF, NCROUG, NCOGEL, X, Y )

C        RETOUR A L'EPAISSEUR NORMALE DU TRACE DES ARETES
 600     CALL XVEPAISSEUR( NEPARF )

      ENDIF

      RETURN
      END
