      SUBROUTINE THDLFX( NCDLF0, NTDL,   NDIM,
     %                   NBTYEL, MNNPEF, NDPGST,
     %                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                   NBDLFX, MONDLFX, MNNDLFX, MNVDLFX, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER TOUS LES DEGRES DE LIBERTE FIXES (CONDITION DIRICHLET)
C -----    ET STOCKER LEUR NUMERO ET VALEUR

C ENTREES:
C --------
C NCDLF0 : CODE 0 SI TOUT DL FIXE NON NUL DOIT DECLENCHER IERR NON NUL
C               1 SINON
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE (DL) DES TEMPERATURES
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3 ou 6
C         (SI AXISYMETRIE NDIM=2 X=>R>=0 Y=>Z Z=0)
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES PLSV
C NUMAOB : NUMERO MAXIMAL DES PLSV
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION

C SORTIES:
C --------
C NBDLFX : NOMBRE DE DL SUPPORT D'UN DEGRE DE LIBERTE FIXE
C MONDLFX: NOMBRE DE MOTS DECLARES DU TABLEAU MC NO et VALEUR DES DL FIXES
C MNNDLFX: ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES, 0 SINON
C MNVDLFX: ADRESSE MCN DU TABLEAU MC DES VALEURS DES DL FIXES, 0 SINON
C IERR   : =0 SI PAS D'ERREUR
C          >1 SINON ET LES ADRESSES MNNDLFX ET MNVDLFX SONT NULLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET   ANALYSE NUMERIQUE UPMC PARIS       AOUT 1998
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/donthe.inc"
      include"./incl/xyzext.inc"

      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)

      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), MNDOEL(4)

      DOUBLE PRECISION  RELMIN, D, XYZP(3), DPARA(10), DLMIN, DLMAX

      IERR   = 0
      MOREE2 = MOTVAR(6)

      IF( NBDLFX .GT. 0 ) THEN
C        DESTRUCTION DES TABLEAUX DECLARES AUPARAVANT
         IF( MNNDLFX .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLFX, MNNDLFX )
         IF( MNVDLFX .GT. 0 ) CALL TNMCDS( 'REEL2',  NBDLFX, MNVDLFX )
      ENDIF

      NBDLFX  = 0
      MONDLFX = 0
      MNNDLFX = 0
      MNVDLFX = 0

CCC   SI CONDITION DE NEUMANN SUR LE BORD DU 6-CUBES
CCC   ALORS DECOMMENTER LA LIGNE CI-DESSOUS!
CCC   IF( NDIM .EQ. 6 ) RETURN

C     CONDITION DE DIRICHLET HOMOGENE SUR LA FRONTIERE
C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      MONDLFX = NTDL
      CALL TNMCDC( 'ENTIER', MONDLFX, MNNDLFX )
      MNNDX = MNNDLFX - 1
      CALL AZEROI( MONDLFX, MCN(MNNDLFX) )

C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
      CALL TNMCDC( 'REEL2', MONDLFX, MNVDLFX )
      MNVDX  = ( MNVDLFX - 1 ) / MOREE2
      DO I=1,MONDLFX
         DMCN( MNVDX + I ) = RELMIN
      ENDDO

C     RECHERCHE DES MIN MAX DES COORDONNEES DE L'OBJET
      NBNOEU = MCN(MNXYZN+WNBNOE)
      NBCOOR = MCN(MNXYZN+WBCOON)
      CALL MIMXPT( NBCOOR, NBNOEU, RMCN(MNXYZN+WYZNOE), COOEXT )


      IF( NDIM .EQ. 6 ) THEN
C        ================================================================
C        TRAITEMENT PARTICULIER DE LA CONDITION DE DIRICHLET D'UN 6-CUBES
C        TOUS LES NOEUDS AYANT UNE COORDONNEE MIN OU MAX SONT
C        SUPPOSES ETRE DES NOEUDS DE DIRICHLET!...
C        LA PRESENCE DE LA FONCTION LU
C        deffonc temperature_bord(t,x,y,z,u,v,w,nutyob,nuobj,temper)
C        IMPOSE LA VALEUR CALCULEE AVEC CETTE FONCTION => NON HOMOGENE
C        L'ABSENCE DE deffonc temperature_bord => CONDITION HOMOGENE
C        ================================================================

C        RECHERCHE DES MIN MAX DES COORDONNEES DE L'OBJET
         NBNOEU = MCN(MNXYZN+WNBNOE)
         NBCOOR = MCN(MNXYZN+WBCOON)
         CALL MIMXPT( NBCOOR, NBNOEU, RMCN(MNXYZN+WYZNOE), COOEXT )

         CALL LXNMNO( NTFONC,'TEMPERATURE_BORD', NOFOTE, I )
         NBDLFX = 0
         MNX    = MNXYZN+WYZNOE-1
         DO 40 J=1,NBNOEU
            DO 30 I=1,NBCOOR
               XYZUVW = RMCN(MNX+I)
               IF( XYZUVW .EQ. COOEXT(I,1) .OR.
     %             XYZUVW .EQ. COOEXT(I,2) ) THEN
C                 CE NOEUD J EST DE DIRICHLET
                  MCN( MNNDLFX + NBDLFX ) = J
                  NBDLFX = NBDLFX + 1
                  IF( NOFOTE .EQ. 0 ) THEN

C                    SI PAS DE FONCTION TEMPERATURE_BORD
C                    LE DEGRE DE LIBERTE AU NOEUD J EST FIXE A 0D0
C                    => CONDITION DE DIRICHLET HOMOGENE
                     DMCN(MNVDX+NBDLFX) = 0D0

                  ELSE
C                    SINON LE DEGRE DE LIBERTE AU NOEUD J EST FIXE
C                    A LA VALEUR EN CE NOEUD DE LA FONCTION
C                    deffonc temperature_bord(t,x,y,z,u,v,w,nutyob,nuobj,temper)
C                    => CONDITION DE DIRICHLET NON HOMOGENE
                     DPARA(1) = 0D0
                     DPARA(2) = RMCN(MNX+1)
                     DPARA(3) = RMCN(MNX+2)
                     DPARA(4) = RMCN(MNX+3)
                     DPARA(5) = RMCN(MNX+4)
                     DPARA(6) = RMCN(MNX+5)
                     DPARA(7) = RMCN(MNX+6)
                     DPARA(8) = 0D0
                     DPARA(9) = 0D0
                     DPARA(10)= 0D0
                     CALL FONVAL( NOFOTE, 10, DPARA,
     %                            NCODEV, DMCN(MNVDX+NBDLFX) )
                  ENDIF
                  GOTO 35
               ENDIF
 30         ENDDO
 35         MNX = MNX + NBCOOR
 40      ENDDO
         GOTO 160
      ENDIF


C     =================================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS 2D ou 3D
C     =================================================
      DO NOTYEL = 1, NBTYEL

C        L'ADRESSE DU TABLEAU NPEF"No TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )

C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF

C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )

C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO NUELEM = 1, NBELEM

C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )

C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, SURFACES ET VOLUMES
            CALL EFTNND( NOOBVC,   NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB,   MNDOEL,
     %                  'CONTACT', MXDOTH, LPCONT,
     %                   NOTYOB )

C           LE RECENSEMENT DES CONTACTS NODAUX
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO J=1,NBNOE

C              LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )

C              LE TYPE OBJET DU NOEUD J DE L'EF
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
               MN1   = NOTYOB(3,J)

               IF( NTYOB .LE. NDIM .AND. NOOB .GT. 0 ) THEN

C                 EXISTE-T-IL UN "CONTACT" EN CE NOEUD ?
C                 --------------------------------------
                  IF( MN1 .GT. 0 ) THEN
C                    CALCUL DE LA TEMPERATURE FIXEE EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XYZP(1) = RMCN(N)
                     XYZP(2) = RMCN(N+1)
                     XYZP(3) = RMCN(N+2)
                     CALL RECONT( NTYOB, NOOB, 3, XYZP,
     %                            MN1,   D )

C                    LA VALEUR DE LA TEMPERATURE FIXEE
                     IF( NCDLF0 .EQ. 0 .AND. ABS(D) .GT. 1D-6 ) THEN

                        NBLGRC(NRERR) = 2
                        IF( LANGAG .EQ. 0 ) THEN
                           KERR(1) = 'TEMPERATURE OU DEPLACEMENT'
                           KERR(2) = 'NON NUL ICI INTERDIT'
                        ELSE
                           KERR(1) = 'TEMPERATURE or DISPLACEMENT'
                           KERR(2) = 'NOT NULL. FORBIDDEN HERE'
                        ENDIF
                        CALL LEREUR
                        IERR = 1
                        GOTO 9999

                     ENDIF

C                    LE DEGRE DE LIBERTE EST FIXE A D
                     DMCN( MNVDX + NONOE ) = D
C                    LE TEMOIN DE VALEUR FIXEE
                     MCN( MNNDX + NONOE ) = 1
                  ENDIF
               ENDIF

            ENDDO
        ENDDO
      ENDDO

C     COMPRESSION DES TABLEAUX DES DL FIXES
C     =====================================
      DLMIN  = -RELMIN
      DLMAX  =  RELMIN
      NBDLFX = 0

      DO I=1,NTDL

C        LE TEMOIN DE VALEUR FIXEE
         IF( MCN( MNNDX + I ) .NE. 0 ) THEN

C           UN DEGRE DE LIBERTE FIXE DE PLUS
            NBDLFX = NBDLFX + 1

C           LE NUMERO DU DL FIXE
            MCN( MNNDX + NBDLFX ) = I

C           LA VALEUR DU DL FIXE
            D = DMCN( MNVDX + I )

C           LE REPOSITIONNEMENT DE LA VALEUR DU DL FIXE
            DMCN( MNVDX + NBDLFX ) = D

            IF( D .LT. DLMIN ) DLMIN = D
            IF( D .GT. DLMAX ) DLMAX = D

         ENDIF

      ENDDO

C     REDUCTION EVENTUELLE DES TABLEAUX DES DEGRES DE LIBERTE FIXES
C     =============================================================
 160  IF( NBDLFX .EQ. 0 ) THEN
C        DESTRUCTION DES 2 TABLEAUX
         CALL TNMCDS( 'ENTIER', MONDLFX, MNNDLFX )
         CALL TNMCDS( 'REEL2',  MONDLFX, MNVDLFX )
         MONDLFX = 0
      ELSE
C        REDUCTION DE LA TAILLE DES 2 TABLEAUX
         MONDLFX = NBDLFX
         CALL TNMCRA('ENTIER', NTDL, MONDLFX, MNNDLFX )
         CALL TNMCRA('REEL2' , NTDL, MONDLFX, MNVDLFX )
      ENDIF

C     AFFICHAGE FINAL
C     ===============
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) NBNOEU,NBDLFX, DLMIN,DLMAX, MNNDLFX,MNVDLFX
      ELSE
         WRITE(IMPRIM,20170) NBNOEU,NBDLFX, DLMIN,DLMAX, MNNDLFX,MNVDLFX
      ENDIF
10170 FORMAT(' thdlfx: TEMPERATURES aux',I9,'NOEUDS.',
     %       ' Nb TEMPERATURES FIXEES=',I9,
     %       ' TEMPERATURE FIXEE MINIMUM=',G14.7,
     %       ' MAXIMUM=',G14.7,' MNNDLFX=',I11,' MNVDLFX=',I11)

20170 FORMAT(' thdlfx: COMPUTED TEMPERATURES at',I9,' NODES.',
     %       ' FIXED TEMPERATURES Nb=',I9,
     %       ' FIXED TEMPERATURE MINIMUM=',G14.7,
     %       ' MAXIMUM=',G14.7,' MNNDLFX=',I11,' MNVDLFX=',I11)

CCC      WRITE(IMPRIM,10180) (MCN(MNNDX+I),DMCN(MNVDX+I),I=1,NBDLFX)
CCC10180 FORMAT(3(' TEMPERATURE',I9,' FIXEE A',G15.7,2X))

      RETURN

C     ERREUR RENCONTREE
C     =================
 9999 IF( MONDLFX .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', MONDLFX, MNNDLFX )
         CALL TNMCDS( 'REEL2',  MONDLFX, MNVDLFX )
         MONDLFX = 0
         NBDLFX  = 0
      ENDIF

      RETURN
      END
