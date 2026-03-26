      SUBROUTINE VITEFX( RELMIN, NTDLVI, NDIM,   MNXYZN,
     %                   NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %                   NBVCFX, NBVIFX, MNNVIFX, MNVVIFX, VIFXMAX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER LES NUMEROS ET VALEURS DES DEGRES DE LIBERTE FIXES
C -----    DE LA VITESSE du FLUIDE
C          LA VITESSE AU BARYCENTRE DES EF N'EST PAS PRISE EN COMPTE
C          POUR LES EF DE BREZZI-FORTIN CAR UNE VALEUR AU BARYCENTRE N'EST
C          PAS SUR LA FRONTIERE DONC N'EST PAS PRISE EN COMPTE
C
C ENTREES:
C --------
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C NTDLVI : NOMBRE TOTAL DE DEGRES DE LIBERTE DES VITESSES
C NDIM   : NOMBRE DE COMPOSANTES DE LA VITESSE=DIMENSION DE L'ESPACE
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NUMIOB : NUMERO MINIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C
C SORTIES:
C --------
C NBVCFX : NOMBRE DE DEGRES DE LIBERTE DE VITESSE CONVECTEE IMPOSEE
C NBVIFX : NOMBRE DE DL VITESSE FIXES
C          ERREUR SI NBVIFX=0 CAR PAS DE CONDITION AUX LIMITES!
C MNNVIFX: ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES VITESSE
C          POUR UN TABLEAU VITESSE(NBNOVI*NDIM)
C MNVVIFX: ADRESSE MCN DU TABLEAU MC DES VALEURS DES DL FIXES
C VIFXMAX: LE MODULE DE LA VITESSE FIXEE MAXIMALE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___blvitesse.inc"
      include"./incl/a___force.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), MNDOEL(4)
      INTEGER           NONOEF(10)
C
      DOUBLE PRECISION  D, RELMIN, XPOIN, YPOIN, ZPOIN, BLVIT(3),
     %                  VIFX, VIFXMAX
C
C     POUR EVITER DE DOUBLER LES TABLEAUX
      IF( NBVIFX .GT. 0 .AND. MNNVIFX .GT. 0 ) THEN
         CALL TNMCDS( 'REEL2',  NBVIFX, MNVVIFX )
         CALL TNMCDS( 'ENTIER', NBVIFX, MNNVIFX )
      ENDIF
C
      MOREE2 = MOTVAR(6)
      NBVCFX = 0
      NBVIFX = 0
      MNNVIFX = 0
      MNVVIFX = 0
      VIFXMAX = 0D0
C
C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      CALL TNMCDC( 'ENTIER', NTDLVI, MNNVIFX )
      MNNDX = MNNVIFX - 1
      CALL AZEROI( NTDLVI, MCN(MNNVIFX) )
C
C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
C     VALEUR DE NON INITIALISATION = RELMIN
      CALL TNMCDC( 'REEL2', NTDLVI, MNVVIFX )
      MNVDX = ( MNVVIFX - 1 ) / MOREE2
      DO I=1,NTDLVI
         DMCN(MNVDX + I) = RELMIN
      ENDDO
C
C     ========================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      NBNOVI = NTDLVI / NDIM
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU ELEMENTS (NPEF)
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
         MNNDEL = MNELE + WUNDEL
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 90 NUELEM = 1, NBELEM
C
C           NO DES NOEUDS DE L'ELEMENT FINI NUELEM
            CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C
C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           TRAITEMENT DES COMPOSANTES FIXEES DE LA VITESSE
C           -----------------------------------------------
C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES
C           SURFACES SUR LE VOLUME
C           LIGNES SUR LES SURFACES ET LE VOLUME
C           POINTS SUR LES LIGNES, LES SURFACES ET LE VOLUME
            CALL EFTNND( NOOBVC,    NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB,    MNDOEL,
     %                  'BLVITESSE', MXDOFL, LPBLVI,
     %                   NOTYOB )
C
C           LE RECENSEMENT DES VITESSES BLOQUEES AUX NOEUDS DE CET EF
C           EXCEPTE POUR LE BARYCENTRE DES EF DE BREZZI-FORTIN
C           QUI NE SONT PAS SUR LA FRONTIERE DU DOMAINE
            DO 40 J=1,NBNOEF
C
C              LE NUMERO DU NOEUD
               NONOE = NONOEF(J)
C
C              LE TYPE OBJET DU NOEUD
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
C
               IF( NTYOB .LE. NDIM .AND. NOOB .GT. 0 ) THEN
C
C                 EXISTE-T-IL UNE FIXATION VITESSE EN CE NOEUD ?
C                 ----------------------------------------------
                  MN1 = NOTYOB(3,J)
                  IF( MN1 .GT. 0 ) THEN
C
C                    CALCUL DES VITESSES FIXES EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XPOIN = RMCN(N)
                     YPOIN = RMCN(N+1)
                     ZPOIN = RMCN(N+2)
                     CALL REBLVI( NTYOB,  NOOB, XPOIN,YPOIN,ZPOIN, MN1,
     %                            NBCOBV, BLVIT )
C
C                    NBCOBV EST LE NOMBRE DE COMPOSANTES FIXEES
                     VIFX = 0D0
                     DO 30 K=1,NBCOBV
C
C                       LE NUMERO DE LA COMPOSANTE K
                        NOCOMP = MCN( MN1 + WUCOBV - 1 + K )

C                       LE NO DU DL EN VITESSE
                        NODLVP = NONOE + (NOCOMP-1) * NBNOVI
C
C                       LE TEMOIN DE FIXATION DU D.L.
                        MCN( MNNDX + NODLVP ) = NOCOMP
C
C                       LA VALEUR DE LA COMPOSANTE DE LA VITESSE FIXEE
                        DMCN( MNVDX + NODLVP ) = BLVIT(K)

                        VIFX = VIFX + BLVIT(K) **2
C
 30                  CONTINUE

                     IF( VIFX .GT. VIFXMAX ) THEN
                        VIFXMAX = VIFX
                     ENDIF

ccc                 if( nonoe .GT. NBNOVI-10 ) then
ccc                     print *,'vitefx: noeud',nonoe,
ccc     %               ' x=',xpoin,' y=',ypoin,' z=',zpoin,
ccc     %               ' V',nocomp,'=',BLVIT(NBCOBV)
ccc                 endif
                  ENDIF
               ENDIF
C
 40         CONTINUE
 90      CONTINUE
100   CONTINUE

C     LA VITESSE FIXEE MAXIMALE
      VIFXMAX = SQRT( VIFXMAX ) 

C     COMPRESSION DU TABLEAU DES DL VITESSE FIXES
C     ===========================================
      NBVCFX=0
      DO I=1,NTDLVI
         IF( MCN( MNNDX + I ) .NE. 0 ) THEN
            NBVIFX = NBVIFX + 1
C           NO DU DL VITESSE FIXE POUR VITESSE(NBNOVI*NDIM)
            MCN( MNNDX + NBVIFX ) = I
C           VALEUR DU DL FIXE
            D = DMCN( MNVDX + I )
            DMCN( MNVDX + NBVIFX ) = D
            IF( D .EQ. 1D222 ) NBVCFX=NBVCFX+1
         ENDIF
      ENDDO
C
C     REDUCTION EVENTUELLE DES TABLEAUX DE DONNEES
C     ============================================
      IF( NBVIFX .EQ. 0 ) THEN
C        DESTRUCTION DES TABLEAUX
         CALL TNMCDS( 'REEL2',  NTDLVI, MNVVIFX )
         CALL TNMCDS( 'ENTIER', NTDLVI, MNNVIFX )
         GOTO 9999
      ELSE
C        REDUCTION DES 2 TABLEAUX
         CALL TNMCRA('REEL2' , NTDLVI, NBVIFX, MNVVIFX )
         CALL TNMCRA('ENTIER', NTDLVI, NBVIFX, MNNVIFX )
      ENDIF
C
cccC     AFFICHAGE FINAL DES 10 DERNIERS DL BLOQUES
cccC     ==========================================
ccc      WRITE(IMPRIM,*)
ccc      MNNDX = MNNVIFX - 1
ccc      MNVDX = ( MNVVIFX - 1 ) / MOREE2
ccc      WRITE(IMPRIM,20170) NBVIFX, NBVC
cccC
ccc      DO I=NBVIFX-10,NBVIFX
cccC        NO DU DL
ccc         NODL = MCN(MNNDX+I)
cccC        NO DU NOEUD
ccc         NOEU = NODL
cccC        NO DE LA COMPOSANTE
ccc         NOCO = 1
ccc 9001    IF( NOEU .GT. NBNOVI ) THEN
ccc            NOCO = NOCO + 1
ccc            NOEU = NOEU - NBNOVI
ccc            GOTO 9001
ccc         ENDIF
cccC        ADRESSE DES COORDONNEES DE CE NOEUD
ccc         N = MNXYZN + WYZNOE + 3 * NOEU - 3
ccc         WRITE(IMPRIM,*)'NOEUD',NOEU,' X=',RMCN(N),' Y=',RMCN(N+1),
ccc     %        ' Z=',RMCN(N+2), ' COMPOSANTE',NOCO,' NoDL=',NODL,
ccc     %        ' FIXE A', DMCN(MNVDX+I)
ccc      ENDDO
cccC
ccc10170 FORMAT(/' NOMBRE DE COMPOSANTES VITESSES  FIXEES =',I10/
ccc     %        ' NOMBRE DE COMPOSANTES VITESSES CONVECTEES FIXEES=',I10)
ccc20170 FORMAT(/' NUMBER of FIXED VELOCITY COMPONENTS =',I10/
ccc     %        ' NUMBER of FIXED CONVECTED VELOCITY COMPONENTS=',I10)
ccc
ccc10180 FORMAT(' VITESSE',I10,' FIXEE A',G15.7,2X)
ccc20180 FORMAT(' VELOCITY COMPONENT',I10,' FIXED to',G15.7,2X)
C
 9999 print *,'vitefx: NDIM=',NDIM,' NBVCFX=',NBVCFX,' NBVIFX=',NBVIFX,
     %        ' MNNVIFX=',MNNVIFX, ' MNVVIFX=', MNVVIFX

      RETURN
      END
