      SUBROUTINE CRLI20( NMLGEX, NUPTVT, OUVFER, NBPTEX, NUPTEX,
     %                   NBLG,   MNNULG, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES ARETES DE L'UNION DES LIGNES EXTRUDEES
C -----    D'UN POLYGONE DEFINI PAR SES NBPTEX SOMMETS ET LES 2 POINTS
C          DE DEFINITION DES EXTREMITES DU VECTEUR TRANSLATION
C
C ENTREES:
C --------
C NMLGEX : NOM    DE LA LIGNE UNION DANS LE LEXIQUE DES LIGNES
C NUPTVT : NUMERO DES 2 POINTS EXTREMITES DU VECTEUR TRANSLATION
C OUVFER : CODE POLYGONE OUVERT(0) ou FERME(1)
C NBPTEX : NOMBRE DE  POINTS SOMMETS DU POLYGONE INITIAL
C NUPTEX : NUMERO DES POINTS SOMMETS DU POLYGONE INITIAL
C
C SORTIES:
C --------
C NBLG   : NOMBRE DE LIGNES DES 2 POLYGONES ET EXTRUDEES
C MNNULG : NUMERO DE LIGNE DANS LE LEXIQUE DES LIGNES DES NBLG LIGNES
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  FEVRIER 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/darete.inc"
C
C     LE SUPER-TABLEAU OU TOUS LES TMC ET TMS OUVERTS SONT STOCKES
C     ICI LE TABLEAU DMCN REEL DOUBLE PRECISION N'EST PAS UTILE
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      ( MCN(1), RMCN(1), DMCN(1) )
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      CHARACTER*24      NOMPOI, NOMLIG, NMLGEX
      CHARACTER*4       KNBRE
      INTEGER           OUVFER, NUPTVT(2), NUPTEX(NBPTEX)
C
      DOUBLE PRECISION  D, XYZTRA(3,2), RAIGEO
      REAL              RXYZ(3)
      INTRINSIC         NINT, SQRT
C
      MOREE2  = MOTVAR(6)
      MNNUPTF = 0
      MNNULG  = 0
      MNPTEX  = 0
C
C     NUPTVT(1:2) nom des 2 points du vecteur translation
      DO K=1,2
         CALL LXNLOU( NTPOIN, NUPTVT(K), NT, MN )
         IF( MN .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'POINT INCONNU POUR LA TRANSLATION'
            ELSE
               KERR(1) = 'UNKNOWN POINT TO DEFINE THE TRANSLATION'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
         CALL LXTSOU( NT, 'XYZSOMMET',  NTSOM, MNSOM )
         IF( MNSOM .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'POINT DE COORDONNEES INCONNUES'
            ELSE
               KERR(1) = 'POINT WITH UNKWOWN COORDINATES'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
C        LES 3 COORDONNEES DU POINT K DU POLYGONE
         MN = MNSOM + WYZSOM - 1
         DO I = 1, 3
            XYZTRA(I,K) = RMCN( MN+I )
         ENDDO
      ENDDO
C
C     NBPTEX nombre de points du polygone
      IF( (OUVFER .EQ. 0 .AND. NBPTEX .LT. 2) .OR.
     %    (OUVFER .NE. 0 .AND. NBPTEX .LT. 3) .OR.
     %     NBPTEX .GE. 100 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE INCORRECT DE POINTS DU POLYGONE'
         ELSE
            KERR(1) = 'INCORRECT NUMBER OF POLYGON POINTS'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     RECUPERATION DES 3 COORDONNEES DES NBPTEX POINTS DU POLYGONE INITIAL
C     --------------------------------------------------------------------
      WRITE(IMPRIM,*)
      CALL TNMCDC( 'REEL2', 3*NBPTEX*2, MNPTEX )
      MNP0 = (MNPTEX-1)/ MOREE2
      DO K=1, NBPTEX
C        LE NUMERO DU POINT K DU POLYGONE INITIAL
         NUMPT = NUPTEX(K)
         CALL LXNLOU( NTPOIN, NUMPT, NT, MN )
         IF( MN .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'POINT INCONNU SUR LE POLYGONE'
            ELSE
               KERR(1) = 'UNKNOWN POINT VERTEX OF THE POLYGON'
            ENDIF
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF
         CALL LXTSOU( NT, 'XYZSOMMET',  NTSOM, MNSOM )
         IF( MNSOM .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'POINT DE COORDONNEES INCONNUES'
            ELSE
               KERR(1) = 'POINT WITH UNKWOWN COORDINATES'
            ENDIF
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
C        LES 3 COORDONNEES DU POINT K DU POLYGONE INITIAL
         MN = MNSOM + WYZSOM - 1
         DO I = 1, 3
            DMCN( MNP0+I ) = RMCN( MN + I )
         ENDDO
C        LE NOM DU POINT K DU POLYGONE
         CALL NMOBNU( 'POINT', NUMPT, NOMPOI )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) K,'INITIAL',NOMPOI,
     %                         (DMCN(MNP0+L),L=1,3)
         ELSE
            WRITE(IMPRIM,20000) 'INITIAL',K,NOMPOI,
     %                         (DMCN(MNP0+L),L=1,3)
         ENDIF
C
C        PASSAGE AU POINT SUIVANT
         MNP0 = MNP0 + 3
      ENDDO
10000 FORMAT('LE POINT',I3,' DU POLYGONE ',A10,' EST ',A,
     %       ' X=',G14.6,' Y=',G14.6,' Z=',G14.6)
20000 FORMAT('THE ',A10,' POLYGON POINT',I3,' IS ',A,
     %       ' X=',G14.6,' Y=',G14.6,' Z=',G14.6)
C
C     LES 3 COORDONNEES DES NBPTEX POINTS DU POLYGONE FINAL PAR TRANSLATION
C     ---------------------------------------------------------------------
      WRITE(IMPRIM,*)
      CALL TNMCDC( 'ENTIER', NBPTEX, MNNUPTF )
      MNP0 = (MNPTEX-1)/ MOREE2
      MNP1 = MNP0 + 3 * NBPTEX
      DO K=1, NBPTEX
C
C        LES 3 COORDONNEES DU POINT K DU POLYGONE FINAL
         DO I = 1, 3
            RXYZ(I) = REAL( DMCN( MNP0+I ) + XYZTRA(I,2)-XYZTRA(I,1) )
         ENDDO
C
C        CE POINT EXISTE T IL DEJA?
         CALL XYZDSPT( RXYZ, NUMPT, RXYZ )
         DO I = 1, 3
            DMCN( MNP1+I ) = RXYZ( I )
         ENDDO
C
         IF( NUMPT .GT. 0 ) THEN
C           LE POINT EXISTE DEJA
C           NUMERO DU POINT DANS LE LX DES POINTS
            MCN( MNNUPTF - 1 + K ) = NUMPT
C           LE NOM DU POINT RETROUVE
            CALL NMOBNU( 'POINT', NUMPT, NOMPOI )
            GOTO 20
         ENDIF
C
C        CONSTRUCTION DU NOUVEAU POINT K DANS LE LX POINTS
C        LE NOM DU POINT K EST CELUI DE LA LIGNE UNION AVEC LE NUMERO K
         NOMPOI = NMLGEX
         N = NUDCNB( NOMPOI )
         IF( K .LT. 10 ) THEN
            WRITE(KNBRE(1:1),'(I1)') K
            N = MIN(N,22)
            NOMPOI = NOMPOI(1:N) // '_' // KNBRE(1:1)
         ELSE
            WRITE(KNBRE(1:2),'(I2)') K
            N = MIN(N,21)
            NOMPOI = NOMPOI(1:N) // '_' // KNBRE(1:2)
         ENDIF
C        SI CE POINT EXISTE, IL EST DETRUIT
         CALL LXLXOU( NTPOIN, NOMPOI, NT1PT, MN1PT )
         IF( MN1PT .GT. 0 ) CALL LXTSDS( NTPOIN, NOMPOI )
C        CONSTRUCTION DU LX POINT
         CALL LXLXDC( NTPOIN, NOMPOI, 24, 8 )
         CALL LXLXOU( NTPOIN, NOMPOI, NT1PT, MN1PT )
C        NUMERO DU POINT DANS LE LX DES POINTS
         CALL NUOBNM( 'POINT', NOMPOI, MCN(MNNUPTF-1+K) )
C
C        CREATION DU TMS DEFINITION AVEC NUTYPO=1: 'COORDONNEES X Y Z'
         CALL LXTNDC( NT1PT, 'DEFINITION', 'MOTS', WOORPO+3 )
         CALL LXTSOU( NT1PT, 'DEFINITION', NT1PDE, MN1PDE )
C        TRANSFORMATION (I POUR IDENTITE)'
         MCN( MN1PDE + WTYTRP ) = 1
C        NUMERO DU TYPE DU POINT
         MCN( MN1PDE + WUTYPO ) = 1
C        COORPO 3 COORDONNEES CARTESIENNES DU POINT XYZ
         RMCN( MN1PDE + WOORPO     ) = REAL( DMCN( MNP1+1 ) )
         RMCN( MN1PDE + WOORPO + 1 ) = REAL( DMCN( MNP1+2 ) )
         RMCN( MN1PDE + WOORPO + 2 ) = REAL( DMCN( MNP1+3 ) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN1PDE + MOTVAR(6) ) = NONMTD( '~>POINT>>DEFINITION' )
C        AJOUT DE LA DATE
         CALL ECDATE( MCN(MN1PDE) )
C
C        CREATION DU TMS XYZSOMMET  POUR SES 3 COORDONNEES X Y Z
         CALL LXTNDC( NT1PT, 'XYZSOMMET', 'MOTS', WYZSOM+3 )
         CALL LXTSOU( NT1PT, 'XYZSOMMET', NT1PXYZ, MN1PXYZ )
C        NBSOM 'Nombre de sommets'
         MCN( MN1PXYZ + WNBSOM ) = 1
C        NBTGS 'Nombre de tangentes'
         MCN( MN1PXYZ + WNBTGS ) = 0
C        NBCOOR 'Nombre coordonnees d'un sommet
         MCN( MN1PXYZ + WBCOOR ) = 3
C        COORPO '3 coordonnees cartesiennes du point' xyz
         RMCN( MN1PXYZ + WYZSOM     ) = REAL( DMCN( MNP1+1 ) )
         RMCN( MN1PXYZ + WYZSOM + 1 ) = REAL( DMCN( MNP1+2 ) )
         RMCN( MN1PXYZ + WYZSOM + 2 ) = REAL( DMCN( MNP1+3 ) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN1PXYZ + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C        AJOUT DE LA DATE
         CALL ECDATE( MCN(MN1PXYZ) )
C
 20      IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) K,'TRANSLATE',NOMPOI,
     %                         (DMCN(MNP1+L),L=1,3)
         ELSE
            WRITE(IMPRIM,20000) 'TRANSLATED',K,NOMPOI,
     %                         (DMCN(MNP1+L),L=1,3)
         ENDIF
C
C        PASSAGE AU SUIVANT
         MNP0 = MNP0 + 3
         MNP1 = MNP1 + 3
      ENDDO
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE' DES ARETES AUTOUR DU POINT
C     ==========================================================================
      NOFOTI = NOFOTIEL()
C     NOFOTI>0 SI CETTE FONCTION 'TAILLE_IDEALE(X,Y,Z)' ou 'EDGE_LENGTH' EXISTE
      NBARLI = 1
      RAIGEO = 1D0
C
C     LES NBPTEX LIGNES DU POLYGONE INITIAL (PLAN NON NECESSAIRE)
C     POUR FERMER LE POLYGONE, LE DERNIER POINT EST LE MEME QUE LE PREMIER
C     ====================================================================
      WRITE(IMPRIM,*)
      CALL TNMCDC( 'ENTIER', 3*NBPTEX, MNNULG )
      NBLG = 0
      MNP0 = (MNPTEX+1)/ MOREE2
      IF( OUVFER .EQ. 0 ) THEN
         NBL = NBPTEX - 1
      ELSE
         NBL = NBPTEX
      ENDIF
      DO K = 1, NBL
C
C        LE NOM DE LA LIGNE EST CELUI DE LA LIGNE UNION AVEC SON NUMERO K
         NOMLIG = NMLGEX
         N = NUDCNB( NOMLIG )
         NBLG = NBLG + 1
         IF( K .LT. 10 ) THEN
            WRITE(KNBRE(1:1),'(I1)') K
            N = MIN(N,22)
            NOMLIG = NOMLIG(1:N) // '_' // KNBRE(1:1)
         ELSE
            WRITE(KNBRE(1:2),'(I2)') K
            N = MIN(N,21)
            NOMLIG = NOMLIG(1:N) // '_' // KNBRE(1:2)
         ENDIF
C        SI CETTE LIGNE EXISTE, ELLE EST DETRUITE
         CALL LXLXOU( NTLIGN, NOMLIG, NT1LG, MN1LG )
         IF( MN1LG .GT. 0 ) CALL LXTSDS( NTLIGN, NOMLIG )
C
C        CONSTRUCTION DE LA LIGNE K
         CALL LXLXDC( NTLIGN, NOMLIG, 24, 8 )
         CALL LXLXOU( NTLIGN, NOMLIG, NT1LG, MN1LG )
C        NUMERO DE LA LIGNE DANS LE LX DES LIGNES
         CALL NUOBNM( 'LIGNE', NOMLIG, MCN(MNNULG-1+NBLG) )
C
C        LA LIGNE A UN TMS DEFINITION DE LIGNE SEGMENT DE DROITE
         CALL LXTNDC( NT1LG, 'DEFINITION', 'MOTS', WUPTFI+1 )
         CALL LXTSOU( NT1LG, 'DEFINITION', NT1LDE, MN1LDE )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MN1LDE + WTYTRL ) = 1
C        TYPE DE LA LIGNE
         MCN( MN1LDE + WUTYLI ) = 2
C        LES 2 POINTS EXTREMITES
         NUPT1 = NUPTEX(K)
         IF( OUVFER .NE. 0 .AND. K .EQ. NBL ) THEN
C           POLYGONE FERME   LIGNE DERNIER POINT -> 1-ER POINT
            NUPT2 = NUPTEX(1)
            MNP1 = (MNPTEX+1)/ MOREE2
         ELSE
            NUPT2 = NUPTEX(K+1)
            MNP1 = MNP0 + 3
         ENDIF
C        LA DEFINITION DE LA LIGNE K DU POLYGONE INITIAL
         MCN( MN1LDE + WBARLI ) = NBARLI
        RMCN( MN1LDE + WAIGEO ) = REAL( RAIGEO )
         MCN( MN1LDE + WUPTIN ) = NUPT1
         MCN( MN1LDE + WUPTFI ) = NUPT2
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN1LDE + MOTVAR(6) ) = NONMTD( '~>LIGNE>>DEFINITION' )
C        AJOUT DE LA DATE
         CALL ECDATE( MCN(MN1LDE) )
C
C        CONSTRUCTION DES TMS NSEF et XYZSOMMET DE LA LIGNE K
         CALL SEGDRO( NT1LG,  DMCN(MNP0), DMCN(MNP1),
     %                NBARLI, RAIGEO,
     %                NTNS1L, MNNS1L, NTXY1L, MNXY1L, IERR )
         IF( IERR .GT. 0 ) GOTO 9999
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10040) K, 'INITIAL', NOMLIG
         ELSE
            WRITE(IMPRIM,20040) 'INITIAL', K, NOMLIG
         ENDIF
         MNP0 = MNP0 + 3
C
      ENDDO
10040 FORMAT('LA LIGNE',I4,' DU POLYGONE ',A10,' EST ',A)
20040 FORMAT('THE ',A10,' POLYGON LINE',I4,' IS ',A)
C
C     LES NBPTEX LIGNES DU POLYGONE FINAL
C     ===================================
      WRITE(IMPRIM,*)
      MNP0 = (MNPTEX+1)/ MOREE2 + 3 * NBPTEX
      DO K = 1, NBL
C
C        LE NOM DE LA LIGNE EST CELUI DE LA LIGNE UNION AVEC SON NUMERO K
         NOMLIG = NMLGEX
         N = NUDCNB( NOMLIG )
         NBLG = NBLG + 1
         IF( NBLG .LT. 10 ) THEN
            WRITE(KNBRE(1:1),'(I1)') NBLG
            N = MIN(N,22)
            NOMLIG = NOMLIG(1:N) // '_' // KNBRE(1:1)
         ELSE
            WRITE(KNBRE(1:2),'(I2)') NBLG
            N = MIN(N,21)
            NOMLIG = NOMLIG(1:N) // '_' // KNBRE(1:2)
         ENDIF
C        SI CETTE LIGNE EXISTE, ELLE EST DETRUITE
         CALL LXLXOU( NTLIGN, NOMLIG, NT1LG, MN1LG )
         IF( MN1LG .GT. 0 ) CALL LXTSDS( NTLIGN, NOMLIG )
C
C        CONSTRUCTION DE LA LIGNE K
         CALL LXLXDC( NTLIGN, NOMLIG, 24, 8 )
         CALL LXLXOU( NTLIGN, NOMLIG, NT1LG, MN1LG )
C        NUMERO DE LA LIGNE DANS LE LX DES LIGNES
         CALL NUOBNM( 'LIGNE', NOMLIG, MCN(MNNULG-1+NBLG) )
C
C        LA LIGNE A UN TMS DEFINITION DE LIGNE SEGMENT DE DROITE
         CALL LXTNDC( NT1LG, 'DEFINITION', 'MOTS', WUPTFI+1 )
         CALL LXTSOU( NT1LG, 'DEFINITION', NT1LDE, MN1LDE )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MN1LDE + WTYTRL ) = 1
C        TYPE DE LA LIGNE
         MCN( MN1LDE + WUTYLI ) = 2
C        LES 2 POINTS EXTREMITES DE LA LIGNE
         NUPT1 = MCN( MNNUPTF-1+K )
         IF( OUVFER .NE. 0 .AND. K .EQ. NBL ) THEN
C           POLYGONE FERME   LIGNE DERNIER POINT -> 1-ER POINT
            NUPT2 = MCN( MNNUPTF )
            MNP1 = (MNPTEX+1)/ MOREE2 + NBPTEX*3
         ELSE
            NUPT2 = MCN( MNNUPTF+K )
            MNP1 = MNP0 + 3
         ENDIF
C        LA DEFINITION DE LA LIGNE K DU POLYGONE FINAL
         MCN( MN1LDE + WBARLI ) = NBARLI
        RMCN( MN1LDE + WAIGEO ) = REAL( RAIGEO )
         MCN( MN1LDE + WUPTIN ) = NUPT1
         MCN( MN1LDE + WUPTFI ) = NUPT2
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN1LDE + MOTVAR(6) ) = NONMTD( '~>LIGNE>>DEFINITION' )
C        AJOUT DE LA DATE
         CALL ECDATE( MCN(MN1LDE) )
C
C        CONSTRUCTION DES TMS NSEF et XYZSOMMET DE LA LIGNE NBLG
         CALL SEGDRO( NT1LG,  DMCN(MNP0), DMCN(MNP1),
     %                NBARLI, RAIGEO,
     %                NTNS1L, MNNS1L, NTXY1L, MNXY1L, IERR )
         IF( IERR .GT. 0 ) GOTO 9999
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10040) K, 'TRANSLATE', NOMLIG
         ELSE
            WRITE(IMPRIM,20040) 'TRANSLATED', K, NOMLIG
         ENDIF
C
         MNP0 = MNP0 + 3
      ENDDO
C     LES NBPTEX LIGNES EXTRUDEES DES SOMMETS DU POLYGONE
C     ===================================================
      WRITE(IMPRIM,*)
C
C     CALCUL DE LA DISTANCE ENTRE LES 2 POINTS DU VECTEUR TRANSLATION
      D = SQRT(  (XYZTRA(1,2)-XYZTRA(1,1))**2
     %         + (XYZTRA(2,2)-XYZTRA(2,1))**2
     %         + (XYZTRA(3,2)-XYZTRA(3,1))**2 )
      IF( D .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LA TRANSLATION EST REDUITE A UN POINT'
         ELSE
            KERR(1) = 'TRANSLATION is REDUCED to a POINT'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     NOMBRE D'ARETES DU VECTEUR TRANSLATION
      IF( DARETE .GT. 0D0 ) THEN
         NBARLI = NINT( D / DARETE )
         IF( NBARLI .LE. 0 ) NBARLI = 1
      ENDIF
C
      MNP0 = (MNPTEX+1)/ MOREE2
      MNP1 =  MNP0 + 3 * NBPTEX
      DO K = 1, NBPTEX
C
C        LE NOM DE LA LIGNE EST CELUI DE LA LIGNE UNION AVEC SON NUMERO K
         NOMLIG = NMLGEX
         N    = NUDCNB( NOMLIG )
         NBLG = NBLG + 1
         IF( NBLG .LT. 10 ) THEN
            WRITE(KNBRE(1:1),'(I1)') NBLG
            N = MIN(N,22)
            NOMLIG = NOMLIG(1:N) // '_' // KNBRE(1:1)
         ELSE
            WRITE(KNBRE(1:2),'(I2)') NBLG
            N = MIN(N,21)
            NOMLIG = NOMLIG(1:N) // '_' // KNBRE(1:2)
         ENDIF
C        SI CETTE LIGNE EXISTE, ELLE EST DETRUITE
         CALL LXLXOU( NTLIGN, NOMLIG, NT1LG, MN1LG )
         IF( MN1LG .GT. 0 ) CALL LXTSDS( NTLIGN, NOMLIG )
C
C        CONSTRUCTION DE LA LIGNE K EXTRUDEE
         CALL LXLXDC( NTLIGN, NOMLIG, 24, 8 )
         CALL LXLXOU( NTLIGN, NOMLIG, NT1LG, MN1LG )
C        NUMERO DE LA LIGNE DANS LE LX DES LIGNES
         CALL NUOBNM( 'LIGNE', NOMLIG, MCN(MNNULG-1+NBLG) )
C
C        LA LIGNE A UN TMS DEFINITION DE LIGNE SEGMENT DE DROITE
         CALL LXTNDC( NT1LG, 'DEFINITION', 'MOTS', WUPTFI+1 )
         CALL LXTSOU( NT1LG, 'DEFINITION', NT1LDE, MN1LDE )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MN1LDE + WTYTRL ) = 1
C        TYPE DE LA LIGNE
         MCN( MN1LDE + WUTYLI ) = 2
C        LES 2 POINTS EXTREMITES
         NUPT1 = NUPTEX(K)
         NUPT2 = MCN( MNNUPTF -1 + K )
C        LA DEFINITION DE LA LIGNE K
         MCN( MN1LDE + WBARLI ) = NBARLI
        RMCN( MN1LDE + WAIGEO ) = REAL( RAIGEO )
         MCN( MN1LDE + WUPTIN ) = NUPT1
         MCN( MN1LDE + WUPTFI ) = NUPT2
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN1LDE + MOTVAR(6) ) = NONMTD( '~>LIGNE>>DEFINITION' )
C        AJOUT DE LA DATE
         CALL ECDATE( MCN(MN1LDE) )
C
C        CONSTRUCTION DES TMS NSEF et XYZSOMMET DE LA LIGNE NBLG
         CALL SEGDRO( NT1LG,  DMCN(MNP0), DMCN(MNP1),
     %                NBARLI, RAIGEO,
     %                NTNS1L, MNNS1L, NTXY1L, MNXY1L, IERR )
         IF( IERR .GT. 0 ) GOTO 9999
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10040) K, 'EXTRUDE', NOMLIG
         ELSE
            WRITE(IMPRIM,20040) 'EXTRUDED', K, NOMLIG
         ENDIF
C
         MNP0 = MNP0 + 3
         MNP1 = MNP1 + 3
      ENDDO
C
C     DESTRUCTION DES TABLEAUX MC DEVENUS INUTILES
 9999 IF( MNNUPTF .GT. 0 ) CALL TNMCDS( 'ENTIER',  NBPTEX,  MNNUPTF )
      IF( MNPTEX  .GT. 0 ) CALL TNMCDS( 'REEL2', 3*NBPTEX*2, MNPTEX )
C
      RETURN
      END
