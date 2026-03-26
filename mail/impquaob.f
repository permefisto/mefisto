      SUBROUTINE IMPQUAOB( KNOMOB, QUAMIN, NBEFMQ, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER L'HISTOGRAMME DES QUALITES DES EF
C -----    LA QUALITE MOYENNE, MINIMALE ET L'ECART TYPE A 1
C          DES QUALITES DES EF DU MAILLAGE D'UN OBJET
C          C-A-D DES VOLUMES  POUR UN OBJET 3D
C                DES SURFACES POUR UN OBJET 2D

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET

C SORTIES:
C --------
C QUAMIN : QUALITE MINIMALE D'UN EF DU MAILLAGE DE L'OBJET
C NBEFMQ : NOMBRE D'EF DE QUALITE INFERIEURE A LA QUALITE MINIMALE
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC Paris & St PIERRE du PERRAY Aout 2012
C....................................................................012
      PARAMETER        (MXTYEL=7, QUAMED=0.005)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/donthe.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"

      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))

      CHARACTER*(*)     KNOMOB
      CHARACTER*128     KNOM
      CHARACTER*36      TITRE
      INTEGER           NONOEF(27)
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
      CHARACTER*4       NOMELE(2)
      INTEGER           NBTYEF(0:9)
      INTEGER           NBEFQU(-3:10)
      REAL              XYZEFS(6,2), XYZARE(3,8), XYZBAR(3)
      INTEGER           NBSOFA(6), NOSOFA(4,6)

      NBEFMQ = 0
      QUAMIN = 0
      MNNBST = 0
      IERR   = 0

C     VERIFICATION DU NOM_DE_L'OBJET et AFFICHAGE
C     VERIFICATION of the OBJECT NAME and DISPLAY
C     ===========================================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
C     FIND the OBJECT NAME in the LEXICON of OBJECTS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF

C     RECHERCHE DU TABLEAU DEFINITION DE L'OBJET
C     FIND the TMS DEFINITION of the OBJECT
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )

C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
C     IF IT DOES NOT EXIST, RETURN TO ASK THE OBJECT NAME
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'SANS DEFINITION => A REDEFINIR'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'WITHOUT DEFINITION => DEFINE AGAIN'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

C     RECHERCHE DU TABLEAU TOPOLOGIE DE L'OBJET
C     FIND the TMS TOPOLOGIE of the OBJECT
      CALL LXTSOU( NTLXOB, 'TOPOLOGIE', NTTOPO, MNTOPO )

C     S'IL N'EXISTE PAS D'INTERPOLATION ALORS
C          RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTTOPO .LE. 0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'OBJET SANS INTERPOLATION'
            KERR(3) = 'RE-EXECUTER MAILLER POUR LE FAIRE'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'OBJECT WITHOUT INTERPOLATION'
            KERR(3) = 'EXECUTE AGAIN MAILLER TO DO IT'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

C     L'ANCIEN HISTORIQUE EST EFFACE
C     PREVIOUS HISTORIC is CANCELLED
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO

C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF" DE l'OBJET
C     FIND THE TMS XYZSOMMET XYZNOEUD XYZPOINT NPEF"xxxx  of the OBJECT
C     Cf $MEFISTO/td/da/a___xyznoeud   a___npef
      NBJEUX = 1
      CALL MIMAOB( NBJEUX, NTLXOB, MXDOTH,
     %             NTTOPO, MNTOPO, NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF, NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) RETURN
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES PLSV
C         NUMAOB          LES 4 NUMEROS MAXIMA DES PLSV
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX MIN MAX DES PLSV
C     MNDOEL THE 4 ADRESSES MCN of THE MIN MAX ARRAY MCN ADDRESSES OF PLSV
C     TABLEAUX DECRIVANT LA THERMIQUE DE L'OBJET COMPLET
C
C     NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
      IF( NBCOOR .NE. 3 ) THEN
         IERR = 2
         RETURN
      ENDIF

C     NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN( MNXYZN + WNBNOE )

C     NDIM LA DIMENSION 1 OU 2 OU 3 OU 6 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBNOEU, MCN(MNXYZN+WYZNOE), NDIM )

C     CALCUL DU NOMBRE DE SOMMETS
      MXNBST = 1+NBNOEU
      CALL TNMCDC( 'ENTIER', MXNBST, MNNBST )
      CALL AZEROI(  MXNBST, MCN(MNNBST) )

C     MISE A ZERO DU NOMBRE D'EF PAR TYPE GEOMETRIQUE
      CALL AZEROI( 10, NBTYEF )

C     INITIALISATIONS POUR LES QUALITES
      DO I=-3,10
         NBEFQU(I) = 0
      ENDDO
      SURVOLEF = 0.
      QUAMOY = 0.0
      QUAMIN = 2.0
      ECATYP = 0.0
      QUAMOB = 0.0
      QUAMIB = 2.0
      ECATYB = 0.0
      NBEFOB = 0
      NBTGS  = 0
      NTYMIN = 1

C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
      DO NOTYEL = 1 , NBTYEL

C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )

C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )

C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )

C        NOMBRE DE SOMMETS DE L'EF
         NBSOEF = NBSOME( NCOGEL )

C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )

C        NOMBRE TOTAL D'EF
         NBEFOB = NBEFOB + NBELEM

C        NOMBRE D'EF SELON LE CODE GEOMETRIQUE
         NBTYEF(NCOGEL) = NBTYEF(NCOGEL) + NBELEM

C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF DE CE TYPE
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
C           POINTS DIFFERENTS DES NOEUDS
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF

C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )

         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'impquaob:',NBELEM,' EF de TYPE',NUTYEL,': ',NOMELE
         ELSE
            PRINT*,'impquaob:',NBELEM,' FE of TYPE',NUTYEL,': ',NOMELE
         ENDIF

C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        DU TABLEAU DES POLY(POINTS INTEGRATION), ...
C        ----------------------------------------------
C        SELON LE TYPE DE L'ELEMENT FINI
         GOTO( 3,3,3,3,1, 1,1,1,1,1,
     %         1,1,3,1,3, 3,1,3,3,3,
     %         3,3,3,3,1, 1,1,3,3,3,
     %         3,3,3,1,1), NUTYEL

C        ERREUR
 1       NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR impquaob: TYPE EF '// NOMELE(1)
     %           // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'ERROR impquaob: FE TYPE '// NOMELE(1)
     %           // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999

 3      DO NUELEM = 1, NBELEM

C           NONOEF(NBNDEL) NUMERO DES NOEUDS DE L'EF NUELEM
C           -----------------------------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, NONOEF )

C           LA QUALITE DE L'EF
C           ------------------
            CALL QUALEF( NCOGEL,   NONOEF, NBNOEU, RMCN(MNXYZN+WYZNOE),
     %                   SURFVOLU, QUALIT, IERR )
            SURVOLEF = SURVOLEF + SURFVOLU

            IF( LORBITE .EQ. 1 ) THEN
               IF( IERR .NE. 0 .OR. QUALIT .LT. QUAMED ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,*) 'impquaob: EF Numero=',NUELEM,
     %                    ' QUALITE=',QUALIT
                     WRITE(IMPRIM,*) '          EF NOEUDS=',
     %                    (NONOEF(I),I=1,NBNDEL)
                  ELSE
                     WRITE(IMPRIM,*) 'impquaob: FE NUMBER=',NUELEM,
     %                    ' QUALITY=',QUALIT
                     WRITE(IMPRIM,*) '          FE NODES =',
     %                    (NONOEF(I),I=1,NBNDEL)
                  ENDIF
               ENDIF
            ENDIF

C           COMPTEUR DU NOMBRE D'EF DE QUALITE MEDIOCRE
            IF( QUALIT .LT. QUAMED ) NBEFMQ = NBEFMQ + 1

C           LA QUALITE MOYENNE
            QUAMOY = QUAMOY + MAX( 0.0, QUALIT )
            ECATYP = ECATYP + (1.0-QUALIT)**2

C           LA QUALITE MINIMALE
            IF( QUALIT .LT. QUAMIN ) THEN
               QUAMIN = QUALIT
               NEFMIN = NUELEM
               NTYMIN = NOTYEL
            ENDIF

C           CONTRIBUTION A L'HISTOGRAMME
            I = INT( QUALIT * 10.0 )
            IF( I .EQ. 0 ) THEN
               IF( QUALIT .GE. 0.01 ) THEN
                  I = 0
               ELSE IF( QUALIT .GE. 0.001 ) THEN
                  I = -1
               ELSE IF( QUALIT .GE. 0.0 ) THEN
                  I = -2
               ELSE
                  I = -3
               ENDIF
            ENDIF
            NBEFQU(I) = NBEFQU(I) + 1

C           DETECTION DU NOMBRE DES SOMMETS
            DO I=1,NBSOEF

C              NUMERO DU SOMMET I DE L'EF NUELEM
               NOST = NONOEF(I)
               IF( NOST .LE. 0 .OR. NOST .GT. NBNOEU ) THEN
                  PRINT*,'impquaob: PB EF',NUELEM,' St:',NOST,
     %                   ' INCORRECT car <0 ou >NBNOEU=',NBNOEU
                  PRINT*,'impquaob: St:',(NONOEF(K),K=1,NBSOEF)
                  NOST = 0
               ENDIF

               MCN(MNNBST+NOST) = NOST
            ENDDO

         ENDDO
      ENDDO

C     LA QUALITE MOYENNE
      QUAMOY = QUAMOY / NBEFOB

C     L'ECART TYPE A 1
      ECATYP = SQRT( ECATYP / NBEFOB )

      NBSOM = 0
      DO NOST=1,NBNOEU
         IF( MCN(MNNBST+NOST) .GT. 0 ) THEN
            NBSOM = NBSOM + 1

ccc         ELSE
ccc           NOST N'EST PAS UN SOMMET DES EF
ccc           PRINT*,'impquaob: NOEUD',NOST,' N''EST PAS un SOMMET=',
ccc                   MCN(MNNBST+NOST),' NBSOM=',NBSOM

         ENDIF
      ENDDO

C     IMPRESSIONS FINALES
C     ===================
      TITRE = 'OBJET'
      I     = NUDCNB( TITRE  )
      N     = NUDCNB( KNOMOB )
      TITRE = TITRE(1:I) // ' ' // KNOMOB(1:N)
      IF( LORBITE .EQ. 1 ) THEN
      IF( LANGAG  .EQ. 0 ) THEN
      WRITE(IMPRIM,10010) TITRE,NBSOM,NBNOEU,NBTGS
      WRITE(IMPRIM,10030) NBEFOB,0,QUAMOY,QUAMIN,NEFMIN,ECATYP

10010 FORMAT(/' QUALITE DU MAILLAGE : ',A,/
     %        ' ---------------------'/
     %' NOMBRE de SOMMETS des EF=',I8,/,
     %' NOMBRE de NOEUDS  des EF=',I8,/,
     %' NOMBRE de TANGENTES     =',I8)
10030 FORMAT(
     %' NOMBRE d''Elements Finis =',I8/,
     %' NOMBRE d''EF AVEC des TGS=',I8/,
     %' QUALITE MOYENNE  des EF =',G10.3,/,
     %' QUALITE MINIMALE des EF =',G10.3,' pour l''EF ',I8,/,
     %' ECART TYPE / 1   des EF =',G10.3)
      IF( NBTYEF(3) .GT. 0 ) WRITE(IMPRIM,10013) NBTYEF(3)
10013 FORMAT(' NOMBRE de TRIANGLES  =',I8)
      IF( NBTYEF(4) .GT. 0 ) WRITE(IMPRIM,10014) NBTYEF(4)
10014 FORMAT(' NOMBRE de QUADRANGLES=',I8)
      IF( NBTYEF(5) .GT. 0 ) WRITE(IMPRIM,10015) NBTYEF(5)
10015 FORMAT(' NOMBRE de TETRAEDRES =',I8)
      IF( NBTYEF(6) .GT. 0 ) WRITE(IMPRIM,10016) NBTYEF(6)
10016 FORMAT(' NOMBRE de PENTAEDRES =',I8)
      IF( NBTYEF(9) .GT. 0 ) WRITE(IMPRIM,10019) NBTYEF(9)
10019 FORMAT(' NOMBRE de PYRAMIDES  =',I8)
      IF( NBTYEF(7) .GT. 0 ) WRITE(IMPRIM,10017) NBTYEF(7)
10017 FORMAT(' NOMBRE de HEXAEDRES  =',I8)
      ELSE
      WRITE(IMPRIM,11010) TITRE,NBSOM,NBNOEU,NBTGS
      WRITE(IMPRIM,11030) NBEFOB,0,QUAMOY,QUAMIN,NEFMIN,ECATYP

11010 FORMAT(/' QUALITY of the MESH : ',A,/
     %        ' ---------------------'/
     %' NUMBER of FE VERTICES     =',I8,/,
     %' NUMBER of FE NODES        =',I8,/,
     %' NUMBER of FE TANGENTS     =',I8)
11030 FORMAT(
     %' NUMBER of Finite Elements =',I8/,
     %' NUMBER of FE WITH TANGENTS=',I8/,
     %' AVERAGE QUALITY of FE     =',G10.3,/,
     %' MINIMUM QUALITY of FE     =',G10.3,' for the FE ',I8,/,
     %' STANDARD DEVIATION/1 of FE=',G12.3)
      IF( NBTYEF(3) .GT. 0 ) WRITE(IMPRIM,11013) NBTYEF(3)
11013 FORMAT(' NUMBER of TRIANGLES  =',I8)
      IF( NBTYEF(4) .GT. 0 ) WRITE(IMPRIM,11014) NBTYEF(4)
11014 FORMAT(' NUMBER of QUADRANGLES=',I8)
      IF( NBTYEF(5) .GT. 0 ) WRITE(IMPRIM,11015) NBTYEF(5)
11015 FORMAT(' NUMBER of TETRAHEDRA =',I8)
      IF( NBTYEF(6) .GT. 0 ) WRITE(IMPRIM,11016) NBTYEF(6)
11016 FORMAT(' NUMBER of PENTAHEDRA =',I8)
      IF( NBTYEF(9) .GT. 0 ) WRITE(IMPRIM,11019) NBTYEF(9)
11019 FORMAT(' NUMBER of PYRAMIDS   =',I8)
      IF( NBTYEF(7) .GT. 0 ) WRITE(IMPRIM,11017) NBTYEF(7)
11017 FORMAT(' NUMBER of HEXAHEDRA  =',I8)
      ENDIF
      ENDIF

      QUALIB = 1.0
      NBEFQU(9) = NBEFQU(9) + NBEFQU(10)

C     TRACE DE LA LEGENDE DES QUALITES SUR ECRAN ET POSTSCRIPT
C     --------------------------------------------------------
      NNX =  35
      NNY = 120
      IF( INTERA .GE. 1 ) THEN

C        EFFACEMENT DE LA LEGENDE DE LA QUALITE SUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
           IF ( LASOPS.EQ.1 ) THEN
             LASOPS = -11
           ELSE
             IF ( LASOPS.EQ.2 ) THEN
               LASOPS = -12
             ELSE
               LASOPS = 0
               NBLGRC(NRERR) = 2
               KERR(1) = 'T1SOBJ: MAUVAISE VALEUR DE LASOPS'
               KERR(2) = '        ARRET DU POSTSCRIPT'
               CALL LEREUR
             ENDIF
           ENDIF
           CALL XVPOSTSCRIPT(LASOPS)
           LASOPS = - LASOPS
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
C        PLACEMENT DE LA LEGENDE SUR L'ECRAN
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = '  QUALITE DES ELEMENTS FINIS'
         ELSE
            KNOM = '  QUALITY of FINITE ELEMENTS'
         ENDIF
         CALL XVCOULEUR( NCROUG )
         CALL XVTEXTE( KNOM(1:28), 28, NNX, NNY )
         NNY = NNY + 16
      ENDIF

      DO 20 N=9,-3,-1
C
C        L'INTERVALLE N DES QUALITES
         QUALIT = N/10.0 + 1E-6
         NN = 0
         IF( N .EQ. 0 ) THEN
            QUALIT = 0.01
            NN = 1
         ELSE IF( N .EQ. -1 ) THEN
            QUALIT = 0.001
            NN = 2
         ELSE IF( N .EQ. -2 ) THEN
            QUALIT = 0.0
            NN = 2
         ELSE IF( N .EQ. -3 ) THEN
            QUALIT = -1.0
            NN = 3
         ENDIF

C        LE NOMBRE D'EF POUR CETTE QUALITE
         IF( NBEFQU(N) .LE. 0 ) GOTO 19

C        LE POURCENTAGE ARRONDI D'EF POUR CETTE QUALITE
         I = NINT( ( NBEFQU(N) * 100.0 ) / NBEFOB )

         IF( LORBITE .EQ. 1 ) THEN
         IF( LANGAG  .EQ. 0 ) THEN
         IF((I.GE.2).OR.(NBEFQU(N).EQ.0)) THEN
           WRITE(IMPRIM,10020) QUALIT,QUALIB,NBEFQU(N),I,('*',J=1,I/2)
         ELSE IF(I.GT.0) THEN
           WRITE(IMPRIM,10020) QUALIT,QUALIB,NBEFQU(N),I,'+'
         ELSE
           WRITE(IMPRIM,10020) QUALIT,QUALIB,NBEFQU(N),I,'.'
         ENDIF
         ELSE
         IF((I.GE.2).OR.(NBEFQU(N).EQ.0)) THEN
           WRITE(IMPRIM,11020) QUALIT,QUALIB,NBEFQU(N),I,('*',J=1,I/2)
         ELSE IF(I.GT.0) THEN
           WRITE(IMPRIM,11020) QUALIT,QUALIB,NBEFQU(N),I,'+'
         ELSE
           WRITE(IMPRIM,11020) QUALIT,QUALIB,NBEFQU(N),I,'.'
         ENDIF
         ENDIF
         ENDIF
10020 FORMAT(F6.3,'=< QUALITE <',F5.3,' :',I6,' EF soit ',
     %       I3,' % : ',100(A1))
11020 FORMAT(F6.3,'=< QUALITY <',F5.3,' :',I6,' FE ',
     %       I3,' % : ',100(A1))
C
         IF( INTERA .GE. 1 .AND. LCRITR .EQ. 1 ) THEN
C           TRACE EN HAUT ET A GAUCHE
            CALL XVCOULEUR( NCOQUA(QUALIT) )
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(KNOM,10021) NBEFQU(N),I,QUALIT,QUALIB
            ELSE
               WRITE(KNOM,11021) NBEFQU(N),I,QUALIT,QUALIB
            ENDIF
            IF( NN .LE. 1 ) KNOM(19:19) = ' '
            IF( NN .LE. 0 ) KNOM(18:18) = ' '
            CALL XVTEXTE( KNOM(1:28+NN), 28+NN, NNX, NNY )
            NNY = NNY + 16
         ENDIF
10021 FORMAT(I6,' EF',I3,' %',F5.3,'< Q <=',F5.3)
10022 FORMAT('  QUALITE MOYENNE =',G10.3)
10023 FORMAT('  QUALITE MINIMALE=',G10.3)
10024 FORMAT('  ECART TYPE A 1  =',G10.3)
10025 FORMAT('  NOMBRE SOMMETS  =',I9)
10026 FORMAT('  NOMBRE TANGENTES=',I9)
10027 FORMAT('  NOMBRE TOTAL EF =',I9)
C
11021 FORMAT(I6,' FE',I3,' %',F5.3,'< Q <=',F5.3)
11022 FORMAT('  FE''s AVERAGE QUALITY=',G10.3)
11023 FORMAT('  FE''s MINIMUM QUALITY=',G10.3)
11024 FORMAT('  ECART TYPE from 1   =',G10.3)
11025 FORMAT('  NUMBER of VERTICES  =',I9)
11026 FORMAT('  NUMBER of TANGENTS  =',I9)
11027 FORMAT('  NUMBER of FE        =',I9)

 19      QUALIB = QUALIT
 20   CONTINUE

C     TRACE DU TEXTE EN GRIS FONCE
      CALL XVCOULEUR( NCGRIS )
      IF( INTERA .GE. 1 .AND. LCRITR .EQ. 1 ) THEN
C        TRACE DE LA QUALITE MOYENNE DU MAILLAGE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10022) QUAMOY
         ELSE
            WRITE(KNOM,11022) QUAMOY
         ENDIF
         NBC = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DE LA QUALITE MINIMALE DU MAILLAGE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10023) QUAMIN
         ELSE
            WRITE(KNOM,11023) QUAMIN
         ENDIF
         NBC = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DE L'ECART TYPE A 1
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10024) ECATYP
         ELSE
            WRITE(KNOM,11024) ECATYP
         ENDIF
         NBC = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 25

C        TRACE DU NOMBRE TOTAL D'EF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10027) NBEFOB
         ELSE
            WRITE(KNOM,11027) NBEFOB
         ENDIF
         CALL SANSDBL( KNOM, NBC )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DU NOMBRE DE SOMMETS
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10025) NBSOM
         ELSE
            WRITE(KNOM,11025) NBSOM
         ENDIF
         CALL SANSDBL( KNOM, NBC )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DU NOMBRE DE TANGENTES
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10026) NBTGS
         ELSE
            WRITE(KNOM,11026) NBTGS
         ENDIF
         CALL SANSDBL( KNOM, NBC )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )

C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
           LASOPS = LASOPS -10
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
      ENDIF
C
C     AFFICHAGE DES COORDONNEES DE L'EF DE QUALITE MINIMALE
C     -----------------------------------------------------
      IF( QUAMIN .LT. 0.3 ) THEN

C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NTYMIN )

C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )

C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )

C        NOMBRE DE SOMMETS DE L'EF
         NBSOEF = NBSOME(NCOGEL)

C        NONOEF(NBNDEL) NUMERO DES NOEUDS DE L'EF NUELEM
         CALL EFNOEU( MNELE, NEFMIN, NBNDEL, NONOEF )

         IF( LORBITE .EQ. 1 ) THEN
            IF( LANGAG  .EQ. 0 ) THEN
               WRITE(IMPRIM,10040) NEFMIN, QUAMIN
10040          FORMAT(/' EF',I8,' de QUALITE MINIMALE ',G10.3)
            ELSE
               WRITE(IMPRIM,11040) NEFMIN, QUAMIN
11040          FORMAT(/' FE',I8,' of MINIMUM QUALITY ',G10.3)
            ENDIF

            DO 50 J=1,NBNDEL
C              LE NUMERO DU NOEUD J DE L'EF
               N = NONOEF(J)
               IF( N .LE. 0 ) GOTO 50
               MN = MNXYZN + WYZNOE + NBCOOR * (N-1)
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10050) N,(RMCN(MN+I),I=0,2)
               ELSE
                  WRITE(IMPRIM,11050) N,(RMCN(MN+I),I=0,2)
               ENDIF
 50         CONTINUE
         ENDIF
10050 FORMAT(' NOEUD',I8,' X=',G15.7,' Y=',G15.7,' Z=',G15.7)
11050 FORMAT(' NODE', I8,' X=',G15.7,' Y=',G15.7,' Z=',G15.7)

C        TRACE DE L'EF DE QUALITE MINIMALE
C        ---------------------------------
         CALL XVEPAISSEUR( 3 )
         MN = MNXYZN + WYZNOE -1
         IF( NDIM .EQ. 2 ) THEN

C           TRACE DU TRIANGLE OU QUADRANGLE
            DO J=1,NBSOEF
               MN1 = MN + 3*NONOEF(J) - 3
               DO I=1,3
                  XYZARE(I,J) = RMCN( MN1 + I )
               ENDDO
            ENDDO
            CALL TRAITS3D( NCROUG, NBSOEF, XYZARE )

         ELSE IF( NDIM .EQ. 3 ) THEN

C           TRACE DU TETRAEDRE ou PYRAMIDE ou PENTAEDRE ou HEXAEDRE
C           NOMBRE DE FACES DE L'EF
            NFACE = NBFACE( NCOGEL )

C           NOMBRE DE SOMMETS DES FACES DE L'EF ET NO DES SOMMETS
            CALL SOFACU( NCOGEL, NBSOFA, NOSOFA )

            DO K=1,NFACE
               DO J=1,NBSOFA(K)
                  MN1 = MN + 3*NONOEF(NOSOFA(J,K)) - 3
                  DO I=1,3
                     XYZARE(I,J) = RMCN( MN1 + I )
                  ENDDO
               ENDDO
               CALL TRAITS3D( NCROUG, NBSOFA(K), XYZARE )
            ENDDO
         ENDIF

C        TRACE DU NO DE NOEUD DES SOMMETS ET DU NO DE L'EF
         DO I=1,3
            XYZBAR(I) = 0
         ENDDO
         DO J=1,NBSOEF

C           TRACE DU NUMERO DU NOEUD SOMMET J DE L'EF
            MN1 = MN + 3 * NONOEF(J) - 3
            CALL ENTIER3D( NCMAGE, RMCN(MN1+1), NONOEF(J) )

C           XYZ DU BARYCENTRE DE L'EF
            DO I=1,3
               XYZBAR(I) = XYZBAR(I) + RMCN(MN1+I)
            ENDDO
         ENDDO
C        XYZ DU BARYCENTRE DE L'EF POUR TRACE LE NO DE l'EF
         DO I=1,3
            XYZBAR(I) = XYZBAR(I) / NBSOEF
         ENDDO

C        TRACE DU NUMERO NEFMIN DE L'EF DE QUALITE MINIMALE
         CALL ENTIER3D( NCNOIR, XYZBAR, NEFMIN )
         CALL XVEPAISSEUR( 0 )

      ENDIF

C     AFFICHAGE DU MIN ET MAX DE CHAQUE COORDONNEE XYZ...
C     ---------------------------------------------------
C     LA BOUCLE SUR LES NBCOOR COORDONNEES DES NBNOEU NOEUDS
C     CALCUL DU MIN ET MAX DE CHAQUE COORDONNEE
      MN = MNXYZN + WYZNOE -1
      DO I=1,NBCOOR
         XYZEFS(I,1) = RMCN( MN + I )
         XYZEFS(I,2) = RMCN( MN + I )
      ENDDO
      DO J=2,NBNOEU
         DO I=1,NBCOOR
C           LA COORDONNEE I DU POINT J
            R = RMCN( MN + I )
            IF( R .LT. XYZEFS(I,1) ) THEN
               XYZEFS(I,1) = R
            ENDIF
            IF( R .GT. XYZEFS(I,2) ) THEN
               XYZEFS(I,2) = R
            ENDIF
         ENDDO
         MN = MN + NBCOOR
      ENDDO

      IF( LORBITE .EQ. 1 ) THEN
         WRITE(IMPRIM,*)
         DO I=1,NBCOOR
            WRITE(IMPRIM,10060) I, XYZEFS(I,1), XYZEFS(I,2)
         ENDDO
      ENDIF
10060 FORMAT(' OBJT COORD',I1,' MIN=',G15.7,' MAX=',G15.7)

 9999 IF( MNNBST .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNBST, MNNBST )

C     AFFICHAGE D'UN AVERTISSEMENT SI LA QUALITE MINIMALE QUAMIN
C     EST TROP PETITE
C     ----------------------------------------------------------
      IF( LORBITE .EQ. 1 .AND. QUAMIN .LE. 0.01 ) THEN
         CALL MEMPXFENETRE
         NBLGRC(NRERR) = 1
         WRITE(KERR(5)(1:15),'(G15.6)') QUAMIN
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='impquaob: OBJET ' // KNOMOB
            KERR(2) ='MAILLAGE avec une TRES MAUVAISE QUALITE MINIMALE='
     %                 // KERR(5)(1:15)
         ELSE
            KERR(1) = 'impquaob: OBJECT ' // KNOMOB
            KERR(2) = 'MESH has a VERY BAD MINIMUM QUALITY='
     %                 // KERR(5)(1:15)
         ENDIF
         CALL SANSDBL( KERR(2), LL )
         CALL SANSDBL( KERR(1), L  )
         KERR(1) = KERR(1)(1:L) // ' : ' // KERR(2)(1:LL)
         CALL LAVERT
      ENDIF

C     SEPARATEUR DE FIN D'OBJET CREE
      IF( LORBITE .EQ. 1 ) WRITE(IMPRIM,19000)
19000 FORMAT(100('-'))

      IF( LORBITE .GE. 1 ) LORBITE = LORBITE + 1

C     MISE A JOUR DU CADRE MAXIMAL DES MAILLAGES DES OBJETS
C     -----------------------------------------------------
      CALL MAJXYZEXT( NBCOOR, XYZEFS )

      RETURN
      END
