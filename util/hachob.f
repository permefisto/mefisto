      SUBROUTINE HACHOB( NMOBJT, NBMOTS, NTFAOB, MNFAOB, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREER SI CE N'EST DEJA FAIT LE TABLEAU DES FACES D'UN OBJET
C -----    A PARTIR DE SA TOPOLOGIE ET CHAINER LES FACES FRONTALIERES
C          cf ~/td/d/a___face
C          Attention: ICI PAS DE FACES INTERFACES RECENSEES ENTRE 2 VOLUMES
C
C ENTREES:
C --------
C NMOBJT : NOM DE L'OBJET
C NBMOTS : NOMBRE DE MOTS D'INFORMATIONS A STOCKER EN PLUS PAR FACE(>=4)
C
C SORTIES:
C --------
C NTFAOB : NUMERO DU TMS DES FACES DE L'OBJET
C          0 SI NON CREE PAR CAUSE D'ERREUR
C MNFAOB : ADRESSE MCN DU TABLEAU DES FACES DE L'OBJET
C IERR   : =0 SI PAS D'ERREUR
C          =1 SI PAS D'EF VOLUMIQUES
C          >0 SI ERREUR
C
C EN SORTIE LE TABLEAU LFACES CONTIENT :
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE OU 0 SI TRIANGLE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C          SI SOMMET 2 < DERNIER SOMMET  => FACE   DIRECTE DANS LE CUBE
C                      >                 => FACE INDIRECTE
C          UNE FACE DIRECTE EST VUE DE L EXTERIEUR AU CUBE
C          SOUS LA FORME DIRECTE
C          (CUBE SIGNIFIE ICI UN EF 3D: tetra, penta, hexaedre)
C          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CE CUBE
C                       <0 SI FACE INDIRECTE DANS CE CUBE
C          SI LA FACE APPARTIENT A 2 CUBES ALORS
C          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CE CUBE
C                       <0 SI FACE INDIRECTE DANS CE CUBE
C          SINON
C          LFACES(7,I)= NUMERO DE LA FACE FRONTALIERE SUIVANTE
C                       0 SI C'EST LA DERNIERE
C          L1FAFR = NUMERO DE LA PREMIERE FACE FRONTALIERE DANS LFACES
C          L1FA2M = NUMERO DE LA PREMIERE FACE INTERFACE   DANS LFACES
C
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
C
C        ( LFACES(9,I)= NUMERO NUTYEL DU TYPE DU DERNIER EF DE CETTE FACE
C                       CETTE INFORMATION EXISTE SEULEMENT POUR UN OBJET
C                       PAS POUR UNE SURFACE OU UN VOLUME!
C          DE PLUS CE NO PEUT ETRE REMPLACE PAR NUTYEL*1 000 000 + NOSURF
C          OU NOSURF EST LE NUMERO DE LA SURFACE DE CETTE FACE FRONTALIERE
C          DANS LE SP openfoam.f )
C
C          NUTGFA(1:8,1:NBFATG)= NUMERO DES 8 TANGENTES DES FACES A TG
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C ...................................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___face.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/ponoel.inc"
      include"./incl/ntmnlt.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMOBJT
      INTEGER           NTELEM(9), MNELEM(9)
      INTEGER           NGS(4), NGS1(4)
      LOGICAL           AVANT
      REAL              X(4),Y(4),Z(4)
C
      IERR = 0
C
C     L'OBJET EST RECHERCHE
C     =====================
C     LE TABLEAU LEXIQUE DE CET OBJET
      CALL LXLXOU( NTOBJE, NMOBJT, NTLXOB, MN )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET INCONNU: ' // NMOBJT
         ELSE
            KERR(1) = 'UNKNOWN OBJECT: ' // NMOBJT
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE DES ADRESSES DES TABLEAUX DE LA TOPOLOGIE DE L'OBJET
C     ==============================================================
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTPOGE : NUMERO      DU TMS 'XYZPOINT'  DE L'OBJET
C     MNPOGE : ADRESSE MCN DU TMS 'XYZPOINT'  DE L'OBJET
C     NTNOEU : NUMERO      DU TMS 'XYZNOEUD'  DE L'OBJET
C     MNNOEU : ADRESSE MCN DU TMS 'XYZNOEUD'  DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE
C     NTELEM : NUMERO      DU TMS 'NPEF' DES NBTYEL TYPES D'EF
C     MNELEM : ADRESSE MCN DU TMS 'NPEF' DES NBTYEL TYPES D'EF
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO ,
     %             NTPOGE, MNPOGE, NTNOEU, MNNOEU ,
     %             NBTYEL, NTELEM, MNELEM, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE TOTAL DE POINTS DE L'OBJET
      NBPOI = MCN( MNPOGE + WNBPOI )
C     LE NOMBRE DE COORDONNEES DES POINTS
      NBCOOR = MCN( MNPOGE + WBCOOR )
C
C     LE TABLEAU 'FACE' DE CET OBJET
      CALL LXTSOU( NTLXOB, 'FACE', NTFAOB, MNFAOB )
      IF( NTFAOB .GT. 0 ) THEN
C        LES FACES SONT ELLES ANTERIEURES A LA TOPOLOGIE DE L'OBJET?
         IF( AVANT( MCN(MNTOPO), MCN(MNFAOB) ) ) THEN
C           LA TOPOLOGIE EST ANTERIEURE AUX FACES
C           LE TABLEAU FACE EST ACCEPTE  => RETOUR SANS RECALCUL
            RETURN
         ENDIF
C        LE TABLEAU FACE EST DETRUIT POUR ETRE ENSUITE REGENERE
         CALL LXTSDS( NTLXOB, 'FACE' )
         NTFAOB = 0
         MNFAOB = 0
      ENDIF
C
C     RECHERCHE DU NOMBRE DES EF VOLUMIQUES DE L'OBJET
      NBTETR = 0
      NBPENT = 0
      NBPYRA = 0
      NBHEXA = 0
      NBTYEL = MCN(MNTOPO+WBTYEL)
      IF( NBCOOR .EQ. 6 ) THEN
C         TRACE RESTREINT AU 3Q1C DES 6-CUBES
          N1TYEL=2
          N2TYEL=2
       ELSE
C         TOUS LES TYPES D'EF
          N1TYEL=1
          N2TYEL=NBTYEL
       ENDIF
      DO 20 I=N1TYEL,N2TYEL
         NUTYEL = MCN(MNELEM(I)+WUTYEL)
         CALL ELNUCG( NUTYEL, NCOGEL )
         GOTO ( 20, 20, 20, 20, 15, 16, 17, 17, 19 ), NCOGEL
C        TETRAEDRES
 15      NBTETR = NBTETR + MCN( MNELEM(I) + WBELEM )
         GOTO 20
C        PENTAEDRES
 16      NBPENT = NBPENT + MCN( MNELEM(I) + WBELEM )
         GOTO 20
C        HEXAEDRES OU 6-CUBES
 17      NBHEXA = NBHEXA + MCN( MNELEM(I) + WBELEM )
         GOTO 20
C        PYRAMIDES
 19      NBPYRA = NBPYRA + MCN( MNELEM(I) + WBELEM )
 20   CONTINUE
      N = NBTETR + NBPENT + NBHEXA + NBPYRA
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'NOMBRE TETRAEDRES =',NBTETR
         WRITE(IMPRIM,*) 'NOMBRE PENTAEDRES =',NBPENT
         WRITE(IMPRIM,*) 'NOMBRE HEXAEDRES  =',NBHEXA
         WRITE(IMPRIM,*) 'NOMBRE PYRAMIDES  =',NBPYRA
         WRITE(IMPRIM,*) 'NOMBRE TOTAL EF 3D=',N
      ELSE
         WRITE(IMPRIM,*) 'NUMBER of TETRAHEDRA =',NBTETR
         WRITE(IMPRIM,*) 'NUMBER of PENTAHEDRA =',NBPENT
         WRITE(IMPRIM,*) 'NUMBER of HEXAHEDRA  =',NBHEXA
         WRITE(IMPRIM,*) 'NUMBER of PYRAMIDS   =',NBPYRA
         WRITE(IMPRIM,*) 'TOTAL NUMBER of 3D FE=',N
      ENDIF
C
      IF( N .LE. 0 ) THEN
C        PAS DE VOLUME => PAS DE FACES FRONTALIERES
         NBLGRC(NRERR) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET ' //  NMOBJT
            KERR(2) = 'SANS ELEMENT FINI VOLUMIQUE'
            KERR(3) = 'REVOYEZ VOS COORDONNEES'
            KERR(4) = 'ATTENTION: EN AXISYMETRIE X est R et Y est Z!'
         ELSE
            KERR(1) = 'OBJECT ' //  NMOBJT
            KERR(2) = 'WITHOUT VOLUMIC FINITE ELEMENT'
            KERR(3) = 'REVIEW THE COORDINATES of VERTICES'
            KERR(4) = 'ATTENTION: IN AXISYMMETRY X is R and Y is Z!'
         ENDIF
         CALL LEREUR
         NTFAOB = 0
         MNFAOB = 0
         IERR   = 1
         RETURN
      ENDIF
C
C     NOMBRE D'ENTIERS POUR UNE FACE
      MOFACE = 5 + NBMOTS
C
C     MAJORATION DU NOMBRE DE FACES : FORMULE A MODIFIER EVENTUELLEMENT
      IF( N .LT. 1000 .OR. NBCOOR .EQ. 6 ) THEN
         MXFACE = 4 * NBTETR + 5 * NBPENT + 6 * NBHEXA + 5 * NBPYRA
      ELSE
C        0.6 CORRECT SI PAS D'ERREUR POUR UN EF SUBDIVISE EN N**3 SOUS-EF
         MXFACE=INT( (4*NBTETR + 5*NBPENT + 6*NBHEXA + 5*NBPYRA) * 0.61)
      ENDIF
C
C     ADRESSAGE ET OUVERTURE DU TABLEAU FACE SANS STOCKAGE DES XYZ DES SOMMETS
C     ------------------------------------------------------------------------
 30   MOTFAC = WFACES + MOFACE*MXFACE
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'hachob DEMANDE',MOTFAC,
     %          ' MOTS du TMS ''FACE'' de l''OBJET ',NMOBJT
      ELSE
        PRINT*,'hachob REQUESTS',MOTFAC,
     %         ' WORDS of TMS ''FACE'' of the OBJECT ',NMOBJT
      ENDIF
      CALL TNMCMX( 'ENTIER', MAXVAR )
      IF( MAXVAR .LT. MOTFAC ) THEN
         NTFAOB = 0
         MNFAOB = 0
         GOTO 33
      ENDIF
 
      CALL LXTNDC( NTLXOB, 'FACE', 'ENTIER', MOTFAC )
      CALL LXTSOU( NTLXOB, 'FACE', NTFAOB, MNFAOB )
 33   IF( NTFAOB .LE.0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'hachob: OBJET ' // NMOBJT
            KERR(2) = 'PAS ASSEZ DE MEMOIRE POUR DECLARER le TMS FACE'
            KERR(3) = '                MOTS DIPONIBLES'
         ELSE
            KERR(1) = 'hachaob: OBJECT ' // NMOBJT
            KERR(2) = 'NOT ENOUGH MEMORY TO STORE TMS FACE'
            KERR(3) = '                WORDS ALLOWED'
         ENDIF
         WRITE(KERR(3)(1:15),'(I15)') MAXVAR
         CALL LEREUR
         IERR   = 2
         NTFAOB = 0
         MNFAOB = 0
         RETURN
      ENDIF
C
C     LE TABLEAU DES FACES EST INITIALISE A ZERO
      MNFACE = MNFAOB + WFACES
      CALL AZEROI( MOFACE*MXFACE, MCN(MNFACE) )
C
C     LA 1-ERE FACE LIBRE
      LIBREF = MXFACE
      NBFAPB = 0
C
      DO 2000 NTYEF=N1TYEL,N2TYEL
         NUTYEL = MCN(MNELEM(NTYEF)+WUTYEL)
         CALL ELNUCG( NUTYEL, NCOGEL )
         GOTO ( 2000, 2000, 2000, 2000, 40, 40, 40, 40, 40 ), NCOGEL
C
C        LE NOMBRE D'EF DE CE TYPE
 40      NBELEM = MCN( MNELEM(NTYEF) + WBELEM )
C        LE NOMBRE DE NOEUDS DE L'EF
         NBNDEL = MCN( MNELEM(NTYEF) + WBNDEL )
C        LE NOMBRE DE POINTS GEOMETRIQUES D'UN EF DE CE TYPE
         NBPGEL = MCN( MNELEM(NTYEF) + WBPGEL )
         IF( NBPGEL .EQ. 0 ) THEN
C           POINTS=NOEUDS   NUNDEL JOUE LE ROLE DE NUPGEL
C           LE DEBUT DU TABLEAU DES NUMEROS DES POINTS DU MAILLAGE
            MNEF   = MNELEM(NTYEF) + WUNDEL - 1
            NBPGEL = NBNDEL
         ELSE
C          POINTS/=NOEUDS
C          LE DEBUT DU TABLEAU DES NUMEROS DES POINTS DU MAILLAGE
           MNEF = MNELEM(NTYEF) + WUNDEL + NBELEM*NBNDEL - 1
         ENDIF
C
C        LE NOMBRE DE FACES DE CE TYPE D'EF
         NFACE = NBFACE( NCOGEL )
C
C        LE NUMERO DES SOMMETS DES FACES DE L'EF
         CALL SOFACU( NCOGEL, NBSOFA, NOSOFA )
C
         DO 1000 N=1,NBELEM
C
C           BOUCLE SUR LES FACES DU CUBE
C           ----------------------------
            DO 200 NF=1,NFACE
C
C              LE NOMBRE DE SOMMETS DE LA FACE
               NBSTF = NBSOFA( NF )
C
C              LE NUMERO NOMIN MIN DES NUMEROS DES SOMMETS DE LA FACE
               ISENS = 0
               NOMIN = NBPOI
               DO J=1,NBSTF
c                 LE NUMERO DU SOMMET J DE LA FACE NF DU CUBE N
                  K      = MCN( MNEF + N + (NOSOFA(J,NF)-1) * NBELEM )
                  NGS(J) = K
                  IF( K .LT. NOMIN ) THEN
                     NOMIN = K
                     ISENS = J
                  ENDIF
               ENDDO
C
C              PERMUTATION CIRCULAIRE DES SOMMETS POUR AMENER NGS(J)
C              LE PLUS PETIT NO EN PREMIERE POSITION DE NGS1
               DO J=1,NBSTF
                  NGS1(J) = NGS(ISENS)
                  ISENS   = ISENS + 1
                  IF( ISENS .GT. NBSTF ) ISENS = 1
               ENDDO
C
C              LA NUMEROTATION DES SOMMETS EST MISE A JOUR POUR
C              QUE LA 1-ERE ARETE DE LA FACE SOIT
C              1:NO LE PLUS FAIBLE DES SOMMETS DE LA FACE
C              2:NO MIN(SOMMET2, SOMMET NBSTF)
C              SI MIN = SOMMET 2  FACE DIRECTE
C                       SOMMET NBSTF FACE INDIRECTE
C              UNE FACE DIRECTE EST VUE DE L EXTERIEUR AU SOLIDE
C              SOUS LA FORME DIRECTE
               J = NGS1( 2 )
               IF( J .GT. NGS1(NBSTF) ) THEN
C                 FACE INDIRECTE
                  ISENS   = -1
                  NGS1(2) = NGS1(NBSTF)
                  NGS1(NBSTF) = J
               ELSE
                  ISENS = 1
               ENDIF
C
C              RECHERCHE OU ADJONCTION DE LA FACE
C              ----------------------------------
               CALL HACHAG( NBSTF, NGS1, MOFACE, MXFACE, MCN(MNFACE), 5,
     &                      LIBREF,   NOFAC )
               IF( NOFAC .EQ. 0 ) THEN
C                 LISTE SATUREE CAR MXFACE SOUS-ESTIMEE
                  CALL LXTSDS( NTLXOB, 'FACE' )
C                 MAJORATION SAUF SI DES EF DISJOINTS
                  MXFACE = 4*NBTETR + 5*NBPENT + 6*NBHEXA + 5*NBPYRA
                  GOTO 30
               ENDIF
               MN = MNFACE + MOFACE * ( ABS(NOFAC) - 1 )
C
C              STOCKAGE ADRESSE DE L'EF DANS NFACE(6,...;NOFAC)
C              ------------------------------------------------
               K = 4
 140           K = K + 1
               IF( K .LT. 7 ) THEN
C
C                 RECHERCHE D'UN NUMERO NUL DE EF
                  IF( MCN(MN + K) .NE. 0) GOTO 140
C
C                 NO EF > 0 SI LA FACE EST   DIRECTE DANS L'EF
C                       < 0 SI LA FACE EST INDIRECTE DANS L'EF
                  MCN( MN + K ) = N * ISENS
C
C                 LE NO (1 A NBTYEL) DU TYPE DE L'EF DANS LE MAILLAGE
                  MCN( MN + 8 ) = NTYEF
C
C                 ajouter ici le calcul des numeros des tangentes
C                 des faces a tg
C
               ELSE
C
C                 ERREUR: IL Y A PLUS DE 2 CUBES CONTENANT CETTE FACE
                  IF( NBFAPB .EQ. 0 ) THEN
C
C                    INITIALISATION DU TRACE DES FACES A PROBLEME
C                    AUCUN ITEM SUR L'ECRAN
                     CALL EFFACE
                     CALL ITEMS0
C                    LES PARAMETRES DU CADRE MAXIMAL
                     CALL VISEE0
C
                  ENDIF
C
C                 UNE FACE A PB DE PLUS
                  NBFAPB = NBFAPB + 1
C
C                 IMPRESSION ET AFFICHAGE DE CETTE FACE A PB
                  MM = MNPOGE + WYZPOI - 4
                  DO 190 L=1,NBSTF
C                    LE NUMERO DU SOMMET ET SES XYZ SONT STOCKES DANS X Y Z
                     NOSM = NGS1(L) * 3
                     X(L) = RMCN(MM+NOSM+1)
                     Y(L) = RMCN(MM+NOSM+2)
                     Z(L) = RMCN(MM+NOSM+3)
 190              CONTINUE
C
C                 TRACE DE LA FACE: UN TRIANGLE OU QUADRANGLE
                  IF( NDIMLI .EQ. 3 ) THEN
                     CALL FAP13D( NCROUG, NCNOIR, 0, NBSTF, X, Y, Z )
                  ELSE
                     CALL FACE2D( NCROUG, NCNOIR, NBSTF, X, Y )
                  ENDIF
C                 IMPRESSION DE CETTE FACE A PB
                  IF( IERR .EQ. 0 ) THEN
                     IF( NBCOOR .NE. 6 ) THEN
                        WRITE(IMPRIM,10200) N
                        NTT1 = ABS(MCN(MN+5))
                        NTT2 = ABS(MCN(MN+6))
                        NBLGRC(NRERR) = 2
                        IF( LANGAG .EQ. 0 ) THEN
                           WRITE(IMPRIM,10201) NF,(NGS1(L),L=1,NBSTF)
                           WRITE(IMPRIM,10202)
     %                          (NGS1(L),X(L),Y(L),Z(L),L=1,NBSTF)
                WRITE(IMPRIM,10203) NTT1,NTT2,N
                WRITE(IMPRIM,10204) NTT1,(MCN(MNEF+NTT1+m*NBELEM),m=0,3)
                WRITE(IMPRIM,10204) NTT2,(MCN(MNEF+NTT2+m*NBELEM),m=0,3)
                WRITE(IMPRIM,10204) N   ,(MCN(MNEF+N   +m*NBELEM),m=0,3)
                      KERR(1) = 'PLUS DE 2 CUBES CONTIENNENT CETTE FACE'
                      KERR(2) =
     %    'LES AUTRES FACES DANS CETTE SITUATION NE SONT PLUS INDIQUEES'
                      ELSE
                           WRITE(IMPRIM,20201) NF,(NGS1(L),L=1,NBSTF)
                           WRITE(IMPRIM,20202)
     %                          (NGS1(L),X(L),Y(L),Z(L),L=1,NBSTF)
                WRITE(IMPRIM,20203) NTT1,NTT2,N
                WRITE(IMPRIM,20204) NTT1,(MCN(MNEF+NTT1+m*NBELEM),m=0,3)
                WRITE(IMPRIM,20204) NTT2,(MCN(MNEF+NTT2+m*NBELEM),m=0,3)
                WRITE(IMPRIM,20204) N   ,(MCN(MNEF+N   +m*NBELEM),m=0,3)
C
                         KERR(1) = 'MORE THAN 2 CUBES CONTAIN this FACE'
                         KERR(2) =
     %                 'After, the other SIMILAR FACES are not reported'
                      ENDIF
                      CALL LEREUR
                      IERR = 6
                     ENDIF
                  ENDIF
C
               ENDIF
 200        CONTINUE
10200    FORMAT(/' HACHOB: CUBE:',I6)
10201    FORMAT(' LA FACE',I3,' DE SOMMETS:',T30,4I6)
10202    FORMAT(' SOMMET',I6,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
10203    FORMAT(' APPARTIENT A PLUS DES 2 CUBES', 3I10)
10204    FORMAT(' EF',I9,' SOMMETS',4I9)
C
20201    FORMAT(' The FACE',I3,' of VERTICES:',T30,4I6)
20202    FORMAT(' VERTEX',I6,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
20203    FORMAT(' BELONGS TO MORE 2 CUBES', 3I10)
20204    FORMAT(' FE',I9,' VERTICES',4I9)
C
 1000    CONTINUE
 2000 CONTINUE
      IF( NBCOOR .NE. 6 .AND. NBFAPB .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) NBFAPB,
     %  ' FACES APPARTIENNENT A PLUS DE 2 CUBES => MAILLAGE avec ERREUR'
         ELSE
            WRITE(IMPRIM,*) NBFAPB,
     %    ' FACES BELONG TO MORE 2 CUBES => MESH with ERROR'
         ENDIF
         RETURN
      ENDIF

C     CHAINAGE DES FACES FRONTALIERES C-A-D APPARTENANT A UN SEUL EF
C     --------------------------------------------------------------
      NBFACES= 0
      MN     = MNFACE
      L1FAFR = 0
      NBFAFR = 0
      DO 3000 N=1,MXFACE

C        SAUT DES FACES NON UTILISEES DANS MNFACE
         IF( MCN( MN ) .EQ. 0 ) GOTO 2900

         NBFACES = NBFACES + 1
C        SAUT DES FACES APPARTENANT A 2 CUBES
         IF( MCN( MN + 6 ) .NE. 0 ) GOTO 2900

C        LA FACE EST FRONTALIERE. ELLE EST CHAINEE AVEC LA PRECEDENTE
         MCN( MN + 6 ) = L1FAFR
         L1FAFR        = N
         NBFAFR        = NBFAFR + 1

C        PASSAGE A LA FACE SUIVANTE
 2900    MN = MN + MOFACE

 3000 CONTINUE

C     LE TABLEAU 'FACE' EST COMPLETE
C     ==============================
C     LE NOMBRE D'ENTIERS PAR FACE
      MCN( MNFAOB + WOFACE ) = MOFACE
C     LA MAJORATION DU NOMBRE DE FACES
      MCN( MNFAOB + WXFACE ) = MXFACE
C     LE NUMERO DE LA PREMIERE FACE FRONTALIERE
      MCN( MNFAOB + W1FAFR ) = L1FAFR
C     LE NOMBRE DE FACES FRONTALIERES
      MCN( MNFAOB + WBFAFR ) = NBFAFR
C
C     ATTENTION: ICI PAS DE FACES INTERFACE ENTRE 2 MATERIAUX RECENSEES!
C     LE NUMERO DE LA PREMIERE FACE INTERFACE (ADJACENTE A 2 MATERIAUX)
      MCN( MNFAOB + W1FA2M ) = 0
C     LE NOMBRE DE FACES INTERFACE (ADJACENTE A 2 MATERIAUX)
      MCN( MNFAOB + WBFA2M ) = 0

C     LE NOMBRE DE FACES AVEC TANGENTES
      MCN( MNFAOB + WBFATG ) = 0
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNFAOB) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFAOB + MOTVAR(6) ) = NONMTD( '~>>>FACE' )

      IF(LANGAG .EQ. 0 ) THEN
         PRINT*,'hachob CORRECT',MXFACE,' FACES DECLAREES pour',
     %           NBFACES,'  FACES EFFECTIVES de l''OBJET ',NMOBJT
      ELSE
         PRINT*,'hachob OK',MXFACE,' DECLARED FACES for',
     %           NBFACES,' USED FACES of the OBJET ',NMOBJT
      ENDIF

      RETURN
      END
