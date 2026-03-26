      SUBROUTINE DFTOP4( NTLXOB, MOELEM, NBELEM, MNELEM,
     %                   MOELTG, NBELTG, MNELTG, NUELTG, LIENEL,
     %                   NOINTE, NBOBPR, NUOBPR, XYZ,
     %                   NTNPEF, MNNPEF, NMNPEF, NUOBIN, NUOBCL,
     %                   IERR   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    INITIALISER LES TABLEAUX NPEF"TYPE_EF DE L'OBJET
C -----    LES NOEUDS ET POINTS NON SOMMETS SERONT AJOUTES DANS DFTOP5
C
C ENTREES :
C ---------
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET
C MOELEM : NOMBRE DE MOTS DE CHAQUE TYPE D'ELEMENT
C NBELEM : NBELEM(1)=NOMBRE DE NOEUDSOMMETS
C          NBELEM(2)=SEGMENTS
C          NBELEM(3)=TRIANGLES
C          NBELEM(4)=QUADRANGLES
C          NBELEM(5)=TETRAEDRES
C          NBELEM(6)=PENTAEDRES
C          NBELEM(7)=HEXAEDRES
C          NBELEM(8)=6-CUBES
C MNELEM : ADRESSE MCN DES EVENTUELS TABLEAUX DES 9 TYPES D'EF
C          0 SI PAS DE D'ELEMENTS FINIS DE CE TYPE
C          MCN(MNELEM(NCOGEL))=LISTE DE HACHAGE DES EF DE CODE NCOGEL
C          CETTE LISTE A ETE CALCULEE PAR UN CALL DFTOP3
C MOELTG : NOMBRE DE TG PAR TYPE D'EF
C NBELTG : NOMBRE D'EF A TG PAR TYPE
C MNELTG : ADRESSE MCN DES EVENTUELS 8 TABLEAUX DES ELEMENTS FINIS A TG
C NUELTG ; NOMBRE D'EF A TG DECLARES DANS DFTOP3
C LIENEL : POSITION DU CHAINAGE POUR LE SP HACHAG
C NOINTE : NUMERO DE L'INTERPOLATION CHOISIE  (CF DFTOPO)
C          1  'AXISYMETRIQUE_DEGRE_1',
C          2  'AXISYMETRIQUE_DEGRE_2',
C          3  'LAGRANGE_DEGRE_1',
C          4  'LAGRANGE_DEGRE_2'
C NBOBPR : NOMBRE D'OBJETS PREMIERS
C NUOBPR : LE TYPEOBJET DE CHAQUE OBJET PREMIER
C XYZ    : 3 COORDONNEES DES SOMMETS DE L'OBJET
C
C SORTIES :
C ---------
C NTNPEF : NUMERO      DES TMS NPEF" PAR TYPE D'ELEMENTS FINIS
C MNNPEF : ADRESSE MCN DES TMS NPEF" PAR TYPE D'ELEMENTS FINIS
C NMNPEF : SUFFIXE DES TABLEAUX NPEF"
C NUOBIN : NUMEROS DES OBJETS INTERNES
C NUOBCL : NUMEROS DES OBJETS AUX LIMITES
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1996
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*1       KSUFIX
      PARAMETER        (KSUFIX='"')
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a___npef.inc"
      include"./incl/ponoel.inc"
      include"./incl/sotgfc.inc"
      CHARACTER*13      KNOM, NOMTMS(9)
      CHARACTER*4       NOM4
      CHARACTER*24      NMPLSV
      INTEGER           MOELEM(9),NBELEM(9),MNELEM(9),
     %                  MOELTG(9),NBELTG(9),MNELTG(9),NUELTG(9),
     %                  LIENEL(9),
     %                  NUOBPR(2,NBOBPR),
     %                  NTNPEF(9),MNNPEF(9),NMNPEF(9),NUOBIN(NBOBPR),
     %                  NUOBCL(NBOBPR)
      INTEGER           NOSO(64),NOST(64),NBEFPU(9),NTG(2),NUTG(2)
      REAL              XYZ(3,*)
      include"./incl/nomele.inc"
C
      IERR = 0
C
C     PROTECTION DU NOMBRE DES ELEMENTS FINIS
      CALL TRTATA( NBELEM, NBEFPU, 9 )
C
C     RECHERCHE DU CODE GEOMETRIQUE MINIMUM PAR DIMENSION
      IF( NBELEM(9) .GT. 0 .OR. NBELEM(7) .GT. 0 .OR.
     %    NBELEM(6) .GT. 0 .OR. NBELEM(5) .GT. 0 ) THEN
         LDIMEF = 3
      ELSE IF( NBELEM(4) .GT. 0 .OR. NBELEM(3) .GT. 0 ) THEN
         LDIMEF = 2
      ELSE IF( NBELEM(2) .GT. 0 ) THEN
         LDIMEF = 1
      ELSE
         LDIMEF = 0
      ENDIF
C
      IF( LDIMEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AUCUN VOLUME, AUCUNE SURFACE, AUCUNE LIGNE'
            KERR(2) = 'OBJET A REDEFINIR'
         ELSE
            KERR(1) = 'NO VOLUME, NO SURFACE, NO LINE'
            KERR(2) = 'DEFINE AGAIN the OBJECT'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
      DO 5 NCOGEL=9,1,-1
         NOMTMS( NCOGEL ) = '    '
 5    CONTINUE
C
C     NECESSAIRE POUR RETROUVER LE NUMERO DES TGS D'UNE FACE DES EF VOLUMIQUES
C     PAS DE TG POUR LES 6-CUBES
      DO 10 NCOGEL=5,9
C        LE NUMERO DES TANGENTES DES FACES DU TETRAEDRE, PENTAEDRE, HEXAEDRE
         CALL TGFACU( NCOGEL, NBTGFA(1,NCOGEL), NOTGFA(1,1,NCOGEL) )
 10   CONTINUE
C
C     ======================================================================
C     BOUCLE SUR LES PYRAMIDES, HEXAEDRES, PENTAEDRES, TETRAEDRES,
C                    QUADRANGLES, TRIANGLES, SEGMENTS, VERTICES
C     ======================================================================
      DO 1000 NCOGEL=9,1,-1
C
C        LE NUMERO DU TABLEAU TMS ~>OBJET>>NPEF"NCOGEL
         NTNPEF(NCOGEL) = 0
C
C        NON TRAITEMENT DES 6-CUBES
         IF( NCOGEL .EQ. 8 ) GOTO 1000
C
C        NBELT : LE NOMBRE ACTUEL D'ELEMENTS FINIS AVEC CE CODE GEOMETRIQUE
C        ATTENTION:  SI UNE FACE D'UN EF VOLUMIQUE EST IDENTIFIEE A UNE FACE
C                    D'UNE SURFACE, ALORS CETTE FACE N'EST PAS UN EF ET
C                    EST SUPPRIMEE DES EF SUSCEPTIBLES DE CREER UN TMS NPEF"
C                    DE MEME POUR UNE ARETE D'UN EF SURFACIQUE OU VOLUMIQUE
C                    DE MEME POUR UN  POINT D'UN EF LINEIQUE OU SURFACIQUE OU VO
         NBELT = NBEFPU(NCOGEL)
         IF( NBELT .LE. 0 ) GOTO 1000
C
         IF( LDIMEF .EQ. 3 .AND. NCOGEL .LE. 4 ) THEN
C           INCOMPATIBILITE DES FACES DES SURFACES ET DES EF VOLUMIQUES
C           AFFICHAGE DES XYZ DES SOMMETS DES EF NON IDENTIFIES
            CALL PERTEF( 4, NBELEM, MOELEM, MNELEM, XYZ, NBEFPE )
            IF( NBEFPE .LE. 0 ) GOTO 1000
C           AFFICHAGE DES EF PERDUS
            NBLGRC(NRERR) = 3
            WRITE(KERR(8)(1:10), '(I10)') NBELT
            IF( LANGAG .EQ. 0 ) THEN
               IF( NCOGEL .EQ. 4 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' QUADRANGLES'
                  KERR(2) = 'DES SURFACES DE L''OBJET'
               ELSE IF( NCOGEL .EQ. 3 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' TRIANGLES'
                  KERR(2) = 'DES SURFACES DE L''OBJET'
               ELSE IF( NCOGEL .EQ. 2 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' ARETES'
                  KERR(2) = 'DES LIGNES DE L''OBJET'
               ELSE
                  KERR(6) = KERR(8)(1:10) // ' SOMMETS'
                  KERR(2) = 'DES POINTS DE L''OBJET'
               ENDIF
               L = NUDCNB( KERR(6) )
               KERR(1) = 'ERREUR DFTOP4: ' // KERR(6)(1:L)
               KERR(3) = 'NON RETROUVES DANS LES EF VOLUMIQUES'
            ELSE
               IF( NCOGEL .EQ. 4 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' QUADRANGLES'
                  KERR(2) = 'of SURFACES of the OBJECT'
               ELSE IF( NCOGEL .EQ. 3 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' TRIANGLES'
                  KERR(2) = 'of SURFACES of the OBJECT'
               ELSE IF( NCOGEL .EQ. 2 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' EDGES'
                  KERR(2) = 'of LINES of the OBJECT'
               ELSE
                  KERR(6) = KERR(8)(1:10) // ' VERTICES'
                  KERR(2) = 'of POINTS of the OBJECT'
               ENDIF
               L = NUDCNB( KERR(6) )
               KERR(1) = 'ERROR DFTOP4: ' // KERR(6)(1:L)
               KERR(3) = 'NOT RETRIEVED in FE of VOLUMES'
            ENDIF
            CALL LEREUR
            IERR = 13
            GOTO 1000
         ENDIF
C
         IF( LDIMEF .EQ. 2 .AND. NCOGEL .LE. 2 ) THEN
C           INCOMPATIBILITE DES  ARETES DES LIGNES ET DES EF SURFACIQUES
C           AFFICHAGE DES XYZ DES SOMMETS DES EF NON IDENTIFIES
            CALL PERTEF( 2, NBELEM, MOELEM, MNELEM, XYZ, NBEFPE )
            IF( NBEFPE .LE. 0 ) GOTO 1000
C           AFFICHAGE DES EF PERDUS
            NBLGRC(NRERR) = 3
            WRITE(KERR(8)(1:10), '(I10)') NBELT
            IF( LANGAG .EQ. 0 ) THEN
               IF( NCOGEL .EQ. 2 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' ARETES'
                  KERR(2) = 'DES LIGNES DE L''OBJET'
               ELSE
                  KERR(6) = KERR(8)(1:10) // ' SOMMETS'
                  KERR(2) = 'DES POINTS DE L''OBJET'
               ENDIF
               L = NUDCNB( KERR(6) )
               KERR(1) = 'ERREUR DFTOP4: ' // KERR(6)(1:L)
               KERR(3) = 'NON RETROUVES DANS LES EF SURFACIQUES'
            ELSE
               IF( NCOGEL .EQ. 2 ) THEN
                  KERR(6) = KERR(8)(1:10) // ' EDGES'
                  KERR(2) = 'of  LINES of the OBJECT'
               ELSE
                  KERR(6) = KERR(8)(1:10) // ' VERTICES'
                  KERR(2) = 'of POINTS of the OBJECT'
               ENDIF
               L = NUDCNB( KERR(6) )
               KERR(1) = 'ERROR DFTOP4: ' // KERR(6)(1:L)
               KERR(3) = 'NOT RETRIEVED in FE of SURFACES'
            ENDIF
            CALL LEREUR
            IERR = 12
            GOTO 1000
         ENDIF
C
         IF( LDIMEF .EQ. 1 .AND. NCOGEL .LE. 1 ) THEN
C           INCOMPATIBILITE DES SOMMETS ET DES EF LINEIQUES
C           AFFICHAGE DES XYZ DES SOMMETS DES EF NON IDENTIFIES
            CALL PERTEF( 1, NBELEM, MOELEM, MNELEM, XYZ, NBEFPE )
            IF( NBEFPE .LE. 0 ) GOTO 1000
C           AFFICHAGE DES EF PERDUS
            WRITE(KERR(8)(1:10), '(I10)') NBELT
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR DFTOP4: '// KERR(8)(1:10) // ' POINTS'
               KERR(2) = 'NON RETROUVES SOMMETS DES EF LINEIQUES'
            ELSE
               KERR(1) = 'ERROR DFTOP4: '// KERR(8)(1:10) // ' POINTS'
               KERR(2) = 'NOT RETRIEVED in FE of LINES'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 1000
         ENDIF
C
C        IL EXISTE DES ELEMENTS FINIS DE CODE GEOMETRIQUE NCOGEL
C        LE NUMERO DE TYPE DE L'ELEMENT FINI POUR L'INTERPOLATION CHOISIE
         NUTYEL = NOTYEL( NCOGEL, NOINTE )
         IF( NUTYEL .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'DFTOP4: TYPE D''ELEMENT FINI NON PROGRAMME'
            ELSE
               KERR(1) = 'DFTOP4: TYPE of FINITE ELEMENT NOT PROGRAMMED'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 1000
         ENDIF
C
C        TRAITEMENT PARTICULIER DU TRIANGLE P1 A CAUSE DES ESTIMATEURS
C        D'ERREUR AVEC 2 POINTS SUR L'ARETE A RENDRE COMPATIBLE AVEC
C        LE QUADRANGLE 2Q1C
         IF( NCOGEL .EQ. 3 .AND. NUTYEL .EQ. 13 .AND.
     %       NBELEM(4) .GT. 0 ) THEN
C           TRIANGLE 2P1D EN PRESENCE DE QUADRANGLES
C           => LE TRIANGLE DEVIENT TRIA 2P1C
C              AVEC 2 POINTS D'INTEGRATION SUR CHAQUE ARETE
            NUTYEL = 29
         ENDIF
C
C        LE SECOND MOT DU NOM DE L'ELEMENT FINI
         NOM4   =  NOMELE( 2, NUTYEL )
         KNOM   = 'NPEF' // KSUFIX // NOM4
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C        NBPOE  : NOMBRE DE POINTS DE L ELEMENT NOMELE
C        NBNOE  : NOMBRE DE NOEUDS DE L ELEMENT NOMELE
C        NOTRAE : CODE TRAITEMENT  DE L ELEMENT NOMELE
C        NBNSOM : NOMBRE DE NOEUDS-SOMMETS DE L ELEMENT
C        NBPOIN : NOMBRE DE POINTS INTERNES DE L ELEMENT
C        NBNOIN : NOMBRE DE NOEUDS INTERNES DE L ELEMENT
C        NARET  : NOMBRE DE SES ARETES
C        NOSOAR : NO DES 2 SOMMETS DE CHACUNE DE SES ARETES
C        NOTYAR : NO DU TYPE       DE CHACUNE DE SES ARETES
C                 ( CF SP TYARCP )
C        NBPOAR : NOMBRE DE POINTS-NON SOMMETS DE CHACUNE DE SES ARETES
C        NOPOAR : NO ELEMENTAIRE DES POINTS NON SOMMETS DE CHACUNE DE SES
C                 ARETES
C        NBNOAR : NOMBRE DE NOEUDS-NON SOMMETS DE CHACUNE DE SES ARETES
C        NONOAR : NO ELEMENTAIRE DES NOEUDS NON SOMMETS DE CHACUNE DE SES
C                 ARETES
C        NFACE  : NOMBRE DE SES FACES
C        NBSOFA : NOMBRE DE SOMMETS DE CHACUNE DE SES FACES
C        NOSOFA : NO ELEMENTAIRE DES SOMMETS DE CHACUNE DE SES FACES
C        NOTYFA : NO DU TYPE       DE CHACUNE DE SES FACES
C        NBPOFA : NOMBRE DE POINTS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C        NOPOFA : NO ELEMENTAIRE DES POINTS-NON SUR LES ARETES DE CHACUNE DE
C                 SES FACES
C        NBNOFA : NOMBRE DE NOEUDS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C        NONOFA : NO ELEMENTAIRE DES NOEUDS-NON SUR LES ARETES DE CHACUNE DE
C                 SES FACES
C
C        RESERVATION DU TABLEAU NPEF"TYPE_EF
C        ICI SOMMETS=NOEUDS=POINTS
         NBSO = NBSOME( NCOGEL )
C        LE TABLEAU DES NUMEROS DES SOMMETS
         CALL TNMCDC( 'ENTIER', NBSO * NBELT, MNSOEL )
C
         NBEL   = 0
C        LE TABLEAU DES POINTS  ( NO ELEMENT, NO SOMMET, NO POINT )
         MXPO   = 0
         NBPO   = 0
         MNPOEL = 0
C        LE TABLEAU DES LIGNES DES ARETES ( NO ELEMENT, NO ARETE, NO LIGNE )
         MXLI   = 0
         NBLI   = 0
         MNLIEL = 0
C        LE TABLEAU DES SURFACES DES FACES (NO ELEMENT, NO FACE, NO SURFACE)
         MXSU   = 0
         NBSU   = 0
         MNSUEL = 0
C        LE TABLEAU DES VOLUMES DES CUBES (NO ELEMENT, NO CUBE, NO VOLUME )
         MXVO   = 0
         NBVO   = 0
         MNVOEL = 0
C
C        BOUCLE SUR LES ELEMENTS FINIS DE TYPE NCOGEL
C        ********************************************
         MOEL = MOELEM( NCOGEL )
         MNEL = MNELEM( NCOGEL ) - 1 - MOEL
         MNSE = MNSOEL - 1 - NBSO
C
         DO 100 I=1,NBELT
C
C           TOUT ELEMENT FINI DEJA VU OU INEXISTANT EST SAUTE
            MNEL = MNEL + MOEL
            IF( MCN( MNEL + 1 ) .LE. 0 ) GOTO 100
C           CET EF EST INITIALISE
            NOBJ = MCN( MNEL + MOEL - 1 )
C           LE NUMERO DE L'OBJET PREMIER QUI LE CONTIENT
            IF( NOBJ .LE. 0 ) GOTO 100
C
C           L'ELEMENT FINI EXISTE : STOCKAGE DU NUMERO DE SES SOMMETS
            NBEL = NBEL + 1
            MNSE = MNSE + NBSO
            DO 12 J=1,NBSO
               MCN( MNSE + J ) = MCN( MNEL + J )
 12         CONTINUE
C           L'ELEMENT FINI EST MARQUE
            MCN( MNEL + MOEL - 1 ) = -NOBJ
C
C           L'ELEMENT FINI EST IL A TG?
            NUEFTG = MCN( MNEL + MOEL )
C           L'ADRESSE MCN DE LA PREMIERE TG DE L'EF A TG
            IF( NUEFTG .GT. 0 ) THEN
               MNEFTG = MNELTG(NCOGEL) + (NUEFTG-1) * MOELTG(NCOGEL)
            ELSE
               MNEFTG = 0
            ENDIF
C
C           LE NUMERO EVENTUEL DE POINT UTILISATEUR DES SOMMETS DE L'EF
C           ===========================================================
            IF( NBELEM(1) .GT. 0 ) THEN
               DO 20 J=1,NBSO
C                 RECHERCHE DE CE SOMMET DANS LE HACHAGE DES POINTS
                  NOSO(1) = MCN(MNEL+J)
                  CALL HACHAR( 1, NOSO, MOELEM(1), NBELEM(1),
     %                         MCN(MNELEM(1)), LIENEL(1), N )
                  IF( N .GT. 0 ) THEN
C
C                    SOMMET RETROUVE
C                    M ADRESSE DANS MCN DE +-NUMERO DU POINT DU SOMMET
                     M = MNELEM(1) - 2 + MOELEM(1) * N
C                    UN ELEMENT FINI POINT DE MOINS SI CE N'ETAIT DEJA FAIT
                     IF( MCN( M ) .GE. 0 ) THEN
                        NBEFPU(1) = NBEFPU(1) - 1
C                       LE POINT EST RETROUVE
                        IF( NUOBPR(2,MCN(M)) .GT. 0 ) THEN
                           CALL NMOBNU('POINT',NUOBPR(2,MCN(M)),NMPLSV)
                           NUOBPR(2,MCN(M)) = -NUOBPR(2,MCN(M))
                           IF( LANGAG .EQ. 0 ) THEN
                              WRITE(IMPRIM,*) 'LE POINT   ',NMPLSV,
     %                                       ' EST UN SOMMET'
                           ELSE
                              WRITE(IMPRIM,*)'THE POINT   ',NMPLSV,
     %                                       ' IS A VERTEX'
                           ENDIF
                        ENDIF
                     ENDIF
C                    LE SOMMET EST UN POINT UTILISATEUR .
C                    STOCKAGE DE SON NUMERO DANS LES OBJETS PREMIERS
                     NN = ABS( MCN( M ) )
                     CALL DFTOPA( I, J, NN, NBELEM(1),
     %                            NBPO, MXPO, MNPOEL )
C                    LE MARQUAGE DU SOMMET POINT
                     MCN( M ) = -NN
C                    LE POINT EST MARQUE 'INTERNE' OU 'AUX LIMITES'
                     IF( NCOGEL .EQ. 1 ) THEN
                        NUOBIN( NN ) = 1
                     ELSE
                        NUOBCL( NN ) = 1
                     ENDIF
C
C                    PAS DE RECHERCHE DE TG CAR SI LE POINT A 1 OU 2 OU 3 TG
C                    ON NE SAIT PAS DANS QUELLE DIRECTION DANS L'EF LES AFFECTER
C
                  ENDIF
 20            CONTINUE
               IF( NCOGEL .LE. 1 ) GOTO 100
            ENDIF
C
C           LE NUMERO DES SOMMETS DES ARETES ET FACES DE L'ELEMENT
            CALL SOARFA( NCOGEL, NOSOAR, NOSOFA )
C
C           LE NUMERO EVENTUEL DE LIGNE DES ARETES DE L'ELEMENT FINI
C           ========================================================
            NBAR = NBARET( NCOGEL )
            IF( NBELEM(2) .GT. 0 ) THEN
               DO 30 J=1,NBAR
C                 LE NUMERO CROISSANT DES SOMMETS DE L'ARETE J
                  IF( NBAR .EQ. 1 ) THEN
                     NOSO(1) = 1
                     NOSO(2) = 2
                  ELSE
                     NOSO(1) = NOSOAR(1,J)
                     NOSO(2) = NOSOAR(2,J)
                  ENDIF
C                 LE NUMERO GLOBAL DES 2 SOMMETS DE L'ARETE
                  NOST(1) = MCN( MNEL + NOSO(1) )
                  NOST(2) = MCN( MNEL + NOSO(2) )
C                 TRI CROISSANT
                  CALL HACREN( 2, 2, 0, NOST )
C                 RECHERCHE DE CETTE ARETE DANS LE HACHAGE DES ARETES
                  CALL HACHAR( 2, NOST, MOELEM(2), NBELEM(2),
     %                         MCN(MNELEM(2)), LIENEL(2), N )
                  IF( N .GT. 0 ) THEN
C                    ARETE RETROUVEE
                     M = MNELEM(2) - 2 + MOELEM(2) * N
C                    ADRESSE MCN DE +-NUMERO DE LA LIGNE UTILISATEUR
C                    UN SEGMENT DE MOINS SI CE N'ETAIT DEJA FAIT
                     IF( MCN( M ) .GE. 0 ) THEN
                        NBEFPU(2) = NBEFPU(2) - 1
C                       UNE ARETE DE LIGNE EST RETROUVEE
                        IF( NUOBPR(2,MCN(M)) .GT. 0 ) THEN
                           CALL NMOBNU('LIGNE',NUOBPR(2,MCN(M)),NMPLSV)
                           NUOBPR(2,MCN(M)) = -NUOBPR(2,MCN(M))
                           IF( LANGAG .EQ. 0 ) THEN
                              WRITE(IMPRIM,*) 'LA LIGNE   ',NMPLSV,
     %                                       ' EST RETROUVEE'
                           ELSE
                              WRITE(IMPRIM,*)'THE LINE    ',NMPLSV,
     %                                       ' IS FOUND'
                           ENDIF
                        ENDIF
                     ENDIF
C                    LE NUMERO DE LA LIGNE UTILISATEUR DANS LES OBJETS PREMIERS
                     NN = ABS( MCN( M ) )
C                    L'ARETE LIGNE EST EMPILEE
                     CALL DFTOPA( I, J, NN, NBELEM(2),
     %                            NBLI, MXLI, MNLIEL )
C                    LE MARQUAGE DE L'ARETE
                     MCN( M ) = -NN
C
C                    LA LIGNE EST MARQUEE 'INTERNE' OU 'AUX LIMITES'
                     IF( NCOGEL .EQ. 2 ) THEN
                        NUOBIN( NN ) = 1
                     ELSE
                        NUOBCL( NN ) = 1
                     ENDIF
C
C                    LE NUMERO D'EF A TG DE L'ARETE J
                     NEFTG = MCN( M+1 )
                     IF( NBAR .GT. 1 .AND. NEFTG .GT. 0 ) THEN
C                       L'EF EST AU MOINS SURFACIQUE
C                       LE NUMERO DES 2 TGS DE L'ARETE
                        MM = MNELTG(2) + (NEFTG-1)*MOELTG(2)
                        NTG(1) = MCN( MM )
                        NTG(2) = MCN( MM+1 )
C                       SI LE HACHAGE A RENVERSE LES SOMMETS
C                       ALORS LES TANGENTES DOIVENT L'ETRE AUSSI
                        IF( NOST(2) .EQ. MCN( MNEL + NOSO(1) ) ) THEN
                           L      = NTG(1)
                           NTG(1) = NTG(2)
                           NTG(2) = L
                        ENDIF
C                       LE NO DE CES TANGENTES DE L'ARETE J DANS L'EF DE TYPE NC
                        CALL TGAREF( NCOGEL, J, NUTG )
                        DO 28 K=1,2
                           IF( NTG(K) .NE. 0 ) THEN
C                             LA TANGENTE EXISTE SUR L'ARETE
                              IF( MNEFTG .GT. 0 ) THEN
C                                IL EXISTE DES TG POUR CET EF DE TYPE NCOGEL
C                                LA TG DE L'ARETE EST IMPOSEE A L'ARETE DE L'EF
C                                AFIN DE FAVORISER LA CONTINUITE C1
                                 MCN( MNEFTG-1+NUTG(K) ) = NTG(K)
                              ELSE
C                               L'EF N'ETAIT PAS A TG => IL LE DEVIENT
                                IF(NUELTG(NCOGEL).GE.NBELTG(NCOGEL))THEN
C                                   LE TABLEAU TROP PETIT EST AUGMENTE
                                    L = NBELTG(NCOGEL)*MOELTG(NCOGEL)
                                    CALL TNMCAU( 'ENTIER', L,
     %                              (NBELTG(NCOGEL)+128)*MOELTG(NCOGEL),
     %                               L, MNELTG(NCOGEL) )
                                    NBELTG(NCOGEL)=NBELTG(NCOGEL)+128
                                 ENDIF
C                                L'ADRESSE DE L'EF A TG DE TYPE NCOGEL
                                 MNEFTG = MNELTG(NCOGEL)+
     %                                    NUELTG(NCOGEL)*MOELTG(NCOGEL)
                                 NUELTG(NCOGEL) = NUELTG(NCOGEL) + 1
                                 CALL AZEROI(MOELTG(NCOGEL),MCN(MNEFTG))
C                                LA TG DE L'ARETE EST IMPOSEE A L'ARETE DE L'EF
                                 MCN( MNEFTG-1+NUTG(K) ) = NTG(K)
C                                L'EF DE TYPE NCOGEL DEVIENT UN EF A TG
                                 MCN( MNEL + MOEL ) = NUELTG(NCOGEL)
                              ENDIF
                           ENDIF
 28                     CONTINUE
                     ENDIF
                  ENDIF
 30            CONTINUE
               IF( NCOGEL .LE. 2 ) GOTO 100
            ENDIF
C
C           LE NUMERO EVENTUEL DE SURFACE DES FACES DE L'ELEMENT FINI
C           =========================================================
            IF( NBELEM(3) .GT. 0 .OR. NBELEM(4) .GT. 0 ) THEN
               NBFA = NBFACE( NCOGEL )
               DO 50 J=1,NBFA
                  IF( NBFA .GT. 1 ) THEN
C                    L'EF DE TYPE NCOGEL EST VOLUMIQUE
C                    LE NOMBRE DE SOMMETS DE LA FACE
                     NB = NBSOFA(J)
                     IF( NOTYFA(J) .LE. 7 ) THEN
C                       TRIANGLE P1 OU P2  => CODE GEOMETRIQUE 3
                        IF( NBELEM(3) .LE. 0 ) GOTO 50
                        NC = 3
                     ELSE
C                       QUADRANGLE Q1 OU Q2  => CODE GEOMETRIQUE 4
                        IF( NBELEM(4) .LE. 0 ) GOTO 50
                        NC = 4
                     ENDIF
C                    LE NUMERO CROISSANT DES NB SOMMETS DE LA FACE J
                     DO 40 K=1,NB
                        NOSO(K) = MCN( MNEL + NOSOFA(K,J) )
 40                  CONTINUE
C                    AJOUT DU NUMERO POSITION DES TGS DE LA FACE DANS L'EF VOLUM
                     DO 42 K=NB+1,NB*3
                        NOSO(K) = NOTGFA(K-NB,J,NCOGEL)
 42                  CONTINUE
C                    POSITION DANS NOSO DE LA PREMIERE TG
                     NU1TG = NB+1
                  ELSE
C                    L'EF DE TYPE NCOGEL EST SURFACIQUE:TRIANGLE OU QUADRANGLE
                     NB = NBSO
                     NC = NCOGEL
                     DO 45 K=1,NB
                        NOSO(K) = MCN( MNEL + K )
 45                  CONTINUE
C                    POSITION DANS NOSO DE LA PREMIERE TG: ICI PAS DE TG
                     NU1TG = 0
                  ENDIF
C
C                 RECHERCHE DE CETTE FACE DANS LE HACHAGE DES FACES
                  CALL HACHTQ( NB, NOSO, NU1TG,
     %                         MOELEM(NC),      NBELEM(NC),
     %                         MCN(MNELEM(NC)), LIENEL(NC), N )
                  IF( N .GT. 0 ) THEN
C                    FACE RETROUVEE
                     M = MNELEM(NC) - 2 + MOELEM(NC) * N
C                    ADRESSE MCN DE +-NUMERO DE LA SURFACE UTILISATEUR
C                    UNE FACE DE MOINS SI CE N'ETAIT DEJA FAIT
                     IF( MCN( M ) .GE. 0 ) THEN
                        NBEFPU(NC) = NBEFPU(NC) - 1
C                       UNE FACE DE SURFACE EST RETROUVEE
                        IF( NUOBPR(2,MCN(M)) .GT. 0 ) THEN
                          CALL NMOBNU('SURFACE',NUOBPR(2,MCN(M)),NMPLSV)
                           NUOBPR(2,MCN(M)) = -NUOBPR(2,MCN(M))
                           IF( LANGAG .EQ. 0 ) THEN
                              WRITE(IMPRIM,*) 'LA SURFACE ',NMPLSV,
     %                                        ' EST RETROUVEE'
                           ELSE
                              WRITE(IMPRIM,*)'THE SURFACE ',NMPLSV,
     %                                       ' IS FOUND'
                           ENDIF
                        ENDIF
                     ENDIF
C                    LE NUMERO DE LA SURFACE UTILISATEUR DANS LES OBJETS PREMIER
                     NN = ABS( MCN( M ) )
C                    LA FACE SURFACE EST EMPILEE
                     CALL DFTOPA( I, J, NN, NBELEM(NC),
     %                            NBSU, MXSU, MNSUEL )
C                    LE MARQUAGE DE LA FACE
                     MCN( M ) = -NN
C
C                    LA SURFACE EST MARQUEE 'INTERNE' OU 'AUX LIMITES'
                     IF( NCOGEL .EQ. 3 .OR. NCOGEL .EQ. 4 ) THEN
                        NUOBIN( NN ) = 1
                     ELSE
                        NUOBCL( NN ) = 1
                     ENDIF
C
C                    LE NUMERO D'EF A TG DE LA FACE
                     NEFTG = MCN( M+1 )
                     IF( NEFTG .GT. 0 .AND. NBFA .GT. 1 ) THEN
C                       LA FACE A TG EXISTE ET L'EF EST VOLUMIQUE
C                       L'ADRESSE MCN DU NUMERO DES MOELTG(NC) TGS DE LA FACE
                        MM = MNELTG(NC) + (NEFTG-1)*MOELTG(NC) -1
C                       NOSO(NB+1:NB+2*NB) CONTIENT LA POSITION DES TGS DANS L'E
                        DO 48 K=1,MOELTG(NC)
                           IF( MCN(MM+K) .NE. 0 ) THEN
C                             LA TANGENTE EXISTE SUR LA FACE
                              IF( MNEFTG .GT. 0 ) THEN
C                                IL EXISTE DES TG POUR CET EF DE TYPE NCOGEL
C                                LA TG DE LA FACE EST IMPOSEE A LA FACE DE L'EF
C                                AFIN DE FAVORISER LA CONTINUITE C1
                                 MCN( MNEFTG-1+NOSO(NB+K) ) = MCN(MM+K)
                              ELSE
C                               L'EF N'ETAIT PAS A TG => IL LE DEVIENT
                                IF(NUELTG(NCOGEL).GE.NBELTG(NCOGEL))THEN
C                                   LE TABLEAU TROP PETIT EST AUGMENTE
                                    L = NBELTG(NCOGEL)*MOELTG(NCOGEL)
                                    CALL TNMCAU( 'ENTIER', L,
     %                              (NBELTG(NCOGEL)+128)*MOELTG(NCOGEL),
     %                               L, MNELTG(NCOGEL) )
                                    NBELTG(NCOGEL)=NBELTG(NCOGEL)+128
                                 ENDIF
C                                L'ADRESSE DE L'EF A TG DE TYPE NCOGEL
                                 MNEFTG = MNELTG(NCOGEL)+
     %                                    NUELTG(NCOGEL)*MOELTG(NCOGEL)
                                 NUELTG(NCOGEL) = NUELTG(NCOGEL) + 1
                                 CALL AZEROI(MOELTG(NCOGEL),MCN(MNEFTG))
C                                LA TG DE L'ARETE EST IMPOSEE A L'ARETE DE L'EF
                                 MCN( MNEFTG-1+NOSO(NB+K) ) = MCN(MM+K)
C                                L'EF DE TYPE NCOGEL DEVIENT UN EF A TG
                                 MCN( MNEL + MOEL ) = NUELTG(NCOGEL)
                              ENDIF
                           ENDIF
 48                     CONTINUE
                     ENDIF
                  ENDIF
 50            CONTINUE
               IF( NCOGEL .LE. 4 ) GOTO 100
            ENDIF
C
C           LE NUMERO EVENTUEL DE VOLUME DU CUBE DE L'ELEMENT FINI
C           ======================================================
C           LE NUMERO CROISSANT DES SOMMETS DU CUBE J
            IF( NBELEM(5) .GT. 0 .OR. NBELEM(6) .GT. 0 .OR.
     %          NBELEM(7) .GT. 0 .OR. NBELEM(9) .GT. 0 ) THEN
C              LES SOMMETS ONT DEJA ETE RENUMEROTES POUR LE HACHAGE
C              RECHERCHE DU NUMERO DE CE CUBE DANS LE HACHAGE DES CUBES
               M = MNEL + 1
               CALL HACHAR( NBSO,MCN(M),MOELEM(NCOGEL),NBELEM(NCOGEL),
     %                      MCN(MNELEM(NCOGEL)), LIENEL(NCOGEL), N )
               IF( N .GT. 0 ) THEN
C                 CUBE RETROUVE
                  M = MNELEM(NCOGEL) - 2 + MOELEM(NCOGEL) * N
C                 ADRESSE MCN DE +-NUMERO DU VOLUME UTILISATEUR
C                 UN CUBE DE MOINS SI CE N'ETAIT DEJA FAIT
                  IF( MCN( M ) .GE. 0 ) THEN
                     NBEFPU(NCOGEL) = NBEFPU(NCOGEL) - 1
C                    UN EF DE VOLUME EST RETROUVE
                     IF( NUOBPR(2,MCN(M)) .GT. 0 ) THEN
                        CALL NMOBNU( 'VOLUME',NUOBPR(2,MCN(M)),NMPLSV )
                        NUOBPR(2,MCN(M)) = -NUOBPR(2,MCN(M))
                        IF( LANGAG .EQ. 0 ) THEN
                           WRITE(IMPRIM,*) 'LE VOLUME  ',NMPLSV,
     %                                     ' EST RETROUVE'
                        ELSE
                           WRITE(IMPRIM,*) 'THE VOLUME  ',NMPLSV,
     %                                     ' IS FOUND'
                        ENDIF
                     ENDIF
                  ENDIF
C                 LE NUMERO DU VOLUME UTILISATEUR DANS LES OBJETS PREMIERS
                  NN = ABS( MCN( M ) )
C                 LE CUBE VOLUME EST AJOUTE
                  IF( MXVO .EQ. 0 ) THEN
                     MXVO = NBELEM(NCOGEL)
                     CALL TNMCDC( 'ENTIER', MXVO, MNVOEL )
                  ENDIF
C                 LE NUMERO DU VOLUME
                  MCN( MNVOEL + NBVO ) = NN
                  NBVO = NBVO + 1
C                 LE MARQUAGE DU CUBE
                  MCN( M ) = -NN
C                 LE VOLUME EST MARQUE 'INTERNE'
                  NUOBIN( NN ) = 1
               ELSE
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'ELEMENT FINI VOLUMIQUE NON RETROUVE'
                     WRITE(IMPRIM,10060) ( MCN(MNEL+K), K=1,NBSO )
                  ELSE
                     KERR(1) = 'FINITE ELEMENT 3D NOT FOUND'
                     WRITE(IMPRIM,20060) ( MCN(MNEL+K), K=1,NBSO )
                  ENDIF
                  CALL LEREUR
10060    FORMAT(' EF VOLUMIQUE NON RETROUVE DE SOMMETS:',8I7)
20060    FORMAT(' FINITE ELEMENT NOT FOUND of VERTICES:',8I7)
                  IERR = 3
               ENDIF
            ENDIF
 100     CONTINUE
C
C        GENERATION DU TABLEAU ~>OBJET>>NPEF"NCOGEL
C        ==========================================
C        LES VARIABLES DE CE TABLEAU
C        NBNOE  : NOMBRE DE NOEUDS DE L ELEMENT DE NOM NOMELE
         NBNDEL = NBNOE
C        NBPOE  : NOMBRE DE POINTS DE L ELEMENT DE NOM NOMELE
         NBPGEL = NBPOE
C        NOTRAE : CODE TRAITEMENT DE L ELEMENT NOMELE
C                 0 : NOEUDS=POINTS=SOMMETS
C                 1 : NOEUDS=POINTS#SOMMETS
C                 2 : NOEUDS#POINTS=SOMMETS
C                 3 : NOEUDS#POINTS#SOMMETS
         IF( NOTRAE .EQ. 0 .OR. NOTRAE .EQ. 1 ) THEN
C           PAS DE TABLEAU DES POINTS GEOMETRIQUES = NOEUDS
            NBPGEL = 0
         ENDIF
C
C        POINT-SOMMET
         NBPSEL = NBPO
         IF( NCOGEL .EQ. 1 ) THEN
            MOPSEL = 0
         ELSE
            MOPSEL = NBPO
         ENDIF
C
C        LIGNE ARETE
         NBLAEL = NBLI
         IF( NCOGEL .EQ. 2 ) THEN
            MOLAEL = 0
         ELSE
            MOLAEL = NBLI
         ENDIF
C
C        SURFACE FACE
         NBSFEL = NBSU
         IF( NCOGEL .EQ. 3 .OR. NCOGEL .EQ. 4 ) THEN
            MOSFEL = 0
         ELSE
            MOSFEL = NBSU
         ENDIF
C
C        VOLUME CUBE
         NBVCEL = NBVO
         IF( NCOGEL .GT. 4 ) THEN
            MOVCEL = 0
         ELSE
            MOVCEL = NBVO
         ENDIF
C
C        LES TG DANS NPEF
         IF( NUELTG(NCOGEL) .GT. 0 ) THEN
            NBELAP = NBELEM( NCOGEL )
            NBEFTG = NUELTG( NCOGEL )
            NBTGEL = MOELTG( NCOGEL )
         ELSE
            NBELAP = 0
            NBEFTG = 0
            NBTGEL = 0
         ENDIF
C
C        LE NOMBRE DE MOTS DU TABLEAU 'NPEF'
         K = WOVCEL + 1 + NBEL   * ( NBNDEL + NBPGEL )
     %                  + NBPSEL + MOPSEL + MOPSEL
     %                  + NBLAEL + MOLAEL + MOLAEL
     %                  + NBSFEL + MOSFEL + MOSFEL
     %                  + NBVCEL + MOVCEL + MOVCEL
     %                  + NBELAP + NBEFTG * NBTGEL
         CALL LXTSOU( NTLXOB, KNOM, NTELE, MNELE )
         IF( NTELE .GT. 0 ) THEN
C           DESTRUCTION DE L'ANCIEN TABLEAU
            CALL LXTSDS( NTLXOB, KNOM )
         ENDIF
         CALL LXTNDC( NTLXOB, KNOM, 'MOTS', K )
         CALL LXTSOU( NTLXOB, KNOM, NTELE, MNELE )
         NTNPEF( NCOGEL ) = NTELE
         MNNPEF( NCOGEL ) = MNELE
C        NOM DU TYPE DES ELEMENTS FINIS
         NMNPEF( NCOGEL ) = ICHARX( NOM4 )
C
C        REMPLISSAGE DES VARIABLES
C        NUMERO DU TYPE DE L'ELEMENT DE CODE GEOMETRIQUE NCOGEL
         MCN( MNELE + WUTYEL ) = NUTYEL
C        NBELEM 'NOMBRE D'EF DE CE TYPE
         MCN( MNELE + WBELEM ) = NBEL
C        NBNDEL 'NOMBRE DE NOEUDS D''UN ELEMENT FINI DE CE TYPE'
         MCN( MNELE + WBNDEL ) = NBNDEL
C        NBPGEL 'NOMBRE DE POINTS GEOMETRIQUES D''UN ELEMENT'
         MCN( MNELE + WBPGEL ) = NBPGEL
C        NBPSEL 'NOMBRE DE POINTS UTILISATEUR SOMMET D''UN ELEMENT'
         MCN( MNELE + WBPSEL ) = NBPSEL
C        MOPSEL 'NOMBRE DE VARIABLES DES TABLEAUX NLPSEL ET NEPSEL'
         MCN( MNELE + WOPSEL ) = MOPSEL
C        NBLAEL 'NOMBRE DE LIGNES UTILISATEUR ARETE D''UN ELEMENT'
         MCN( MNELE + WBLAEL ) = NBLAEL
C        MOLAEL 'NOMBRE DE VARIABLES DES TABLEAUX NLLAEL ET NELAEL'
         MCN( MNELE + WOLAEL ) = MOLAEL
C        NBSFEL 'NOMBRE DE SURFACES UTILISATEUR FACE D''UN ELEMENT'
         MCN( MNELE + WBSFEL ) = NBSFEL
C        MOSFEL 'NOMBRE DE VARIABLES DES TABLEAUX NLSFEL ET NESFEL'
         MCN( MNELE + WOSFEL ) = MOSFEL
C        NBVCEL 'NOMBRE DE VOLUMES UTILISATEUR CUBE D''UN ELEMENT'
         MCN( MNELE + WBVCEL ) = NBVCEL
C        MOVCEL 'NOMBRE DE VARIABLES DES TABLEAUX NLVCEL ET NEVCEL'
         MCN( MNELE + WOVCEL ) = MOVCEL
C        NOMBRE D'EF AVEC POINTEUR SUR EF A TG
         MCN( MNELE + WBELAP ) = NBELAP
C        NOMBRE D'EF A TG
         MCN( MNELE + WBELTG ) = NBEFTG
C        NOMBRE DE TANGENTES D'UN EF
         MCN( MNELE + WBTGEL ) = NBTGEL
C
C        LE REMPLISSAGE DES TABLEAUX
C        NUNDEL(1..NBELEM,1..NBNDEL) 'NUMERO DES NOEUDS DE CHAQUE EF'
         MN1 = MNSOEL
         MN2 = MNELE + WUNDEL
C        SOEL RANGE PAR ELEMENT
C        NDEL RANGE PAR NUMERO LOCAL DES NOEUDS
         DO 510 I=0,NBEL-1
            MN3 = MN2 + I
C           POUR L'INSTANT SEULS LES SOMMETS SONT TRAITES
            DO 500 N=1,NBSO
               MCN( MN3 ) = MCN( MN1 )
               MN1 = MN1 + 1
               MN3 = MN3 + NBEL
 500        CONTINUE
 510     CONTINUE
C
C        NUPGEL(1..NBELEM,1..NBPGEL)  EST VIDE ICI
C        NUMERO DES POINTS GEOMETRIQUES DE CHAQUE ELEMENT FINI
C
C
C        NUPSEL(1..NBPSEL) 'NOM DU POINT UTILISATEUR SOMMET D''UN ELEMENT'
C        NLPSEL(1..MOPSEL) 'NUMERO LOCAL A L''ELEMENT DU SOMMET POINT'
C        NEPSEL(1..MOPSEL) 'NUMERO DE L''ELEMENT DU SOMMET POINT'
         MN  = MNPOEL
         MN1 = MNELE + WUNDEL + NBEL * ( NBNDEL + NBPGEL )
         MN2 = MN1 + NBPSEL
         MN3 = MN2 + MOPSEL
         DO 710 I=0,NBPSEL-1
C           LE NUMERO DU POINT-SOMMET
            MCN( MN1 + I ) = ABS( NUOBPR(2,ABS(MCN(MN))) )
            IF( MOPSEL .GT. 0 ) THEN
               MCN( MN2+I ) = MCN( MN+1 )
               MCN( MN3+I ) = MCN( MN+2 )
            ENDIF
            MN = MN + 3
 710     CONTINUE
C
C        NULAEL(1..NBLAEL) 'NOM DE LA LIGNE UTILISATEUR ARETE D''UN ELEMENT'
C        NLLAEL(1..MOLAEL) 'NUMERO LOCAL A L''ELEMENT DE LA LIGNE ARETE'
C        NELAEL(1..MOLAEL) 'NUMERO DE L''ELEMENT DE LA LIGNE ARETE'
         MN  = MNLIEL
         MN1 = MN3 + MOPSEL
         MN2 = MN1 + NBLAEL
         MN3 = MN2 + MOLAEL
         DO 720 I=0,NBLAEL-1
C           NUMERO DE LIGNE SUPPORT DE L ARETE
            MCN( MN1 + I ) = ABS( NUOBPR(2,ABS(MCN(MN))) )
            IF( MOLAEL .GT. 0 ) THEN
C              NUMERO LOCAL DE L'ARETE DANS L'EF MCN(MN+2)
               MCN( MN2+I ) = MCN( MN+1 )
C              NUMERO DE L'ELEMENT FINI GLOBAL DANS LE MAILLAGE DES EF
               MCN( MN3+I ) = MCN( MN+2 )
            ENDIF
            MN = MN + 3
 720     CONTINUE
C
C        NUSFEL(1..NBSFEL) 'NOM DE LA SURFACE UTILISATEUR FACE D''UN ELEMENT'
C        NLSFEL(1..MOSFEL) 'NUMERO LOCAL A L''ELEMENT DE LA SURFACE FACE'
C        NESFEL(1..MOSFEL) 'NUMERO DE L''ELEMENT DE LA SURFACE FACE'
         MN  = MNSUEL
         MN1 = MN3 + MOLAEL
         MN2 = MN1 + NBSFEL
         MN3 = MN2 + MOSFEL
         DO 730 I=0,NBSFEL-1
C           NUMERO DE SURFACE SUPPORT DE LA FACE
            MCN( MN1 + I ) = ABS( NUOBPR(2,ABS(MCN(MN))) )
            IF( MOSFEL .GT. 0 ) THEN
C              NUMERO LOCAL DE LA FACE DANS L'EF MCN(MN+2)
               MCN( MN2+I ) = MCN( MN+1 )
C              NUMERO DE L'ELEMENT FINI GLOBAL DANS LE MAILLAGE DES EF
               MCN( MN3+I ) = MCN( MN+2 )
            ENDIF
            MN = MN + 3
 730     CONTINUE
C
C        NUVCEL(1..NBVCEL) 'NOM DU VOLUME UTILISATEUR D''UN ELEMENT'
C        NLVCEL(1..MOVCEL) 'NUMERO LOCAL A L''ELEMENT DU VOLUME'
C        NEVCEL(1..MOVCEL) 'NUMERO DE L''ELEMENT DU VOLUME'
         MN  = MNVOEL
         MN1 = MN3 + MOSFEL
         DO 740 I=0,NBVCEL-1
C           LE NUMERO DU VOLUME
            MCN( MN1 + I ) = ABS( NUOBPR( 2, ABS( MCN(MN+I) ) ) )
C           LES 2 AUTRES INFORMATIONS EGALES A 1 ET NUMERO EF NE SONT
C           PAS STOCKEES CAR EVIDENTES
 740     CONTINUE
C
C        LPELTG(1..NBELAP) NUMERO>0 DE L'EF A TG SINON 0
C        NUTGEL(1..NBELTG,1..NBTGEL) +-NO DES TANGENTES DE L'EF A TG
         IF( NBELAP .GT. 0 ) THEN
C           LE TABLEAU LPELTG EXISTE
            MN1    = MN1 + NBVCEL - 1
            MN2    = MN1 + NBELAP
            MOEL   = MOELEM( NCOGEL )
            MNEL   = MNELEM( NCOGEL ) - MOEL
            NBEL   = 0
C
            DO 760 I=1,NBELT
C
C              TOUT ELEMENT FINI DEJA VU OU INEXISTANT EST SAUTE
               MNEL = MNEL + MOEL
               IF( MCN( MNEL ) .LE. 0 ) GOTO 760
C              CET EF EST INITIALISE
               NBEL = NBEL + 1
C              LE NUMERO DE L'EF A TG
               NUEFTG = MCN( MNEL + MOEL - 1 )
C              LE NUMERO POINTEUR SUR L'EF A TG
               MCN(MN1+NBEL) = NUEFTG
               IF( NUEFTG .GT. 0 ) THEN
C                 EF A TG : L'ADRESSE MCN DE L'EF A TG
                  MNEFTG = MNELTG(NCOGEL) + (NUEFTG-1) * NBTGEL - 1
C                 LES NUMEROS DES TANGENTES SONT AJOUTEES
                  MN3 = MN2 + NUEFTG
                  DO 750 J=1,NBTGEL
                     MCN(MN3) = MCN(MNEFTG+J)
C                    ATTENTION A L'ORDRE DE RANGEMENT: 1..NBELTG, 1..NBTGEL
                     MN3 = MN3 + NBEFTG
 750              CONTINUE
               ENDIF
 760        CONTINUE
         ENDIF
C
C        LA DATE
         CALL ECDATE( MCN(MNELE) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNELE + MOTVAR(6) ) = NONMTD( '~>>>NPEF' )
C
         NOMTMS(NCOGEL) = KNOM
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'dftop4: FIN DE CREATION DU TMS: ' // KNOM
         ELSE
            WRITE(IMPRIM,*) 'dftop4: END of CREATION of TMS: ' // KNOM
         ENDIF
C
C        DESTRUCTION DES TABLEAUX DEVENUS INUTILES
         CALL TNMCDS( 'ENTIER', NBSO*NBELT, MNSOEL )
         IF( MXPO .GT. 0 ) CALL TNMCDS( 'ENTIER', MXPO*3, MNPOEL )
         IF( MXLI .GT. 0 ) CALL TNMCDS( 'ENTIER', MXLI*3, MNLIEL )
         IF( MXSU .GT. 0 ) CALL TNMCDS( 'ENTIER', MXSU*3, MNSUEL )
         IF( MXVO .GT. 0 ) CALL TNMCDS( 'ENTIER', MXVO  , MNVOEL )
C
 1000 CONTINUE
      IF( IERR .NE. 0 ) GOTO 9000
C
C     REMISE A POSITIF DE NUOBPR(2,.)
      DO I=1,NBOBPR
         NUOBPR(2,I) = ABS( NUOBPR(2,I) )
      ENDDO
C
C     EXISTE-T-IL DES FACES  DES SURFACES NON RETROUVEES DANS LES EF?
C     EXISTE-T-IL DES ARETES DES LIGNES   NON RETROUVEES DANS LES EF?
C     EXISTE-T-IL DES SOMMETS DES POINTS  NON RETROUVES  DANS LES EF?
      IF( LDIMEF .EQ. 3 ) THEN
         I = 4
      ELSE IF( LDIMEF .EQ. 2 ) THEN
         I = 2
      ELSE
         I = 1
      ENDIF
      INTER0 = INTERA
      INTERA = 3
      DO 1055 N=I,1,-1
         IF( NBEFPU(N) .GT. 0 ) THEN
C           UN EF D'UN PLS N'A PAS ETE RETROUVE DANS LES EF DE L'OBJET
            IERR = IERR + N
            NBLGRC(NRERR) = 1
            WRITE(KERR(8)(1:10), '(I10)') NBEFPU(N)
            IF( LANGAG .EQ. 0 ) THEN
               IF( LDIMEF .EQ. 3 ) THEN
                  IF( N .EQ. 4 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' QUADRANGLES DES SURFACES OUBLIES DANS LES FACES DES EF 3D'
                  ELSE IF( N .EQ. 3 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' TRIANGLES DES SURFACES OUBLIES DANS LES FACES DES EF 3D'
                  ELSE IF( N .EQ. 2 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' ARETES DES LIGNES OUBLIEES DANS LES ARETES DES EF 3D'
                  ELSE IF( N .EQ. 1 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' POINTS OUBLIES DANS LES SOMMETS DES EF 3D'
                  ENDIF
               ELSE IF( LDIMEF .EQ. 2 ) THEN
                  IF( N .EQ. 2 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' ARETES DES LIGNES OUBLIEES DANS LES ARETES DES EF 2D'
                  ELSE IF( N .EQ. 1 ) THEN
                     KERR(1) = KERR(8)(1:10)
     % // ' POINTS OUBLIES DANS LES SOMMETS DES EF 2D'
                  ENDIF
               ELSE
                  IF( N .EQ. 1 ) THEN
                     KERR(1) = KERR(8)(1:10)
     % // ' POINTS OUBLIES DANS LES SOMMETS DES EF 1D'
                  ENDIF
               ENDIF
            ELSE
               IF( LDIMEF .EQ. 3 ) THEN
                  IF( N .EQ. 4 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN QUADRANGLES of SURFACES among the FACES of 3D FE'
                  ELSE IF( N .EQ. 3 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN TRIANGLES of SURFACES among the FACES of 3D FE'
                  ELSE IF( N .EQ. 2 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN EDGES of LINES among the EDGES of 3D FE'
                  ELSE IF( N .EQ. 1 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN POINTS among the VERTICES of 3D FE'
                  ENDIF
               ELSE IF( LDIMEF .EQ. 2 ) THEN
                  IF( N .EQ. 2 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN EDGES of LINES among the EDGES of 2D FE'
                  ELSE IF( N .EQ. 1 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN POINTS among the VERTICES of 2D FE'
                  ENDIF
               ELSE
                  IF( N .EQ. 1 ) THEN
                     KERR(1) = KERR(8)(1:10)
     %  // ' FORGOTTEN POINTS among the VERTICES of 1D FE'
                  ENDIF
               ENDIF
            ENDIF
            CALL LEREUR
         ENDIF
 1055 CONTINUE
      INTERA = INTER0
      IF( IERR .NE. 0 ) GOTO 9000
      RETURN

C     INCOMPATIBILITES DES FACES ARETES SOMMETS
 9000 DO 9010 NCOGEL=9,1,-1
         IF( NOMTMS(NCOGEL)(1:4) .NE. '    ' ) THEN
            CALL LXTSDS( NTLXOB, NOMTMS(NCOGEL) )
         ENDIF
 9010 CONTINUE

      RETURN
      END
