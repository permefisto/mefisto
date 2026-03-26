        SUBROUTINE VOEX16( NTLXVC, LADEFI, RADEFI,
     %                     NTFAVC, MNFAVC, NTSOVC, MNSOVC, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  GENERER LA PENTA-HEXAEDRISATION D'UN VOLUME
C -----  OBTENUE PAR COUCHES LIMITES ET AUTRES A PARTIR D'UNE SURFACE
C        ORIENTEE MAILLEE AVEC DES TRIANGLES et/ou QUADRANGLES
C        REMARQUE: EPAISSEUR1COUCHE<0 DONNE L'AUTRE SENS DE LA NORMALE!
C
C        ATTENTION: SI LA SURFACE N'EST PAS CONVEXE LES HEXAEDRES
C                   PEUVENT SE CHEVAUCHER
C                   TRACER POUR VERIFIER SI LE RESULTAT EST CORRECT
C ENTREES:
C --------
C NTLXVC : NUMERO DU TMS DU LEXIQUE DU VOLUME A CREER
C LADEFI : TABLEAU ENTIER DE DEFINITION DU VOLUME
C RADEFI : TABLEAU REEL   DE DEFINITION DU VOLUME
C          CF '~/td/d/a_volume__definition'
C
C SORTIES:
C --------
C NTFAVC : NUMERO      DU TMS 'NSEF' DES NUMEROS DES CUBES
C MNFAVC : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES CUBES
C          CF '~/td/d/a___nsef'
C NTSOVC : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOVC : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF '~/td/d/a___xyzsommet'
C IERR   : = 0 SI PAS D'ERREUR
C          <>0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTHOR : Dai-Ni Hsieh    National Taiwan University TAIPEI TAIWAN 1/2010
C MODIFS : Alain Perronnet National Taiwan University TAIPEI TAIWAN 1/2010
C23456...............................................................012
      REAL        VALTEM
      PARAMETER  (VALTEM=1.23E37)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           EFTYPE
      INTEGER           LADEFI(0:*), NOSOEL(4+8), MNST(8)
      EQUIVALENCE      (MNST1, MNST(1)), (MNST2, MNST(2)),
     %                 (MNST3, MNST(3)), (MNST4, MNST(4)),
     %                 (MNST5, MNST(5)), (MNST6, MNST(6)),
     %                 (MNST7, MNST(7)), (MNST8, MNST(8))
      REAL              RADEFI(0:*)
      DOUBLE PRECISION  VECN(3), RRR
      REAL              XYZTEMP(3,4)
      CHARACTER*24      KNOMSU
C
      IERR   = 0
      MNSOSU = 0
      MNNBTQ = 0
      MNVECN = 0
      NSTG1  = 0
      NSTG2  = 0
      NSTG3  = 0
      NSTG4  = 0
      NSTG5  = 0
      NSTG6  = 0
      NSTG7  = 0
      NSTG8  = 0
C
C     VERIFICATION DES DONNEES
C     data verification
C     number of first boundary layers
      NBCOLA = LADEFI( WBCOLA )
      IF( NBCOLA .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER),'(I12)') NBCOLA
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE INCORRECT de PREMIERES COUCHES LIMITES='
     %                 //KERR(MXLGER)(1:12)
         ELSE
            KERR(1) = 'INCORRECT NUMBER of FIRST BOUNDARY LAYERS='
     %                 //KERR(MXLGER)(1:12)
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     EPAISSEUR D'UNE DES PREMIERES COUCHES LIMITES
C     thickness of ONE first boundary layer
      EPCOLA = RADEFI( WPCOLA )
      IF( EPCOLA .EQ. 0.0 ) THEN
C        L'EPAISSEUR PEUT ETRE POSITIVE OU NEGATIVE
C        thickness may be positive or negative
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER),'(G15.7)') EPCOLA
         IF( LANGAG .EQ. 0 ) THEN
           KERR(1)='EPAISSEUR INCORRECTE des PREMIERES COUCHES LIMITES='
     %              //KERR(MXLGER)(1:15)
         ELSE
           KERR(1)='INCORRECT THICKNESS of FIRST BOUNDARY LAYERS='
     %              //KERR(MXLGER)(1:15)
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     NOMBRE DES DERNIERES COUCHES LIMITES
C     number of LAST layers
      NBLACO = LADEFI( WBLACO )
      IF( NBLACO .LT. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER),'(I12)') NBLACO
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE INCORRECT de DERNIERES COUCHES='
     %                 //KERR(MXLGER)(1:12)
         ELSE
            KERR(1) = 'INCORRECT NUMBER of LAST LAYERS='
     %                 //KERR(MXLGER)(1:12)
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     EPAISSEUR DE LA DERNIERE COUCHE LIMITE
C     thickness of the LAST boundary layer
      EPLACO = RADEFI( WPLACO )
      IF( NBLACO .GT. 0 .AND. EPLACO .EQ. 0.0 ) THEN
C        L'EPAISSEUR PEUT ETRE POSITIVE OU NEGATIVE
C        thickness may be positive or negative
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER),'(G15.7)') EPLACO
         IF( LANGAG .EQ. 0 ) THEN
           KERR(1)='EPAISSEUR INCORRECTE des DERNIERES COUCHES='
     %              //KERR(MXLGER)(1:15)
         ELSE
           KERR(1)='INCORRECT THICKNESS of LAST BOUNDARY LAYERS='
     %              //KERR(MXLGER)(1:15)
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LES 2 EPAISSEURS DOIVENT AVOIR LE MEME SIGNE
C     the 2 thicknesses must have the same sign
C     UN SIGNE - FORCE LE SIGNE - POUR LES 2 EPAISSEURS
C     the sign - OF ONE INDUCED the sign - for the 2 THICKNESSES
      IF( EPCOLA * EPLACO .LT. 0 ) THEN
         EPLACO = -ABS( EPLACO )
         EPCOLA = -ABS( EPCOLA )
      ENDIF
C
C     RESTAURATION DU MAILLAGE DE LA SURFACE
C     retrieve the mesh of the surface
C     --------------------------------------
C     name of the surface to be magnified
      NUSFAG = LADEFI( WUSFAG )
      IF( NUSFAG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM INCORRECT de la SURFACE A GROSSIR'
         ELSE
            KERR(1) = 'UNKNOWN NAME of the SURFACE TO BE MAGNIFIED'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     OUVERTURE DE LA SURFACE NUSFAG
C     open the surface NUSFAG
      CALL LXNLOU( NTSURF, NUSFAG, NT, MN )
      IF( NT .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM INCORRECT de la SURFACE A GROSSIR'
         ELSE
            KERR(1) = 'UNKNOWN NAME of the SURFACE TO BE MAGNIFIED'
         ENDIF
         IERR = 1
         RETURN
      ENDIF
C
C     LE NOM KNOMSU  DE LA SURFACE
C     KNOMSU the name of the surface
      CALL NMOBNU( 'SURFACE', NUSFAG, KNOMSU )
C
C     ASSURER L'ORIENTATION DU MAILLAGE DE LA SURFACE
C     PAR PARCOURS DES EF A TRAVERS LES ARETES COMMUNES
C     ET PERMUTATION DES SOMMETS 2-NBSTEF DU SECOND EF
C     D'UNE ORIENTATION DIFFERENTE DU PREMIER EF
      CALL SUORIENT( KNOMSU,
     %               NTTQSU, MNTQSU, NTSOSU, MNSOSU, IERR )
C
C     RESTAURATION DU TABLEAU XYZSOMMET
C     retrieve the array XYZSOMMET
      CALL LXTSOU( NT, 'XYZSOMMET', NTSOSU, MNSOSU )
      IF( NTSOSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMSU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE SANS XYZSOMMET'
         ELSE
            KERR(2) = 'SURFACE WITHOUT XYZSOMMET'
         ENDIF
         IERR = 2
         RETURN
      ENDIF
C     LE NOMBRE DE SES SOMMETS ET TANGENTES
C     the numbers of vertices and tangents
      NBSOSU = MCN( MNSOSU+WNBSOM )
      NBTGSU = MCN( MNSOSU+WNBTGS )
C
C     RESTAURATION DU TABLEAU NSEF
C     retrieve the array NSEF
      CALL LXTSOU( NT, 'NSEF', NTTQSU, MNTQSU )
      IF( NTTQSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMSU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE SANS NSEF'
         ELSE
            KERR(2) = 'SURFACE WITHOUT NSEF'
         ENDIF
         IERR = 3
         RETURN
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DE LA SURFACE A DEPLACER
      CALL NSEFPA( MCN(MNTQSU),
     %             NUTYSU, NBSOEL, NBSOEF, NBTGEL,
     %             LDAPEF, LDNGEF, LDTGEF, NBTQSU,
     %             NX  , NY  , NZ  ,
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE A GROSSIR INCORRECTE'
         ELSE
            KERR(1) = 'INCORRECT SURFACE TO BE MAGNIFIED'
         ENDIF
         CALL LEREUR
         IERR = 7
         RETURN
      ENDIF
C
C     NOMBRE TRIANGLE QUADRANGLE A TG DE LA SURFACE
C     number of triangles or quadrangles with tangents of the surface
      NBTQTG = MCN( MNTQSU + WBEFTG )
C
C     ADRESSE AVANT 1-ER NO DE SOMMET DE LA SURFACE
      MNEFSU = MNTQSU + WUSOEF - 1
C
C     ADRESSE AVANT 1-ER POINTEUR SUR EF A TG
      MNEFAPSU = MNTQSU + LDAPEF -1
C     ADRESSE AVANT 1-ER CODE GEOMETRIQUE DES EF A TG
      MNEFCGSU = MNTQSU + LDTGEF -1
C     ADRESSE AVANT NO DE TG DES ARETES A TG
      MNEFSUTG = MNTQSU + LDTGEF -1
C
C     LE TYPE DU MAILLAGE DU VOLUME
C     LE MAILLAGE DU VOLUME EST NON STRUCTURE
      NUTYVC = 0
C
C     LE NOMBRE DE SOMMETS DU VOLUME MAILLEE
C     number of vertices of the volume
      NBSOVC = NBSOSU * ( 1 + NBCOLA + NBLACO )
C     LE NOMBRE DE TANGENTES DU VOLUME EST CELUI DE LA SURFACE INITIALE
      NBTGVC = NBTGSU
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
C     construct the array 'XYZSOMMET' of the volume
C     ----------------------------------------------------
      CALL LXTNDC( NTLXVC, 'XYZSOMMET', 'MOTS',
     %             WYZSOM + 3 * ( NBSOVC + NBTGVC ) )
      CALL LXTSOU( NTLXVC, 'XYZSOMMET',  NTSOVC, MNSOVC )
C     MISE A VALEUR TEMOIN DE NON UTILISATION DES X DES SOMMETS DU VOLUME
C     set a flag value on the x-coordinate of vertices of the volume
C     to indicate unvisited or not
      MN = MNSOVC + WYZSOM
      DO 10 I=1, NBSOVC
         RMCN( MN ) = VALTEM
         MN = MN + 3
 10   CONTINUE
C
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
C     construct the array 'NSEF' of the volume
C     -----------------------------------------------
C     LE NOMBRE D'EF CUBES DU VOLUME MAILLE
C     number of cubes of the mesh of the volume
      NBFAVC = NBTQSU * ( NBCOLA + NBLACO )
C
C     DECLARATION DES NO DE SOMMETS DES QUADRANGLES DU VOLUME
C     LE NOMBRE DE QUADRANGLES A TG
      IF( NBTGSU .GT. 0 ) THEN
         NBTGEF = 24
C        MEME STRUCTURE DU VOLUME QUE CELLE DE LA SURFACE
         NBEFTG = NBTQTG
         NBEFAP = NBFAVC
      ELSE
         NBTGEF = 0
         NBEFTG = 0
         NBEFAP = 0
      ENDIF
      CALL LXTNDC( NTLXVC, 'NSEF', 'ENTIER',
     %             WUSOEF + 8*NBFAVC + NBEFAP + NBEFTG * (1+NBTGEF) )
      CALL LXTSOU( NTLXVC, 'NSEF',  NTFAVC, MNFAVC )
C
C     LA PREMIERE SURFACE DE CUBES=PREMIERE COUCHE LIMITE
C     the first surface of cubes = first boundary layer
C     LES COORDONNEES DES SOMMETS ET TANGENTES DES CUBES
C     coordinates of vertices and tangents of cubes
      MNSTVC = MNSOVC + WYZSOM - 1
C
C     LA PREMIERE SURFACE DE SOMMET EST LA SURFACE INITIALE
C     vertices of the first surface are vertices of the original surface
      MN = MNSOSU + WYZSOM - 1
      DO 20 I = 1, 3*NBSOSU
         RMCN( MNSTVC + I ) = RMCN( MN + I )
 20   CONTINUE
C
C     LA SEULE SURFACE A TANGENTES DES CUBES EST CELLE DE LA SURFACE INITIALE
      MNVCTG = MNSOVC + WYZSOM + 3 * NBSOVC - 1
      MN     = MN + 3 * NBSOSU
      DO 30 I = 1, 3*NBTGSU
         RMCN( MNVCTG + I ) = RMCN( MN + I )
 30   CONTINUE
C
C     LES ADRESSAGES DANS MCN DANS LE TMS NSEF
      MNEFVC   = MNFAVC + WUSOEF - 1
      MNEFAPVC = MNEFVC + 8*NBFAVC
      MNEFVCTG = MNEFAPVC + NBEFAP + NBEFTG
C
C     TRAITEMENT DES QUADRANGLES A TG = CEUX DE LA PREMIERE COUCHE
      IF( NBTGSU .GT. 0 ) THEN
C
C        POINTEUR SUR LES PREMIERS CUBES A TG
         DO 40 I=1,NBTQTG
            MCN(MNEFAPVC+I) = MCN(MNEFAPSU+I)
 40      CONTINUE
C
C        POINTEUR NUL SUR LES CUBES DES COUCHES 2 3 ...
         CALL AZEROI( NBEFAP-NBTQTG, MCN(MNEFAPVC+NBTQTG+1) )
C
C        CODE GEOMETRIQUE C1 Degre 3  ICI zero
         CALL AZEROI( NBEFTG, MCN(MNEFAPVC+NBEFAP+1) )
C
C        MISE A ZERO DU NUMERO DES TANGENTES
         CALL AZEROI( NBTGEF*NBEFTG, MCN(MNEFVCTG+1) )
      ENDIF
C
C     define and initialize the counting array
      CALL TNMCDC('ENTIER', NBSOSU, MNNBTQ)
      DO 50 I = 1, NBSOSU
         MCN(MNNBTQ - 1 + I) = 0
 50   CONTINUE
C
      CALL TNMCDC('REEL', 3*NBSOSU, MNVECN)
C
C     CONSTRUCTION DES CUBES ET DES SOMMETS DE LA COUCHE 1
C     construct cubes and vertices of layer 1
      MNSTVC = MNSOVC + WYZSOM - 4
      DO 60 I = 1, NBTQSU
C
C        LE NUMERO DES 4 SOMMETS ET DES 4 TANGENTES DU TRIANGLE QUADRANGLE I
C        numbers of 4 vertices and numbers of 4 tangents of triangle or quadrang
         CALL NSEFNS( I    ,  NUTYSU, NBSOEF, NBTGEL,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNTQSU, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C
C        NCOGEL: NUMERO DU CODE GEOMETRIQUE DE L'EF NUELEM
C                2:SEGMENT
C        NUGEEF: >=0 NUMERO DE DEFINITION GEOMETRIQUE DE L'EF
C        NUEFTG: >0  NUMERO DE CET EF NUELEM PARMI LES EF A TG
C                =0  SI L'EF I N'A PAS DE TG
C        NOSOEL:  NUMERO DES 2 SOMMETS
C                 SUIVI EVENTUELLEMENT DES NBTGEF +-TANGENTES
C        LE - INDIQUE QUE LES COMPOSANTES DE LA TANGENTES SONT A INVERSER
C        SI PAS DE TANGENTES (NBTGEF=0):
C           ARETE: NO SOMMET1 , NS2
C        S'IL EXISTE DES TANGENTES (NBTGEF>0) ELLES SONT RANGEES PAR SOMMET:
C           ARETE: NO SOMMET1, NS2,
C                  NO TANGENTE1(S1S2), NT2(S2S1)
C
C        a value indicates the type of the finite element
C        zero -> triangle -> pentahedron, nonzero -> quadrangle -> hexahedron
         EFTYPE = NOSOEL(4)
C
         IF( EFTYPE .EQ. 0 ) THEN
C
C           TRIANGLE -> PENTAHEDRON
C
C           LES 3 NUMEROS DES SOMMETS DU TRIANGLE I
C           the 3 numbers of vertices of the triangle I
            IF( EPCOLA .GT. 0 ) THEN
               NST1 = NOSOEL(1)
               NST2 = NOSOEL(2)
               NST3 = NOSOEL(3)
            ELSE
               NST1 = NOSOEL(1)
               NST2 = NOSOEL(3)
               NST3 = NOSOEL(2)
            ENDIF
C
C           accumulate the counting of each vertex
            MCN( MNNBTQ-1 + NST1 ) = MCN( MNNBTQ-1 + NST1 ) + 1
            MCN( MNNBTQ-1 + NST2 ) = MCN( MNNBTQ-1 + NST2 ) + 1
            MCN( MNNBTQ-1 + NST3 ) = MCN( MNNBTQ-1 + NST3 ) + 1
C
C           numbers of the other vertices of the cube
            NST4 = NST1 + NBSOSU
            NST5 = NST2 + NBSOSU
            NST6 = NST3 + NBSOSU
            NST7 = 0
            NST8 = 0
C
            IF( NUEFTG .GT. 0 ) THEN
C              TRIANGLE A 6 TANGENTES => PENTAEDRE A 6 TANGENTES AUTRES 0
               IF( EPCOLA .GT. 0 ) THEN
C                 SENS DIRECT
                  NSTG1 = NOSOEL(4+1)
                  NSTG2 = NOSOEL(4+2)
                  NSTG3 = NOSOEL(4+3)
                  NSTG4 = NOSOEL(4+4)
                  NSTG5 = NOSOEL(4+5)
                  NSTG6 = NOSOEL(4+6)
               ELSE
C                 SENS INDIRECT
                  NSTG1 = NOSOEL(4+2)
                  NSTG2 = NOSOEL(4+1)
                  NSTG3 = NOSOEL(4+6)
                  NSTG4 = NOSOEL(4+5)
                  NSTG5 = NOSOEL(4+4)
                  NSTG6 = NOSOEL(4+3)
               ENDIF
               NSTG7 = 0
               NSTG8 = 0
           ENDIF
C
         ELSE
C
C           QUADRANGLE -> HEXAHEDRON
C
C           LES 4 NUMEROS DES SOMMETS DU QUADRANGLE I
C           the 4 numbers of vertices of the quadrangle I
            IF( EPCOLA .GT. 0 ) THEN
               NST1 = NOSOEL(1)
               NST2 = NOSOEL(2)
               NST3 = NOSOEL(3)
               NST4 = NOSOEL(4)
            ELSE
               NST1 = NOSOEL(1)
               NST2 = NOSOEL(4)
               NST3 = NOSOEL(3)
               NST4 = NOSOEL(2)
            ENDIF
C
C           accumulate the counting of each vertex
            MCN( MNNBTQ-1 + NST1 ) = MCN( MNNBTQ-1 + NST1 ) + 1
            MCN( MNNBTQ-1 + NST2 ) = MCN( MNNBTQ-1 + NST2 ) + 1
            MCN( MNNBTQ-1 + NST3 ) = MCN( MNNBTQ-1 + NST3 ) + 1
            MCN( MNNBTQ-1 + NST4 ) = MCN( MNNBTQ-1 + NST4 ) + 1
C
C           numbers of the other vertices of the cube
            NST5 = NST1 + NBSOSU
            NST6 = NST2 + NBSOSU
            NST7 = NST3 + NBSOSU
            NST8 = NST4 + NBSOSU
C
            IF( NUEFTG .GT. 0 ) THEN
C              QUADRANGLE A 8 TANGENTES => HEXAEDRE A 8 TANGENTES AUTRES 0
               IF( EPCOLA .GT. 0 ) THEN
C                 SENS DIRECT
                  NSTG1 = NOSOEL(4+1)
                  NSTG2 = NOSOEL(4+2)
                  NSTG3 = NOSOEL(4+3)
                  NSTG4 = NOSOEL(4+4)
                  NSTG5 = NOSOEL(4+5)
                  NSTG6 = NOSOEL(4+6)
                  NSTG7 = NOSOEL(4+7)
                  NSTG8 = NOSOEL(4+8)
               ELSE
C                 SENS INDIRECT
                  NSTG1 = NOSOEL(4+2)
                  NSTG2 = NOSOEL(4+1)
                  NSTG3 = NOSOEL(4+8)
                  NSTG4 = NOSOEL(4+7)
                  NSTG5 = NOSOEL(4+6)
                  NSTG6 = NOSOEL(4+5)
                  NSTG7 = NOSOEL(4+4)
                  NSTG8 = NOSOEL(4+3)
               ENDIF
            ENDIF
C
         ENDIF
C
         MCN( MNEFVC + 1 ) = NST1
         MCN( MNEFVC + 2 ) = NST2
         MCN( MNEFVC + 3 ) = NST3
         MCN( MNEFVC + 4 ) = NST4
         MCN( MNEFVC + 5 ) = NST5
         MCN( MNEFVC + 6 ) = NST6
         MCN( MNEFVC + 7 ) = NST7
         MCN( MNEFVC + 8 ) = NST8
         MNEFVC = MNEFVC + 8
C
C        LES 8 EVENTUELLES TANGENTES DE L'EF DE SURFACE SONT
C        POSITIONNEES PARMI LES 24 TANGENTES DE L'EF DE VOLUME
         IF( NBTGSU .GT. 0 ) THEN
            IF( NUEFTG .GT. 0 ) THEN
C              TRIANGLE OU QUADRANGLE A TG
C              LES 8 TANGENTES DU QUADRANGLE SONT POSITIONNES DANS LE CUBE
C              ORDRE DES TG DU QUADRANGLE:
C              T1:S1S2 T2:S1S4  T3:S2S3 T4:S2S1  T5:S3S4 T6:S3S2  T7:S4S1 T8:S4S
C              PENTAEDRE : NO SOMMET1 , NS2 , NS3, NS4, NS5, NS6, 0 , 0 ,
C                          NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                          NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S5),
C                          NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S6),
C                          NO TANGENT10(S4S5), NT11(S4S6), NT12(S4S1),
C                          NO TANGENT13(S5S6), NT14(S5S4), NT15(S5S2),
C                          NO TANGENT16(S6S4), NT17(S6S5), NT18(S6S3),
C                              0,               0,          0,
C                              0,               0,          0
C              HEXAEDRE  : NO SOMMET1 , NS2 , NS3, NS4, NS5, NS6, NS7 , NS8
C                          NO TANGENTE1(S1S2), NT2(S1S4),  NT3(S1S5),
C                          NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S6),
C                          NO TANGENTE7(S3S4), NT8(S3S2),  NT9(S3S7),
C                          NO TANGENT10(S4S1), NT11(S4S3), NT12(S4S8),
C                          NO TANGENT13(S5S6), NT14(S5S8), NT15(S5S1),
C                          NO TANGENT16(S6S7), NT17(S6S5), NT18(S6S2),
C                          NO TANGENT19(S7S8), NT20(S7S6), NT21(S7S3),
C                          NO TANGENT22(S8S5), NT23(S8S7), NT24(S8S4)
               MCN( MNEFVCTG + 1 ) = NSTG1
               MCN( MNEFVCTG + 2 ) = NSTG2
C
               MCN( MNEFVCTG + 4 ) = NSTG3
               MCN( MNEFVCTG + 5 ) = NSTG4
C
               MCN( MNEFVCTG + 7 ) = NSTG5
               MCN( MNEFVCTG + 8 ) = NSTG6
C
               MCN( MNEFVCTG + 10 ) = NSTG7
               MCN( MNEFVCTG + 11 ) = NSTG8
C              TOUTES LES AUTRES TANGENTES SONT NULLES ET C'EST DEJA FAIT
               MNEFVCTG = MNEFVCTG + 24
            ENDIF
         ENDIF
C
C        CALCUL DES 3 COMPOSANTES DU VECTEUR NORMAL AU TRIANGLE ou
C        QUADRANGLE I DE LA SURFACE
C        calculate the 3 components of the normal vector to triangle or
C        quadrangle I of the surface
C        DE SOMMETS NST1, NST2, et NST3
C        the vertices NST1, NST2 and NST3
         MNST1 = MNSTVC + 3 * NST1
         MNST2 = MNSTVC + 3 * NST2
         MNST3 = MNSTVC + 3 * NST3
         MNST4 = MNSTVC + 3 * NST4
         MNST5 = MNSTVC + 3 * NST5
         MNST6 = MNSTVC + 3 * NST6
         MNST7 = MNSTVC + 3 * NST7
         MNST8 = MNSTVC + 3 * NST8
C
C        LE VECTEUR NORMAL
C        the normal vector
         CALL VECNOR3( RMCN(MNST1+1), RMCN(MNST2+1),
     %                 RMCN(MNST3+1), VECN )
C
C        NORMALISATION EN LONGUEUR SELON L'EPAISSEUR
C        normalization in length by thickness
C        REMARQUE: EPCOLA<0 DONNE L'AUTRE SENS DU VECTEUR NORMAL!
C        remark:   EPCOLA<0 gives the other direction of the normal vector!
         RRR =abs(EPCOLA) / SQRT( VECN(1)**2 + VECN(2)**2 + VECN(3)**2 )
         VECN(1) = VECN(1) * RRR
         VECN(2) = VECN(2) * RRR
         VECN(3) = VECN(3) * RRR
C
C        record the current coordinates
         DO 56 J = 1, 4
            IF(EFTYPE .EQ. 0) THEN
               XYZTEMP(1,J) = RMCN(MNST(J+3) + 1)
               XYZTEMP(2,J) = RMCN(MNST(J+3) + 2)
               XYZTEMP(3,J) = RMCN(MNST(J+3) + 3)
            ELSE
               XYZTEMP(1,J) = RMCN(MNST(J+4) + 1)
               XYZTEMP(2,J) = RMCN(MNST(J+4) + 2)
               XYZTEMP(3,J) = RMCN(MNST(J+4) + 3)
            ENDIF
 56      CONTINUE
C
C        CALCUL DES COORDONNEES DES SOMMETS SUR LA NORMALE A UNE EPAISSEUR
C        calculate the coordinates of vertices along the normal vector
C        with the given thickness
         IF( EFTYPE .EQ. 0 ) THEN
C
C           PENTAHEDRON
C           ST1 -> ST4, ST2 -> ST5, ST3 -> ST6
            DO 57 J = 1, 3
               IF( RMCN( MNST(J+3) + 1 ) .EQ. VALTEM ) THEN
C
C                 SOMMET VU POUR LA PREMIERE FOIS. DEPLACEMENT NORMAL
C                 vertex visited for the first time.
C                 movement along the normal direction
C
C                 record the normal
                  RMCN(MNVECN-4 + 3*MNST(J) + 1) = REAL( VECN(1) )
                  RMCN(MNVECN-4 + 3*MNST(J) + 2) = REAL( VECN(2) )
                  RMCN(MNVECN-4 + 3*MNST(J) + 3) = REAL( VECN(3) )
C
C                 compute the coordinates
                  RMCN(MNST(J+3) + 1)=REAL(RMCN(MNST(J) + 1) + VECN(1))
                  RMCN(MNST(J+3) + 2)=REAL(RMCN(MNST(J) + 2) + VECN(2))
                  RMCN(MNST(J+3) + 3)=REAL(RMCN(MNST(J) + 3) + VECN(3))
C
               ELSE
C
C                 SOMMET VU POUR LA DEUXIEME FOIS. DEPLACEMENT NORMAL DE LA MOYE
C                 vertex visited not for the first time
                  RMCN(MNST(J+3) + 1) = REAL(RMCN(MNST(J+3) +1)+
     %                                       RMCN(MNST(J)   +1)+VECN(1))
                  RMCN(MNST(J+3) + 2) = REAL(RMCN(MNST(J+3) + 2) +
     %                                  RMCN(MNST(J)   + 2) + VECN(2))
                  RMCN(MNST(J+3) + 3) = REAL(RMCN(MNST(J+3) + 3) +
     %                                  RMCN(MNST(J)   + 3) + VECN(3))
               ENDIF
 57         CONTINUE
C
         ELSE
C
C           HEXAHEDRON
C           ST1 -> ST5, ST2 -> ST6, ST3 -> ST7, ST4 -> ST8
            DO 58 J = 1, 4
               IF( RMCN( MNST(J+4) + 1 ) .EQ. VALTEM ) THEN
C
C                 SOMMET VU POUR LA PREMIERE FOIS. DEPLACEMENT NORMAL
C                 vertex visited for the first time. move along the normal direc
C
C                 record the normal
                  RMCN(MNVECN-4 + 3*MNST(J) + 1) = REAL( VECN(1) )
                  RMCN(MNVECN-4 + 3*MNST(J) + 2) = REAL( VECN(2) )
                  RMCN(MNVECN-4 + 3*MNST(J) + 3) = REAL( VECN(3) )
C
C                 compute the coordinates
                  RMCN(MNST(J+4) + 1) = REAL(RMCN(MNST(J)+1) + VECN(1) )
                  RMCN(MNST(J+4) + 2) = REAL(RMCN(MNST(J)+2) + VECN(2) )
                  RMCN(MNST(J+4) + 3) = REAL(RMCN(MNST(J)+3) + VECN(3) )
C
               ELSE
C
C                 SOMMET VU POUR LA DEUXIEME FOIS. DEPLACEMENT NORMAL DE LA MOYE
C                 vertex visited not for the first time
                  RMCN(MNST(J+4) + 1) = REAL(RMCN(MNST(J+4) + 1) +
     %                                  RMCN(MNST(J)   + 1) + VECN(1) )
                  RMCN(MNST(J+4) + 2) = REAL(RMCN(MNST(J+4) + 2) +
     %                                  RMCN(MNST(J)   + 2) + VECN(2) )
                  RMCN(MNST(J+4) + 3) = REAL(RMCN(MNST(J+4) + 3) +
     %                                  RMCN(MNST(J)   + 3) + VECN(3) )
               ENDIF
 58         CONTINUE
C
         ENDIF
 60   CONTINUE
C
C     compute the "average" of the coordinates of the vertices
      MN = MNSOVC + WYZSOM + 3*NBSOSU - 1
      DO 70 I = 1, NBSOSU
         RMCN(MN + 1) = RMCN(MN + 1) / MCN( MNNBTQ-1 + I )
         RMCN(MN + 2) = RMCN(MN + 2) / MCN( MNNBTQ-1 + I )
         RMCN(MN + 3) = RMCN(MN + 3) / MCN( MNNBTQ-1 + I )
         MN = MN + 3
 70   CONTINUE
C
C     CUBES DES COUCHES LIMITES 2 a NBCOLA
      MNEFVC0 = MNFAVC + WUSOEF - 1
      DO 90 NC = 2, NBCOLA
         DO 80 I = 1, NBTQSU
C
            EFTYPE = MCN( MNEFVC0 + 8 )
C
            IF( EFTYPE .EQ. 0 ) THEN
C
C              PENTAHEDRON
C
C              LES NUMEROS DES SOMMETS 1, 2, et 3 DU NOUVEAU TRIANGLE
C              SONT LES SOMMETS 4, 5, et 6 DU PENTAHEDRON DE LA COUCHE PRECEDENT
C              the numbers 1, 2, and 3 of the vertices of the new triangle are
C              the numbers 4, 5, and 6 of the vertices of the previous pentahedr
               NST1 = MCN( MNEFVC0 + 4 )
               NST2 = MCN( MNEFVC0 + 5 )
               NST3 = MCN( MNEFVC0 + 6 )
               MNEFVC0 = MNEFVC0 + 8
C
C              numbers of the other vertices of the triangle I
               NST4 = NST1 + NBSOSU
               NST5 = NST2 + NBSOSU
               NST6 = NST3 + NBSOSU
               NST7 = 0
               NST8 = 0
C
            ELSE
C
C              PENTAHEDRON
C
C              LES NUMEROS DES SOMMETS 1, 2, 3, et 4 DU NOUVEAU QUADRANGLE
C              SONT LES SOMMETS 5, 6, 7, et 8 DU HEXAHEDRON DE LA COUCHE PRECEDE
C              the numbers 1, 2, 3, and 4 of the vertices of the new quadrangle
C              the numbers 5, 6, 7, and 8 of the vertices of the previous hexahe
               NST1 = MCN( MNEFVC0 + 5 )
               NST2 = MCN( MNEFVC0 + 6 )
               NST3 = MCN( MNEFVC0 + 7 )
               NST4 = MCN( MNEFVC0 + 8 )
               MNEFVC0 = MNEFVC0 + 8
C
C              numbers of the other vertices of the triangle I
               NST5 = NST1 + NBSOSU
               NST6 = NST2 + NBSOSU
               NST7 = NST3 + NBSOSU
               NST8 = NST4 + NBSOSU
C
            ENDIF
C
            MCN( MNEFVC + 1 ) = NST1
            MCN( MNEFVC + 2 ) = NST2
            MCN( MNEFVC + 3 ) = NST3
            MCN( MNEFVC + 4 ) = NST4
            MCN( MNEFVC + 5 ) = NST5
            MCN( MNEFVC + 6 ) = NST6
            MCN( MNEFVC + 7 ) = NST7
            MCN( MNEFVC + 8 ) = NST8
            MNEFVC = MNEFVC + 8
C
C           CALCUL DES 3 COMPOSANTES DU VECTEUR NORMAL A TRIANGLE I DE LA SURFAC
C           calculate the 3 components of the normal vector to triangle I
C           of the surface
            MNST1 = MNSTVC + 3 * NST1
            MNST2 = MNSTVC + 3 * NST2
            MNST3 = MNSTVC + 3 * NST3
            MNST4 = MNSTVC + 3 * NST4
            MNST5 = MNSTVC + 3 * NST5
            MNST6 = MNSTVC + 3 * NST6
            MNST7 = MNSTVC + 3 * NST7
            MNST8 = MNSTVC + 3 * NST8
C
C           LE VECTEUR NORMAL
C           the normal vector
            CALL VECNOR3( RMCN(MNST1+1), RMCN(MNST2+1),
     %                    RMCN(MNST3+1), VECN )
C           NORMALISATION EN LONGUEUR SELON L'EPAISSEUR
C           normalization in length by thickness
            RRR = abs(EPCOLA) / SQRT(VECN(1)**2+VECN(2)**2+VECN(3)**2)
            VECN(1) = VECN(1) * RRR
            VECN(2) = VECN(2) * RRR
            VECN(3) = VECN(3) * RRR
C
C           CALCUL DES COORDONNEES DES SOMMETS SUR LA NORMALE A UNE EPAISSEUR
C           calculate the coordinates of vertices along the normal vector
C           with the given thickness
C
            IF( EFTYPE .EQ. 0 ) THEN
C
C              PENTAHEDRON
C              ST1 -> ST4, ST2 -> ST5, ST3 -> ST6
C
               DO 71 J = 1, 3
                  IF( RMCN( MNST(J+3) + 1 ) .EQ. VALTEM ) THEN
C                    SOMMET VU POUR LA PREMIERE FOIS. DEPLACEMENT NORMAL
C                    vertex visited for the first time. move along the
C                    normal direction
                     RMCN(MNST(J+3) + 1) = REAL(RMCN(MNST(J)+1)+VECN(1))
                     RMCN(MNST(J+3) + 2) = REAL(RMCN(MNST(J)+2)+VECN(2))
                     RMCN(MNST(J+3) + 3) = REAL(RMCN(MNST(J)+3)+VECN(3))
                  ELSE
C                    SOMMET VU POUR LA DEUXIEME FOIS. DEPLACEMENT NORMAL DE
C                    LA MOYENNE
C                    vertex visited not for the first time
                     RMCN(MNST(J+3)+1) = REAL(RMCN(MNST(J+3)+1) +
     %                                   RMCN(MNST(J)  +1) + VECN(1))
                     RMCN(MNST(J+3)+2) = REAL(RMCN(MNST(J+3)+2) +
     %                                   RMCN(MNST(J)  +2) + VECN(2))
                     RMCN(MNST(J+3)+3) = REAL(RMCN(MNST(J+3)+3) +
     %                                   RMCN(MNST(J)  +3) + VECN(3))
                  ENDIF
 71            CONTINUE
C
            ELSE
C
C              HEXAHEDRON
C              ST1 -> ST5, ST2 -> ST6, ST3 -> ST7, ST4 -> ST8
C
               DO 72 J = 1, 4
                  IF( RMCN( MNST(J+4) + 1 ) .EQ. VALTEM ) THEN
C                    SOMMET VU POUR LA PREMIERE FOIS. DEPLACEMENT NORMAL
C                    vertex visited for the first time. move along the
C                    normal direction
                     RMCN(MNST(J+4) + 1) = REAL(RMCN(MNST(J)+1)+VECN(1))
                     RMCN(MNST(J+4) + 2) = REAL(RMCN(MNST(J)+2)+VECN(2))
                     RMCN(MNST(J+4) + 3) = REAL(RMCN(MNST(J)+3)+VECN(3))
                  ELSE
C                    SOMMET VU POUR LA DEUXIEME FOIS. DEPLACEMENT NORMAL DE
C                    LA MOYENNE
C                    vertex visited not for the first time
                     RMCN(MNST(J+4)+1) = REAL(RMCN(MNST(J+4)+1) +
     %                                   RMCN(MNST(J)  +1) + VECN(1))
                     RMCN(MNST(J+4)+2) = REAL(RMCN(MNST(J+4)+2) +
     %                                   RMCN(MNST(J)  +2) + VECN(2))
                     RMCN(MNST(J+4)+3) = REAL(RMCN(MNST(J+4)+3) +
     %                                   RMCN(MNST(J)  +3) + VECN(3))
                  ENDIF
 72            CONTINUE
            ENDIF
 80      CONTINUE
C
C        compute the average of the coordinates of the vertices
         MN = MNSOVC + WYZSOM + NC*3*NBSOSU - 1
         DO 85 I = 1, NBSOSU
            RMCN(MN + 1) = RMCN(MN + 1) / MCN( MNNBTQ-1 + I )
            RMCN(MN + 2) = RMCN(MN + 2) / MCN( MNNBTQ-1 + I )
            RMCN(MN + 3) = RMCN(MN + 3) / MCN( MNNBTQ-1 + I )
            MN = MN + 3
 85      CONTINUE
C
 90   CONTINUE
C
C     CUBES DES DERNIERES COUCHES a NBLACO
      DO 110 NC = 1, NBLACO
C        EPAISSEUR DE LA COUCHE  VARIATION LINEAIRE ENTRE EPCOLA ET EPLACO
         EPAIS = ( EPCOLA * (NBLACO-NC) + EPLACO * NC )/ NBLACO
         DO 100 I = 1, NBTQSU
C
            EFTYPE = MCN( MNEFVC0 + 8 )
C
            IF( EFTYPE .EQ. 0 ) THEN
C
C              PENTAHEDRON
C
C              LES NUMEROS DES SOMMETS 1, 2, et 3 DU NOUVEAU TRIANGLE
C              SONT LES SOMMETS 4, 5, et 6 DU PENTAHEDRON DE LA COUCHE PRECEDENT
C              the numbers 1, 2, and 3 of the vertices of the new triangle are
C              the numbers 4, 5, and 6 of the vertices of the previous pentahedr
               NST1 = MCN( MNEFVC0 + 4 )
               NST2 = MCN( MNEFVC0 + 5 )
               NST3 = MCN( MNEFVC0 + 6 )
               MNEFVC0 = MNEFVC0 + 8
C
C              numbers of the other vertices of the triangle I
               NST4 = NST1 + NBSOSU
               NST5 = NST2 + NBSOSU
               NST6 = NST3 + NBSOSU
               NST7 = 0
               NST8 = 0
C
            ELSE
C
C              PENTAHEDRON
C
C              LES NUMEROS DES SOMMETS 1, 2, 3, et 4 DU NOUVEAU QUADRANGLE
C              SONT LES SOMMETS 5, 6, 7, et 8 DU HEXAHEDRON DE LA COUCHE PRECEDE
C              the numbers 1, 2, 3, and 4 of the vertices of the new quadrangle
C              the numbers 5, 6, 7, and 8 of the vertices of the previous hexahe
               NST1 = MCN( MNEFVC0 + 5 )
               NST2 = MCN( MNEFVC0 + 6 )
               NST3 = MCN( MNEFVC0 + 7 )
               NST4 = MCN( MNEFVC0 + 8 )
               MNEFVC0 = MNEFVC0 + 8
C
C              numbers of the other vertices of the triangle I
               NST5 = NST1 + NBSOSU
               NST6 = NST2 + NBSOSU
               NST7 = NST3 + NBSOSU
               NST8 = NST4 + NBSOSU
C
            ENDIF
C
            MCN( MNEFVC + 1 ) = NST1
            MCN( MNEFVC + 2 ) = NST2
            MCN( MNEFVC + 3 ) = NST3
            MCN( MNEFVC + 4 ) = NST4
            MCN( MNEFVC + 5 ) = NST5
            MCN( MNEFVC + 6 ) = NST6
            MCN( MNEFVC + 7 ) = NST7
            MCN( MNEFVC + 8 ) = NST8
            MNEFVC = MNEFVC + 8
C
C           CALCUL DES 3 COMPOSANTES DU VECTEUR NORMAL A TRIANGLE I DE LA SURFAC
C           calculate the 3 components of the normal vector to triangle I of the
C
            MNST1 = MNSTVC + 3 * NST1
            MNST2 = MNSTVC + 3 * NST2
            MNST3 = MNSTVC + 3 * NST3
            MNST4 = MNSTVC + 3 * NST4
            MNST5 = MNSTVC + 3 * NST5
            MNST6 = MNSTVC + 3 * NST6
            MNST7 = MNSTVC + 3 * NST7
            MNST8 = MNSTVC + 3 * NST8
C
C           LE VECTEUR NORMAL
C           the normal vector
            CALL VECNOR3( RMCN(MNST1+1), RMCN(MNST2+1),
     %                    RMCN(MNST3+1), VECN )
C           NORMALISATION EN LONGUEUR SELON L'EPAISSEUR
C           normalization in length by thickness
            RRR = abs(EPAIS) / SQRT(VECN(1)**2+VECN(2)**2+VECN(3)**2)
            VECN(1) = VECN(1) * RRR
            VECN(2) = VECN(2) * RRR
            VECN(3) = VECN(3) * RRR
C
C           CALCUL DES COORDONNEES DES SOMMETS SUR LA NORMALE A UNE EPAISSEUR
C           calculate the coordinates of vertices along the normal vector
C           with the given thickness
C
            IF( EFTYPE .EQ. 0 ) THEN
C
C              PENTAHEDRON
C              ST1 -> ST4, ST2 -> ST5, ST3 -> ST6
C
               DO 91 J = 1, 3
                  IF( RMCN( MNST(J+3) + 1 ) .EQ. VALTEM ) THEN
C                    SOMMET VU POUR LA PREMIERE FOIS. DEPLACEMENT NORMAL
C                    vertex visited for the first time. move along the
C                    normal direction
                     RMCN(MNST(J+3) + 1) = REAL(RMCN(MNST(J)+1)+VECN(1))
                     RMCN(MNST(J+3) + 2) = REAL(RMCN(MNST(J)+2)+VECN(2))
                     RMCN(MNST(J+3) + 3) = REAL(RMCN(MNST(J)+3)+VECN(3))
                  ELSE
C                    SOMMET VU POUR LA DEUXIEME FOIS. DEPLACEMENT NORMAL
C                    DE
C                    LA MOYENNE
C                    vertex visited not for the first time
                     RMCN(MNST(J+3)+1) = REAL(RMCN(MNST(J+3)+1) +
     %                                   RMCN(MNST(J)  +1) + VECN(1))
                     RMCN(MNST(J+3)+2) = REAL(RMCN(MNST(J+3)+2) +
     %                                   RMCN(MNST(J)  +2) + VECN(2))
                     RMCN(MNST(J+3)+3) = REAL(RMCN(MNST(J+3)+3) +
     %                                   RMCN(MNST(J)  +3) + VECN(3))
                  ENDIF
 91            CONTINUE
C
            ELSE
C
C              HEXAHEDRON
C              ST1 -> ST5, ST2 -> ST6, ST3 -> ST7, ST4 -> ST8
C
               DO 92 J = 1, 4
                  IF( RMCN( MNST(J+4) + 1 ) .EQ. VALTEM ) THEN
C                    SOMMET VU POUR LA PREMIERE FOIS. DEPLACEMENT NORMAL
C                    vertex visited for the first time. move along the
C                    normal direction
                     RMCN(MNST(J+4) + 1) =REAL(RMCN(MNST(J)+1) +VECN(1))
                     RMCN(MNST(J+4) + 2) =REAL(RMCN(MNST(J)+2) +VECN(2))
                     RMCN(MNST(J+4) + 3) =REAL(RMCN(MNST(J)+3) +VECN(3))
                  ELSE
C                    SOMMET VU POUR LA DEUXIEME FOIS. DEPLACEMENT NORMAL
C                    DE
C                    LA MOYENNE
C                    vertex visited not for the first time
                     RMCN(MNST(J+4)+1) = REAL(RMCN(MNST(J+4)+1) +
     %                                   RMCN(MNST(J)  +1) + VECN(1))
                     RMCN(MNST(J+4)+2) = REAL(RMCN(MNST(J+4)+2) +
     %                                   RMCN(MNST(J)  +2) + VECN(2))
                     RMCN(MNST(J+4)+3) = REAL(RMCN(MNST(J+4)+3) +
     %                                   RMCN(MNST(J)  +3) + VECN(3))
                  ENDIF
 92            CONTINUE
            ENDIF
 100     CONTINUE
C
C        compute the average of the coordinates of the vertices
         MN = MNSOVC + WYZSOM + (NC+NBCOLA)*3*NBSOSU - 1
         DO 105 I = 1, NBSOSU
            RMCN(MN + 1) = RMCN(MN + 1) / MCN( MNNBTQ-1 + I )
            RMCN(MN + 2) = RMCN(MN + 2) / MCN( MNNBTQ-1 + I )
            RMCN(MN + 3) = RMCN(MN + 3) / MCN( MNNBTQ-1 + I )
            MN = MN + 3
 105     CONTINUE
C
 110  CONTINUE
C
C     DESTRUCTION DES TABLEAUX
      CALL TNMCDS('ENTIER', NBSOSU, MNNBTQ)
      CALL TNMCDS('REEL', 3*NBSOSU, MNVECN)
C
C     MISE A JOUR DU TABLEAU 'XYZSOMMET'
C     ----------------------------------
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOVC + WBCOOR ) = 3
C     LE NOMBRE DE SOMMETS
      MCN( MNSOVC + WNBSOM ) = NBSOVC
C     LE NOMBRE DE TANGENTES
      MCN( MNSOVC + WNBTGS ) = NBTGVC
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOVC) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOVC + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     MISE A JOUR DU TABLEAU 'NSEF'
C     -----------------------------
C     TYPE DE L'OBJET : VOLUME
      MCN( MNFAVC + WUTYOB ) = 4
      MCN( MNFAVC + WUTYMA ) = NUTYVC
C     NOMBRE DE SOMMETS PAR EF
      MCN( MNFAVC + WBSOEF ) = 8
C     NOMBRE D'EF DU MAILLAGE
      MCN( MNFAVC + WBEFOB ) = NBFAVC
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNFAVC + WUTFMA ) = -1
C     LES TANGENTES STOCKEES
      MCN( MNFAVC + WBTGEF ) = NBTGEF
      MCN( MNFAVC + WBEFAP ) = NBEFAP
      MCN( MNFAVC + WBEFTG ) = NBEFTG
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFAVC) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFAVC + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
      RETURN
      END
