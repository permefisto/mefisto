      SUBROUTINE XYZ1VXYZ( XYZPOI, NBEF,   NBNOEF, NUNOTE,
     %                     MOFACE, MXFACE, LFACES,
     %                     XYZ1, NTE1, XYZ,   NTEXYZ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EN PARTANT DU POINT XYZ1 CONTENU DANS LE TETRAEDRE NTE1
C -----    RETROUVER LE NUMERO NTEXYZ DU TETRAEDRE CONTENANT LE POINT XYZ
C          OU LE PLUS PROCHE MAIS PAS TROP LOIN (EN CAS DE TETRAEDRE DEGENERE)
C          SI LE PARCOURS XYZ1->XYZ SORT DU MAILLAGE ALORS NTEXYZ=0

C ENTREES:	
C --------
C XYZPOI : 3 COORDONNEES DES SOMMETS=POINTS DES TETRAEDRES
C NBEF   : NOMBRE DE TETRAEDRES DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN TETRAEDRE
C          (SEULS LES 4 SOMMETS SONT UTILES ICI et ONT POUR No 1 2 3 4)
C NUNOTE : NUMERO DES NBNOEF NOEUDS DES NBEF TETRAEDRES

C MOFACE : NOMBRE DE MOTS DE CHAQUE FACE DU TABLEAU LFACES
C MXFACE : NOMBRE MAXIMAL de FACES DU TABLEAU LFACES
C LFACES : LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C          SI SOMMET 2 < DERNIER SOMMET  => FACE   DIRECTE DANS L'ELEMENT FINI
C                      >                 => FACE INDIRECTE
C          UNE FACE DIRECTE EST VUE DE L EXTERIEUR DE L'ELEMENT FINI
C          SOUS LA FORME DIRECTE
C          LFACES(6,I)= NUMERO DU 1-ER ELEMENT FINI CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CET ELEMENT FINI
C                       <0 SI FACE INDIRECTE DANS CET ELEMENT FINI

C          SI LA FACE APPARTIENT A 2 ELEMENTS FINIS ALORS
C          LFACES(7,I)= NUMERO DU 2-EME ELEMENT FINI CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CET ELEMENT FINI
C                       <0 SI FACE INDIRECTE DANS CET ELEMENT FINI
C          SINON
C          LFACES(7,I)= NUMERO DE LA FACE FRONTALIERE SUIVANTE
C                       0 SI C'EST LA DERNIERE
C          L1FAFR = NUMERO DE LA PREMIERE FACE FRONTALIERE DANS LFACES
C                   CHAINAGE SUIVANT DANS LFACES(7,I)
C          L1FA2M = NUMERO DE LA PREMIERE FACE INTERFACE 2 MATERIAUX DANS LFACES
C                   CHAINAGE SUIVANT DANS LFACES(8,I)
C          SI VOLUME MULTI-MATERIAUX
C          LFACES(8,I)= NO DE LA FACE INTERFACE ENTRE 2 MATERIAUX SUIVANTE
C                       0 SI FACE SANS INFORMATION SUPPLEMENTAIRE
C          SINON VOLUME MONO-MATERIAU
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE
C                       0 SI FACE SANS TANGENTE

C NTE1   : NUMERO du TETRAEDRE CONTENANT LE POINT XYZ1
C XYZ1   : POINT DE DEPART DE LA RECHERCHE VERS NTEXYZ
C XYZ    : POINT DE TETRAEDRE NTEXYZ A RETROUVER

C SORTIES:
C --------
C NTEXYZ : >0  NUMERO DU TETRAEDRE CONTENANT XYZ INTERNE
C          <0 -NUMERO DU TETRAEDRE CONTENANT XYZ SUR UNE FACE FRONTALIERE
C              OU UN SOMMET FRONTALIER
C              XYZ EST DEPLACE A CE POINT SUR LA FRONTIERE
C          =0  PAS DE TETRAEDRE CONTENANT XYZ
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2020
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"

      COMMON / TRTETR /   STOPTE, TRACTE
      LOGICAL             STOPTE, TRACTE
ccc     %                               , TRACTE0
ccc      CHARACTER*80        KTITRE
      INTEGER             MXTEINT,      MXTEPAR
      PARAMETER         ( MXTEINT=1024, MXTEPAR=128 )

      REAL               XYZPOI(3,*)
ccc     %                  , XYZPT1(3), XYZPT2(3)

      INTEGER            NTEXYZ, NTE1, NTE, NTEMIN,
     %                   NBEF,   NBNOEF, NUNOTE(NBEF,NBNOEF),
     %                   MOFACE, MXFACE, LFACES(MOFACE,MXFACE),
     %                   NBTEINT, NOTEINT(MXTEINT), NOTEPAR(MXTEPAR)

      DOUBLE PRECISION   XYZ1(3), XYZ(3), XYZC(3),
     %                   Xe(3,4), DELTAe, CBTR(3), PTI(3), D,
     %                   CBMIN, CBMAX, CB1234,
     %                   CB1234MIN, CBTE1234(MXTEPAR), DELTATE(MXTEPAR)

      DOUBLE PRECISION   CBTE(4), CB1TE, CB2TE, CB3TE, CB4TE
      EQUIVALENCE       (CBTE(1),CB1TE), (CBTE(2),CB2TE),
     %                  (CBTE(3),CB3TE), (CBTE(4),CB4TE)

      INTEGER            I, K, NS, IMIN, IMAX, LINTER, NBFANE, NONOUI,
     %                   NTESUIV, NTEPREC, NTEPAR, NBTEPAR, NBOUCL,
     %                   NSINTE, NSAI1, NSAI2, NF1, NF2, NA, NSA1, NSA2

      INTEGER            NOSOFATE(3,4)
      DATA               NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      INTEGER            NOSOARTE(2,6)
      DATA               NOSOARTE / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /

      INTEGER            NOAR2FAC(4,4)
      DATA               NOAR2FAC / 0, 2, 3, 1,
     %                              2, 0, 6, 5,
     %                              3, 6, 0, 4,
     %                              1, 5, 4, 0 /
      INTRINSIC          DBLE, ABS

      TRACTE = .FALSE.

C     NOMBRE DE TETRAEDRES PARCOURUS POUR TROUVER XYZ
      NBTEPAR = 0
C     NOMBRE DE BOUCLES DETECTEES
      NBOUCL = 0

C     LE POINT COURANT XYZC POUR CHEMINER DE XYZ1 dans NTE1
C     VERS LE TETRAEDRE NTEXYZ CONTENANT XYZ
      XYZC( 1 ) = XYZ1( 1 )
      XYZC( 2 ) = XYZ1( 2 )
      XYZC( 3 ) = XYZ1( 3 )
      NTEPREC = 0
      NTEPAR  = 0
      NTESUIV = NTE1
      GOTO 15

C     TEST POUR DETECTER TROP DE TETRAEDRES PARCOURUS
 10   IF( NBTEPAR .GE. MXTEPAR ) THEN
         NTEXYZ = 0
         PRINT*,'xyz1vxyz: XYZ=',XYZ,' NTE1=',NTE1,
     %          ' TROP de TETRAEDRES PARCOURUS=',NBTEPAR,
     %          ' BOUCLE a DETECTER...  NTEXYZ=',NTEXYZ
         GOTO 9999

      ENDIF

      IF( NBTEPAR .GT. 16 ) THEN
C        pour voir la boucle ...
         PRINT*,'xyz1vxyz: NTE1=',NTE1,' Nb TETRA PARCOURUS=',NBTEPAR,
     %         'NTEPREC=',NTEPREC,' NTEPAR=',NTEPAR,' NTESUIV=',NTESUIV
         PRINT*,'xyz1vxyz: CB=',CBTE
      ENDIF

      IF( NTESUIV .EQ. NTE1 ) THEN

C        BOUCLE DETECTEE DANS LE PARCOURS
C        --------------------------------
         NBOUCL = NBOUCL + 1

C        LE PARCOURS BOUCLE DE NTE1 A NTESUIV EN PASSANT PAR NTEPREC et NTEPAR
ccc         PRINT*
ccc         PRINT*,'xyz1vxyz: BOUCLE',NBOUCL,' de',NBTEPAR,
ccc     %          ' TETRAEDRES de NTE1=NTESUIV',NTE1,
ccc     %          ' par NTEPREC=',NTEPREC,' NTEPAR=',NTEPAR
ccc         PRINT*,'xyz1vxyz: Tetra PARCOURUS:',(NOTEPAR(K),K=1,NBTEPAR)
ccc         PRINT*,'xyz1vxyz: CBTE1234:',(CBTE1234(K),K=1,NBTEPAR)
ccc         PRINT*,'xyz1vxyz: DELTAe  :',(DELTATE( K),K=1,NBTEPAR)

cccC        TRACE DES NBTEPAR TETRAEDRES DU PARCOURS et des POINTS XYZPTI XYZ
ccc         TRACTE0=TRACTE
cccccc         IF( NBTEPAR .GE. 4 ) TRACTE = .TRUE.
ccc      KTITRE='xyz1vxyz: NTE1=          BOUCLE de      TETRAEDRES pour AT
ccc     %TEINDRE +XYZ2'
ccc         WRITE(KTITRE(16:24),'(I9)') NTE1
ccc         WRITE(KTITRE(36:38),'(I3)') NBTEPAR
ccc         CALL SANSDBL( KTITRE, K )

ccc         XYZPT1(1) = REAL( XYZ1( 1 ) )
ccc         XYZPT1(2) = REAL( XYZ1( 2 ) )
ccc         XYZPT1(3) = REAL( XYZ1( 3 ) )

ccc         XYZPT2(1) = REAL( XYZ( 1 ) )
ccc         XYZPT2(2) = REAL( XYZ( 2 ) )
ccc         XYZPT2(3) = REAL( XYZ( 3 ) )
ccc         CALL TRANBTET( 0, KTITRE(1:K), XYZPOI, NBEF, NBNOEF, NUNOTE,
ccc     %                  NBTEPAR, NOTEPAR, XYZPT1, XYZPT2 )
ccc         TRACTE=TRACTE0

         IF( NBOUCL .GE. 2 ) THEN

C           SECONDE BOUCLE DETECTEE: IMPOSSIBLE ALLER PLUS LOIN
C           CHOIX DU TETRAEDRE NTEMIN DE PLUS FAIBLE CB1234
C           PARMI LES TETRAEDRES PARCOURUS
C           ---------------------------------------------------
            CB1234MIN = 1D100
            NTEMIN    = 0
            DO 11 I = 1, NBTEPAR
               NTE = NOTEPAR( I )
C              XYZ EST IL INTERNE AU TETRAEDRE NTE?
               CALL XYZD1TE( NBEF, NBNOEF, NUNOTE, XYZPOI,
     %                       XYZ,  NTE,    DELTAe, CBTE,   NONOUI )
               IF( NONOUI .EQ. -1 ) THEN
C                 RENCONTRE DU TETRAEDRE NTE DE VOLUME<0
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT *,'xyz1vxyz: RENCONTRE du TETRAEDRE',NTE,
     %                       ' de VOLUME*6=',Deltae
                  ELSE
                     PRINT *,'xyz1vxyz: PATH by TETRAHEDRON',NTE,
     %              ' of VOLUME*6=',Deltae,' <=0'
                  ENDIF
                  GOTO 11
               ENDIF

               CB1234 = ABS(CBTE(1)) +ABS(CBTE(2))
     %                 +ABS(CBTE(3)) +ABS(CBTE(4))
               IF( CB1234 .LT. CB1234MIN ) THEN
                  CB1234MIN = CB1234
                  NTEMIN    = NTE
                  IMIN      = I
               ENDIF

               IF( NONOUI .EQ. 1 ) THEN
C                 OUI: XYZ est INTERNE au TETRAEDRE NTE
                  NTEXYZ = NTE
                  GOTO 9999
               ENDIF
 11         ENDDO

C           PAS DE TETRAEDRE DE NOTEPAR CONTENANT XYZ
C           CHOIX DU TETRAEDRE NTEMIN DE CB1234 MINIMUM
            NTEXYZ = NTEMIN
ccc            PRINT*,'xyz1vxyz: FIN de BOUCLE CHOIX NTEMIN=NTEXYZ=',NTEXYZ
ccc     %            ,' CB1234=',CBTE1234(IMIN),' DELTAe=',DELTATE(IMIN)
            GOTO 9999

         ENDIF

         GOTO( 30, 40, 50, 55, 9900 ), NBFANE+1

      ENDIF


C     PARCOURS A PARTIR DU TETRAEDRE NTESUIV SUSCEPTIBLE DE CONTENIR XYZ
C     ==================================================================
C     UN TETRAEDRE PARCOURU DE PLUS
 15   NBTEPAR = NBTEPAR + 1

C     LE NO DU TETRAEDRE PARCOURU
      NOTEPAR( NBTEPAR ) = NTESUIV

      NTEPREC = NTEPAR
      NTEPAR  = NTESUIV

C     XYZ EST IL INTERNE AU TETRAEDRE NTEPAR?
C     ---------------------------------------
      CALL XYZD1TE( NBEF, NBNOEF, NUNOTE, XYZPOI, XYZ, NTEPAR,
     %              DELTAe, CBTE, NONOUI )
C      NONOUI=1  OUI XYZ   EST     INTERNE AU TETRAEDRE NTEPAR  CBTE CALCULE
C            =0  NON XYZ N'EST PAS INTERNE AU TETRAEDRE NTEPAR  CBTE CALCULE
C            =-1 TETRAEDRE NTE DE VOLUME<=0  CBTE NON CALCULE

      IF( NONOUI .EQ. -1 ) THEN
C        RENCONTRE DU TETRAEDRE NTEPAR DE VOLUME<0
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *,'xyz1vxyz: ATTENTION RENCONTRE du TETRAEDRE',NTEPAR,
     %              ' de VOLUME*6=',Deltae,' <=0'
         ELSE
            PRINT *,'xyz1vxyz: ATTENTION PATH by TETRAHEDRON',NTEPAR,
     %              ' of VOLUME*6=',Deltae,' <=0'
         ENDIF
         GOTO 9900
      ENDIF

      IF( NONOUI .EQ. 1 ) THEN

C        OUI: XYZ est INTERNE au TETRAEDRE NTEPAR=>NTEXYZ
C        ------------------------------------------------
         NTEXYZ = NTEPAR
         GOTO 9999

      ENDIF

C     LE TETRAEDRE NTEPAR NE CONTIENT PAS XYZ ALORS
C     RECHERCHE DE LA POSITION DE XYZ PAR RAPPORT A SES FACES, ARETES, SOMMETS
C     ------------------------------------------------------------------------
      CBTE1234( NBTEPAR ) = ABS(CBTE(1)) + ABS(CBTE(2))
     %                    + ABS(CBTE(3)) + ABS(CBTE(4))
      DELTATE(  NBTEPAR ) = DELTAe

C     NOMBRE DE FACES DE CB NEGATIVES DE XYZ DANS LE TETRAEDRE NTEPAR?
      NBFANE = 0
      IMIN   = 0
      CBMIN  = 1D100
      IMAX   = 0
      CBMAX  =-1D100
      NBFANE = 0
      DO I=1,4
         D = CBTE(I)
         IF( D .LE. 0D0 ) THEN
            NBFANE = NBFANE + 1
         ENDIF
         IF( D .LT. CBMIN ) THEN
            CBMIN = D
            IMIN = I
         ENDIF
         IF( D .GT. CBMAX ) THEN
            CBMAX = D
            IMAX = I
         ENDIF
      ENDDO

      GOTO( 30, 40, 60, 60, 9900 ), NBFANE+1

C     AUCUNE CB NEGATIVE: XYZ EST INTERNE A NTEPAR
C     --------------------------------------------
 30   NTEXYZ = NTEPAR
      GOTO 9999

C     UNE CB NEGATIVE CBMIN: XYZ EST DERRIERE ou SUR LA FACE IMIN+1 DE NTEPAR
C     -----------------------------------------------------------------------
C     RECHERCHE DU TETRAEDRE DERRIERE CETTE FACE IMIN+1 DE NTEPAR
 40   IF( IMIN .EQ. 4 ) THEN
         IMIN = 1
      ELSE
         IMIN = IMIN + 1
      ENDIF

C     RECUPERATION DES COORDONNEES Xe(3,4) DES 4 SOMMETS DU TETRAEDRE NTEPAR
      DO I=1,4
         NS = NUNOTE( NTEPAR, I )
         DO K=1,3
            Xe( K, I ) = DBLE( XYZPOI( K, NS ) )
         ENDDO
      ENDDO

C     XYZC LE BARYCENTRE DU TETRAEDRE NTEPAR
      DO K=1,3
         XYZC( K ) = ( Xe(K,1) + Xe(K,2) + Xe(K,3) + Xe(K,4) ) / 4
      ENDDO

C     CALCUL DU POINT D'INTERSECTION DE LA DROITE XYZC-XYZ et
C     DE LA FACE IMIN DE NTEPAR
      CALL INTARTR0( XYZC,  XYZ,
     %               Xe( 1, NOSOFATE(1,IMIN) ),
     %               Xe( 1, NOSOFATE(2,IMIN) ),
     %               Xe( 1, NOSOFATE(3,IMIN) ),
     %               LINTER, PTI, CBTR )
C LINTER: -3 SI XYZC=XYZ PAS DE CALCUL DE PTI
C         -2 SI XYZC-XYZ EST DANS  LE PLAN DU TRIANGLE
C         -1 SI XYZC-XYZ PARALLELE AU PLAN DU TRIANGLE SANS ETRE DEDANS
C         -1 SI XYZC-XYZ PARALLELE AU PLAN DU TRIANGLE  PAS DE CALCUL DE PTI
C          0 SI XYZC-XYZ N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
C          1 SI XYZC-XYZ   INTERSECTE     LE TRIANGLE ENTRE CES 3 SOMMETS
C            C-A-D  PTI EST INTERNE AU TRIANGLE OU SUR UNE DE SES 3 ARETES
C PTI   : LES 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C CBT   : LES 3 COORDONNEES BARYCENTRIQUES DE PTI DANS LA FACE

      IF( LINTER .EQ. -3 ) THEN
C        XYZC = XYZ  est CONTENU DANS NTEPAR
         NTEXYZ = NTEPAR
         GOTO 9999
      ENDIF

      IF( LINTER .EQ. -1 ) THEN

C        XYZC-XYZ PARALLELE AU PLAN 
C        PAS DE PTI POINT D'INTERSECTION CALCULE
C        LA FACE EST AU NIVEAU DU BARYCENTRE => NTEPAR EST ECRASE

ccc         PRINT *
ccc         PRINT *,'xyz1vxyz:  XYZC-XYZ PARALLELE a la FACE',
ccc     %            IMIN,' de NTEPAR=',NBTEPAR
ccc         PRINT *,'xyz1vxyz: XYZC=',XYZC
ccc         PRINT *,'xyz1vxyz: XYZ =',XYZ
ccc         PRINT *,'xyz1vxyz: XS1F=',(Xe(K,NOSOFATE(1,IMIN)),K=1,3)
ccc         PRINT *,'xyz1vxyz: XS2F=',(Xe(K,NOSOFATE(2,IMIN)),K=1,3)
ccc         PRINT *,'xyz1vxyz: XS3F=',(Xe(K,NOSOFATE(3,IMIN)),K=1,3)
ccc         PRINT *,'xyz1vxyz: PTI =',PTI
cccccc         PRINT *,'xyz1vxyz: CBTR=',CBTR     NON CALCULE

cccC        TRACE DES NBTEPAR TETRAEDRES DU PARCOURS et du POINT XYZ
ccc         TRACTE0 = TRACTE
cccccc         TRACTE = .TRUE.
ccc         KTITRE ='xyz1vxyz: ETRANGE XYZC-XYZ N''INTERSECTE PAS LA FACE'
ccc         CALL SANSDBL( KTITRE, K )
ccc         XYZPT1(1) = REAL( XYZC( 1 ) )
ccc         XYZPT1(2) = REAL( XYZC( 2 ) )
ccc         XYZPT1(3) = REAL( XYZC( 3 ) )

ccc         XYZPT2(1) = REAL( XYZ( 1 ) )
ccc         XYZPT2(2) = REAL( XYZ( 2 ) )
ccc         XYZPT2(3) = REAL( XYZ( 3 ) )
ccc         CALL TRANBTET( 0, KTITRE(1:K), XYZPOI, NBEF, NBNOEF, NUNOTE,
ccc     %                  NBTEPAR, NOTEPAR, XYZPT1, XYZPT2 )
ccc         TRACTE = TRACTE0

         GOTO 9900

      ELSE IF( LINTER .EQ. -2 ) THEN

C        XYZC-XYZ DANS LE PLAN   A AMELIORER...
         print*,'xyz1vxyz: XYZC-XYZ DANS LE PLAN de la FACE IMIN. A TRAI
     %TER...'
         GOTO 9900

      ELSE

C        XYZC DEVIENT LE POINT D'INTERSECTION PTI de la DROITE XYZC-XYZ
C        et DE LA FACE IMIN DE NTEPAR
         XYZC(1) = PTI(1)
         XYZC(2) = PTI(2)
         XYZC(3) = PTI(3)

      ENDIF

C     RECHERCHE DU TETRAEDRE NTESUIV DERRIERE LA FACE IMIN DE NTEPAR
      CALL NTED1FAC( NBEF,   NBNOEF, NUNOTE,
     %               MOFACE, MXFACE, LFACES,
     %               IMIN,   NTEPAR, NTESUIV )

      IF( NTESUIV .EQ. 0 ) THEN

C        FACE FRONTIERE: PAS DE TETRAEDRE DERRIERE LA FACE IMIN DE NTEPAR
C        IMPOSSIBLE D'ALLER AU DELA
         NTEXYZ = -NTEPAR
         XYZ( 1 ) = PTI( 1 )
         XYZ( 2 ) = PTI( 2 )
         XYZ( 3 ) = PTI( 3 )
         GOTO 9999

      ELSE

C        TETRAEDRE OPPOSE A LA FACE IMIN DE NTEPAR
         NTEXYZ = NTESUIV
         GOTO 10

      ENDIF


C     2 CB NEGATIVES: XYZ EST DERRIERE L'ARETE COMMUNE AUX 2 FACES DE CB<0
C     --------------------------------------------------------------------
C     TROUVER A PARTIR DE LA FACE DU TETRAEDRE NTEPAR QUI VOIT
C     LE MIEUX LE POINT XYZ A ATTEINDRE, LE TETRAEDRE NTESUIV
C     QUI EST DERRIERE CETTE FACE DE NTEPAR

C     RECHERCHE DU NO DES 2 SOMMETS NSAI1-NSAI2 DE CETTE ARETE de NTEPAR
C     APPARTENANT A 2 FACES DE CB NEGATIVES
 50   NF1 = 0
      NF2 = 0
      DO I=1,4
         D = CBTE(I)
         IF( D .LE. 0D0 ) THEN
            IF( NF1 .EQ. 0 ) THEN
C              NO DU 1-ER CB<0
               NF1 = I
            ELSE
C              NO DU 2-EME CB<0
               NF2 = I
               GOTO 51
            ENDIF
         ENDIF
      ENDDO

C     QUELLE EST L'ARETE K de NTEPAR COMMUNE AUX 2 FACES NF1+1 NF2+1?
 51   IF( NF1 .EQ. 4 ) THEN
         NF1 = 1
      ELSE
         NF1 = NF1 + 1
      ENDIF
      IF( NF2 .EQ. 4 ) THEN
         NF2 = 1
      ELSE
         NF2 = NF2 + 1
      ENDIF
      NA = NOAR2FAC( NF1, NF2 )

C     ARETE NA DE SOMMETS NSAI1 NSAI2 DANS LE TETRAEDRE NTEPAR
      NSAI1 = NUNOTE( NTEPAR, NOSOARTE(1,NA) )
      NSAI2 = NUNOTE( NTEPAR, NOSOARTE(2,NA) )

C     RECHERCHE DES TETRAEDRES NTEPAR D'ARETE NSAI1 NSAI2
      CALL TETR1AR(NBEF, NUNOTE, NSAI1,NSAI2, MXTEINT, NBTEINT, NOTEINT)

C     RETRAIT DU TETRAEDRE NTEPAR et NTE1 DE CETTE LISTE
      K = 0
      DO I=1,NBTEINT
         NTE = NOTEINT( I )
         IF( NTE .NE. NTEPAR .AND. NTE .NE. NTE1 ) THEN
            K = K + 1
            NOTEINT( K ) = NTE
         ENDIF
      ENDDO
      NBTEINT = K

      IF( NBTEINT .EQ. 0 ) THEN
C        PAS DE TETRAEDRE D'ARETE NSAI1-NSAI2 POUR CONTINUER
C        2 TETRAEDRES POUR CETTE ARETE. ELLE DOIT ETRE FRONTIERE
C        CALCUL DE PTI POINT PROJETE DE XYZ SUR L'ARETE
         DO I=1,4
            NS = NUNOTE( NTEPAR, I )
            DO K=1,3
               Xe( K, I ) = DBLE( XYZPOI( K, NS ) )
            ENDDO
         ENDDO

         NSA1 = NOSOARTE(1,NA)
         NSA2 = NOSOARTE(2,NA)
         CALL PTPRDRD( XYZ, Xe(1,NSA1), Xe(1,NSA2),  PTI )

C        PTI EST IL ENTRE LES 2 SOMMETS DE L'ARETE
         D = ( ( PTI(1)-Xe(1,NSA1) )**2 + ( PTI(2)-Xe(2,NSA1) )**2
     %       + ( PTI(3)-Xe(3,NSA1) )**2 ) /
     %       ( ( Xe(1,NSA2)-Xe(1,NSA1) )**2 + (Xe(2,NSA2)-Xe(2,NSA1))**2
     %       + ( Xe(3,NSA2)-Xe(3,NSA1) )**2 )

         IF( D .LE. 1D0 ) THEN

C           RECHERCHE DU SIGNE NEGATIF SI PTI EST EXTERIEUR A NSAI1-NSAI2
            IF( ( PTI(1)-Xe(1,NSA1) )*(Xe(1,NSA2)-Xe(1,NSA1) ) .LT. 0D0
     %     .OR. ( PTI(2)-Xe(2,NSA1) )*(Xe(2,NSA2)-Xe(2,NSA1) ) .LT. 0D0
     %     .OR. ( PTI(3)-Xe(3,NSA1) )*(Xe(3,NSA2)-Xe(3,NSA1) ) .LT. 0D0)
     %      THEN
C             PTI EST EXTERIEUR A L'ARETE NSAI1-NSAI2 FRONTIERE
C             XYZ = MILIEU DE L'ARETE NSAI1-NSAI2 LE NOUVEAU POINT XYZ
              XYZ(1)=DBLE( XYZPOI( 1, NSAI1 ) + XYZPOI( 1, NSAI2 ) ) / 2
              XYZ(2)=DBLE( XYZPOI( 2, NSAI1 ) + XYZPOI( 2, NSAI2 ) ) / 2
              XYZ(3)=DBLE( XYZPOI( 3, NSAI1 ) + XYZPOI( 3, NSAI2 ) ) / 2
              GOTO 53
           ENDIF

         ENDIF

C        PTI INTERIEUR A L'ARETE NSAI1-NSAI2 => XYZ=PTI SUR L'ARETE FRONTALIERE
         XYZ( 1 ) = PTI( 1 )
         XYZ( 2 ) = PTI( 2 )
         XYZ( 3 ) = PTI( 3 )

ccc         PRINT*,'xyz1vxyz: NTE1=',NTE1,' NTEPAR=',NTEPAR,
ccc     %          ' SEULS TETRAEDRES D''ARETE',NSAI1,NSAI2,
ccc     %          ' SORTIE avec NTEXYZ=',NTEXYZ,' et NPI=XYZ=',XYZ

 53      NTEXYZ = -NTEPAR
         GOTO 9999

      ENDIF

C     RECHERCHE DU TETRAEDRE CONTENANT XYZ DANS CES NBTEINT TETRAEDRES
      GOTO 70


C     3 CB NEGATIVES: XYZ EST DERRIERE LE SOMMET IMAX DU TETRAEDRE NTEPAR
C                     ET DANS LE TETRAEDRE OPPOSE AUX 3 FACES DE CB<0
C     -------------------------------------------------------------------
 55   NSINTE = NUNOTE( NTEPAR, IMAX )

C     XYZC = XYZPOI(NSINTE) LE NOUVEAU POINT DE DEPART
      XYZC( 1 ) = DBLE( XYZPOI( 1, NSINTE ) )
      XYZC( 2 ) = DBLE( XYZPOI( 2, NSINTE ) )
      XYZC( 3 ) = DBLE( XYZPOI( 3, NSINTE ) )

C     RECHERCHE DES TETRAEDRES DE SOMMET NSINTE.  TRES COUTEUX!...
      CALL TETR1ST( NBEF, NUNOTE, NSINTE, MXTEINT, NBTEINT, NOTEINT )

C     RETRAIT DU TETRAEDRE NTEPAR et NTE1 DE CETTE LISTE
      K = 0
      DO I=1,NBTEINT
         NTE = NOTEINT( I )
         IF( NTE .NE. NTEPAR .AND. NTE .NE. NTE1 ) THEN
            K = K+1
            NOTEINT( K ) = NTE
         ENDIF
      ENDDO
      NBTEINT = K

 59   IF( NBTEINT .EQ. 0 ) THEN
C        PAS DE TETRAEDRE DE SOMMET NSINTE POUR CONTINUER
         NTEXYZ = -NTEPAR
         XYZ( 1 ) = XYZC( 1 )
         XYZ( 2 ) = XYZC( 2 )
         XYZ( 3 ) = XYZC( 3 )
ccc         PRINT*,'xyz1vxyz: NTE1=',NTE1,' NTEPAR=',NTEPAR,
ccc     %          ' SEULS TETRAEDRES DE SOMMET NSINTE=',NSINTE,
ccc     %          ' SORTIE avec NTEXYZ=',NTEXYZ,' XYZ SOMMET=',XYZ
         GOTO 9999
      ENDIF

C     RECHERCHE DU TETRAEDRE CONTENANT XYZ DANS CES NBTEINT TETRAEDRES
      GOTO 70


C     TROUVER A PARTIR DE LA FACE DU TETRAEDRE NTEPAR QUI VOIT
C     LE MIEUX LE POINT XYZ A ATTEINDRE, LE TETRAEDRE NTESUIV
C     QUI EST DERRIERE CETTE FACE DE NTEPAR
C     XYZC DEVIENT LE BARYCENTRE DE CETTE FACE DE NTEPAR et NTESUIV
C     -------------------------------------------------------------
 60   CALL NFAQVXYZ( NBEF,   NBNOEF, NUNOTE,
     %               MOFACE, MXFACE, LFACES,
     %               XYZPOI, XYZ,    NTEPAR,   XYZC, NTESUIV )

      IF( NTESUIV .EQ. 0 ) THEN

C        PAS DE FACE DU TETRAEDRE NTEPAR QUI VOIE XYZ
         IF( NBFANE .EQ. 2 ) THEN
            GOTO 50
         ELSE IF( NBFANE .EQ. 3 ) THEN
            GOTO 55
         ELSE
            NBTEINT = 0
            GOTO 59
         ENDIF

      ELSE

C        RECHERCHE A PARTIR DU TETRAEDRE NTESUIV  
         NTEXYZ = NTESUIV
         GOTO 10

      ENDIF


C     RECHERCHE DU TETRAEDRE CONTENANT XYZ DANS CES NBTEINT TETRAEDRES
C     ----------------------------------------------------------------
 70   CB1234MIN = 1D100
      NTEMIN    = 0
      DO I = 1, NBTEINT
         NTE = NOTEINT( I )
C        XYZ EST IL INTERNE AU TETRAEDRE NTE?
         CALL XYZD1TE( NBEF, NBNOEF, NUNOTE, XYZPOI,
     %                 XYZ,  NTE,    DELTAe, CBTE,   NONOUI )

         CB1234 = ABS(CBTE(1)) +ABS(CBTE(2)) +ABS(CBTE(3)) +ABS(CBTE(4))
         IF( CB1234 .LT. CB1234MIN ) THEN
            CB1234MIN = CB1234
            NTEMIN    = NTE
         ENDIF

         IF( NONOUI .EQ. -1 ) THEN
C           RENCONTRE DU TETRAEDRE NTE DE VOLUME<0
            IF( LANGAG .EQ. 0 ) THEN
            PRINT *,'xyz1vxyz: RENCONTRE du TETRAEDRE',NTE,
     %              ' de VOLUME*6=',Deltae
            ELSE
            PRINT *,'xyz1vxyz: PATH by TETRAHEDRON',NTE,
     %              ' of VOLUME*6=',Deltae,' <=0'
            ENDIF
            GOTO 9900
         ENDIF

         IF( NONOUI .EQ. 1 ) THEN
C           OUI: XYZ est INTERNE au TETRAEDRE NTE
            NTEXYZ = NTE
            GOTO 9999
         ENDIF
      ENDDO

C     PAS DE TETRAEDRE DE NOTEINT CONTENANT XYZ
C     ESSAI A PARTIR DU TETRAEDRE NTEMIN DE CB1234 MINIMUM
      NTESUIV = NTEMIN
      GOTO 10


C     ERREUR: TETRAEDRE NTEPAR AVEC DELTAe<=0 RENCONTRE
 9900 NTEXYZ = 0
      print*,'xyz1vxyz: DELTAe=VOLUME<0 de TETRAEDRE XYZ=',XYZ,
     %       ' NTE1=',NTE1,' sortie avec NTEXYZ=', NTEXYZ

 9999 RETURN
      END
