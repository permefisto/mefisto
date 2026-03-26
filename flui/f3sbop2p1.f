      SUBROUTINE F3SBOP2P1( NEF00,   NBEF,   NUNOTE,
     %                      MOFACE,  MXFACE, LFACES,
     %                      DT,      CoefTr, XYZNOE, NDDLNO,
     %                      TempouVit, TEMPER, VXYZPNtn, VITMXtn,
     %                      BE,      IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU VECTEUR TRANSPORT NAVIER STOKES DU TETRAEDRE TAYLOR HOOD
C ----- INTEGRATION RETROGRADE DES CARACTERISTIQUES cf O. PIRONNEAU
C       POLYNOME LAGRANGE DE DEGRE 2 POUR LA VITESSE
C       POLYNOME LAGRANGE DE DEGRE 1 POUR LA PRESSION
C       TETRAEDRE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C       INTEGRATION EN 15 POINTS DU TETRAEDRE P5-EXACTE de
C       Si TempouVit=0 calcul de
C          CoefTr Integrale tP2 TEMPER(tn,X(tn;tn+1,x) dx  (TEMPERATURE)
C       Sinon Si TempouVit=1 calcul de
C          CoefTr Integrale tP2 ui(tn,X(tn;tn+1,x) dx   pour i=1,2
C       Sinon Si TempouVit=2 calcul de
C          CoefTr Integrale tP2 ( TEMPER(tn,x)
C                                -u1(tn,x) d TEMPER(tn,x)/dx1
C                                -u2(tn,x) d TEMPER(tn,x)/dx2
C                                -u3(tn,x) d TEMPER(tn,x)/dx3 ) dx

C ENTREES:
C --------
C NEF00  : NUMERO DANS LE MAILLAGE DU TETRAEDRE TAYLOR-HOOD A TRAITER
C NBEF   : NOMBRE DE TETRAEDRES DU MAILLAGE
C NUNOTE : NUMERO DES 4 SOMMETS ET 6 MILIEUX DES ARETES DES TETRAEDRES

C MOFACE : NOMBRE DE MOTS DE CHAQUE FACE DU TABLEAU LFACES
C MXFACE : NOMBRE MAXIMAL D'FACES DU TABLEAU LFACES
C LFACES : TABLEAU DES FACES DU MAILLAGE  (cf hachag.f)
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C                       1 CUBE EST UN TETRA ou PENTA ou HEXAEDRE
C          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                       0 SI PAS DE 1-ER  CUBE
C          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                       ou CHAINAGE SUR LA FACE FRONTALIERE SUIVANTE
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
C          LFACES(9,I)= NUMERO DU TYPE DU DERNIER EF DE CETTE FACE
C                       CETTE INFORMATION EXISTE SEULEMENT POUR UN OBJET
C                       PAS POUR UNE SURFACE OU UN VOLUME!

C DT     : LE PAS CONSTANT DE TEMPS ENTRE tn ET tn+1
C CoefTr : COEFFICIENT DE LA VITESSE TRANSPORTEE
C XYZNOE : 3 COORDONNEES DES SOMMETS ET MILIEUX DES ARETES DES TETRAEDRES
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE A 1+NBNOEU

C TempouVit: =0 calcul de CoefTr Integrale tP2 TEMPER(tn,X(tn;tn+1,x) dx
C            =1 CoefTr Integrale tP2 Som ui(tn,x) d TEMPER(tn)/dxi dx

C TEMPER : VECTEUR GLOBAL des DL TEMPERATURE au TEMPS tn
C VXYZPNtn:VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDL)
C          NOEUD PAR NOEUD  au TEMPS tn
C VITMXtn: NORME DE LA VITESSE MAXIMALE AU TEMPS PRECEDANT

C SORTIES:
C --------
C BE     : VECTEUR ELEMENTAIRE DE 3x10 COEFFICIENTS + 4 PRESSIONS=0 ICI
C          LES 4 DERNIERS COEFFICIENTS BE(31:34) SONT MIS A ZERO
C IERR   : CODE D'ERREUR 0 PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Mars 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Mai 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C MODIFS : ALAIN PERRONNET St Pierre du Perray              Octobre 2020
C MODIFS : ALAIN PERRONNET Veulettes sur mer              Septembre 2022
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2023
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      INTEGER            LECTEU, IMPRIM, INTERA, NUNITE
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)

      REAL               DT, XYZNOE(3,*)
      INTEGER            MOFACE, MXFACE, LFACES(MOFACE,MXFACE),
     %                   NEF00, NBEF, NUNOTE(NBEF,10)
      INTEGER            NDDLNO(0:*), IERR, TempouVit, MOTBE
      INTEGER            NEF, NEF0, NS, I, K, L, ITERT, LINTER,
     %                   IK,  NDL, NBSSPATE, NOFAMX

      DOUBLE PRECISION   TEMPER(1:*), VXYZPNtn(1:*), VITMXtn
      DOUBLE PRECISION   BE(34), DETM33, TEMPERA
      DOUBLE PRECISION   DDT, D, CoefTr, RDPOIDS, PRCBL1(10), PRCBL0(10)
ccc   DOUBLE PRECISION   VNORM
      DOUBLE PRECISION   VXNUEF, CBTR(3), PTITR(3),
     %                   V(3), VA(3), Xe(3,4), DELTAe, DELTAe1,
     %                   X1, Y1, Z1, CB1BL, CB2BL, CB3BL, CB4BL,
     %                   XEBL0, YEBL0, ZEBL0, XEBL1, YEBL1, ZEBL1,
     %                   XYZ0(3), XYZ1(3)
      EQUIVALENCE      (XYZ0(1),XEBL0), (XYZ0(2),YEBL0), (XYZ0(3),ZEBL0)
      EQUIVALENCE      (XYZ1(1),XEBL1), (XYZ1(2),YEBL1), (XYZ1(3),ZEBL1)

C     LE NUMERO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER            NOSOFATE(3,4)
      DATA               NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      INTRINSIC          NINT, SQRT, DBLE, ABS

C     POIDS ET POINTS D'INTEGRATION SUR LE TETRAEDRE UNITE P5 EXACTE
C     DOUBLE PRECISION   POTEP5( NBPTI=15 ), COTEP5( 3, NBPTI )
      include"./incl/tep5poco.inc"

      IERR = 0
C     INITIALISATION A ZERO DU VECTEUR ELEMENTAIRE BE
      IF( TempouVit .EQ. 0 ) THEN
         MOTBE = 10
      ELSE
         MOTBE = 30
      ENDIF
      DO I=1,MOTBE
         BE(I) = 0D0
      ENDDO

C     CALCUL DU NOMBRE NBSSPATE DE SOUS PAS DE TEMPS DDT DANS LE PAS DE TEMPS DT
C     --------------------------------------------------------------------------
      CALL SPTP2P13D( DT,  XYZNOE, NEF00, NBEF, NUNOTE, NDDLNO,
     %                VITMXtn, VXYZPNtn,  VXNUEF, NBSSPATE, DDT )

C     INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C     RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS DU TETRAEDRE NEF
      DO I=1,4
         NS = NUNOTE(NEF00,I)
         DO K=1,3
            Xe(K,I) = DBLE( XYZNOE( K, NS ) )
         ENDDO
      ENDDO

C     CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF INITIAL
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      Z1 = Xe(3,1)
      DELTAe1 = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                  Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                  Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )

      IF( DELTAe1 .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'f3sbop2p1: ATTENTION EF',NEF00,
     %                 ' de VOLUME*6=',DELTAe1,' Be NON PRIS EN COMPTE'
         ELSE
            WRITE(IMPRIM,*) 'f3sbop2p1: ATTENTION FE',NEF00,
     %                    ' of VOLUME*6=',DELTAe1,' Be is NOT COMPUTED'
         ENDIF
         GOTO 9999
      ENDIF

C     CALCUL DU VECTEUR tPchapeau(bl) ( Pchapeau(X(tn;tn+1,Fe(bl)) Vitesse(k) )
C     =========================================================================
C     Boucle sur les points d'integration dans le tetraedre NEF00
      DO 100 L=1,NBPTI

C        INITIALISATION DE LA VITESSE ET DE L'EF INITIAL DE RECHERCHE
         IERR = 0
         NEF0 = 0
         NEF  = NEF00

C        RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS DU TETRAEDRE NEF
         DO I=1,4
            NS = NUNOTE(NEF,I)
            DO K=1,3
               Xe(K,I) = DBLE( XYZNOE( K, NS ) )
            ENDDO
         ENDDO

C        CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF INITIAL
         X1 = Xe(1,1)
         Y1 = Xe(2,1)
         Z1 = Xe(3,1)
         DELTAe = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                    Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                    Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )
         IF( DELTAe .LE. 0D0 ) THEN
C           JACOBIEN de l'EF NEF NUL ou NEGATIF
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'f3sbop2p1: attention EF',NEF,
     %         ' de VOLUME*6=',DELTAe,' Be NON CALCULE pour L=',L
            ELSE
               WRITE(IMPRIM,*) 'f3sbop2p1: attention FE',NEF,
     %         ' of VOLUME*6=',DELTAe,' Be is NOT COMPUTED for L=',L
            ENDIF
            GOTO 100
         ENDIF

C        COORDONNEES DU POINT bl D'INTEGRATION DANS EF REFERENCE
C        LES 3 DERNIERES COORDONNEES BARYCENTRIQUES SONT AUSSI X Y Z
C        DANS L'ELEMENT DE REFERENCE
         CB2BL = COTEP5( 1, L )
         CB3BL = COTEP5( 2, L )
         CB4BL = COTEP5( 3, L )
         CB1BL = 1.D0 - CB2BL - CB3BL - CB4BL

C        COORDONNEES (XEBL1,YEBL1,ZEBL1)=Fe(bl) SUR LE TETRAEDRE NEF a tn+1
C       (XEBL1,YEBL1,ZEBL1)=Fe(CB2BL,CB3BL,CB4BL) SUR LE TETRAEDRE NEF
         XEBL1=CB1BL*Xe(1,1)+CB2BL*Xe(1,2)+CB3BL*Xe(1,3)+CB4BL*Xe(1,4)
         YEBL1=CB1BL*Xe(2,1)+CB2BL*Xe(2,2)+CB3BL*Xe(2,3)+CB4BL*Xe(2,4)
         ZEBL1=CB1BL*Xe(3,1)+CB2BL*Xe(3,2)+CB3BL*Xe(3,3)+CB4BL*Xe(3,4)

C        CALCUL DES Pchapeau(bl)  POLYNOMES LAGRANGE DE DEGRE 2 AU PT L
         PRCBL1(1) = CB1BL * ( 2.D0 * CB1BL - 1.D0 )
         PRCBL1(2) = 2.D0 * CB2BL * CB2BL - CB2BL
         PRCBL1(3) = 2.D0 * CB3BL * CB3BL - CB3BL
         PRCBL1(4) = 2.D0 * CB4BL * CB4BL - CB4BL
         D         = 4.D0 * CB1BL
         PRCBL1(5) = D    * CB2BL
         PRCBL1(6) = 4.D0 * CB2BL * CB3BL
         PRCBL1(7) = D    * CB3BL
         PRCBL1(8) = D    * CB4BL
         PRCBL1(9) = 4.D0 * CB2BL * CB4BL
         PRCBL1(10)= 4.D0 * CB3BL * CB4BL

C        CALCUL DE LA VITESSE V(tn,(CB2BL,CB3BL,CB4BL)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         V(3) = 0D0
         DO I=1,10
C           NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUNOTE(NEF,I) - 1 )
            V(1) = V(1) +  PRCBL1(I) * VXYZPNtn( NDL+1 )
            V(2) = V(2) +  PRCBL1(I) * VXYZPNtn( NDL+2 )
            V(3) = V(3) +  PRCBL1(I) * VXYZPNtn( NDL+3 )
         ENDDO

         IF( NBSSPATE .LE. 0 ) THEN
C           SI VITESSE QUASI NULLE PAS DE REMONTEE DE LA CARACTERISTIQUE
C           Les  Pchapeau(XYRCBL) des POLYNOMES LAGRANGE de DEGRE 2
C           SONT DEJA CALCULES DANS PRCBL1
            DO I=1,10
               PRCBL0(I) = PRCBL1(I)
            ENDDO
            GOTO 80
         ENDIF

C        ---------------------------------------------------------------
C        CALCUL DU POINT X(tn;tn+1,Fe(bl)) c-a-d du POINT a L'INSTANT tn
C        QUI SERA AU POINT D'INTEGRATION Fe(Bl) a L'INSTANT tn+1
C        (CB2BL,CB3BL) = Fe-1(XEBL1,YEBL1)
C        ---------------------------------------------------------------
C        COORDONNEES du POINT INITIAL sur l'EF e pour les ITERATIONS n
         XEBL0 = XEBL1
         YEBL0 = YEBL1
         ZEBL0 = ZEBL1

C        NBSSPATE PAS DE TEMPS ENTRE t-DT et t pour calculer
C        X(tn;tn+1,(XEBL1,YEBL1,ZEBL1))
         DO 10 ITERT = 1, NBSSPATE

C           SAUVEGARDE DE LA VITESSE AVANT LA NOUVELLE ITERATION
            VA(1) = V(1)
            VA(2) = V(2)
            VA(3) = V(3)

C           CALCUL DES Pchapeau(XCB3BL)   POLYNOMES LAGRANGE de DEGRE 2
            PRCBL0(1) = CB1BL * ( 2.D0 * CB1BL - 1.D0 )
            PRCBL0(2) = 2.D0 * CB2BL * CB2BL - CB2BL
            PRCBL0(3) = 2.D0 * CB3BL * CB3BL - CB3BL
            PRCBL0(4) = 2.D0 * CB4BL * CB4BL - CB4BL
            D         = 4.D0 * CB1BL
            PRCBL0(5) = D    * CB2BL
            PRCBL0(6) = 4.D0 * CB2BL * CB3BL
            PRCBL0(7) = D    * CB3BL
            PRCBL0(8) = D    * CB4BL
            PRCBL0(9) = 4.D0 * CB2BL * CB4BL
            PRCBL0(10)= 4.D0 * CB3BL * CB4BL

C           CALCUL DE LA VITESSE V(tn,(CB2BL,CB3BL)) dans l'EF NEF
            V(1) = 0D0
            V(2) = 0D0
            V(3) = 0D0
            DO I=1,10
C              NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C              AU POINT X(tn;tn+1,Fe(bl))
               NDL  = NDDLNO( NUNOTE(NEF,I) - 1 )
               V(1) = V(1) +  PRCBL0(I) * VXYZPNtn( NDL+1 )
               V(2) = V(2) +  PRCBL0(I) * VXYZPNtn( NDL+2 )
               V(3) = V(3) +  PRCBL0(I) * VXYZPNtn( NDL+3 )
            ENDDO

C           LE POINT (XEBL0,YEBL0,ZEBL0) SERA TRANSPORTE en (XEBL1,YEBL1,ZEBL1)
C           SUR L'EF COURANT NEF00 de l'instant tn a tn+1 par pas ddt
C           -------------------------------------------------------------------
            XEBL0 = XEBL0 - DDT * V(1)
            YEBL0 = YEBL0 - DDT * V(2)
            ZEBL0 = ZEBL0 - DDT * V(3)

C           SAUVEGARDE DU NO DU DERNIER EF PARCOURU
            NEF0 = NEF

C           DANS QUEL ELEMENT FINI NEF du MAILLAGE EST CE POINT XYZEBL0?
            CALL XYZDSTE( NBEF,   10,     NUNOTE,
     %                    MOFACE, MXFACE, LFACES,
     %                    XYZNOE, NEF0,   XEBL0, YEBL0, ZEBL0,
     %                    NEF,    NOFAMX, IERR )
C NEF    : SI IERR=0 NUMERO DU TETRAEDRE CONTENANT LE POINT XPT,YPT,ZPT
C          SI IERR>0 NUMERO DU DERNIER TETRAEDRE PARCOURU
C NOFAMX : NUMERO (1 a 4) DE LA FACE EN CAS DE POINT EXTERIEUR AU MAILLAGE
C                  0 SINON
C IERR   : 0 PAS D'ERREUR LE POINT XEBL0, YEBL0, ZEBL0 EST DANS LE TETRAEDRE NEF
C          1 POINT EXTERIEUR AU MAILLAGE
C            PAS de TETRAEDRE derriere LA FACE NOFAMX de
C            NEF LE DERNIER TETRAEDRE DU PARCOURS POUR ATTEINDRE XYZEBL0
C          2 POINT PROCHE D'UNE FACE NON RETROUVEE DANS LE MAILLAGE
C          3 RENCONTRE D'UNE FACE DANS 3 TETRAEDRES CE QUI EST INCORRECT
C          4 NOMBRE DE TETRAEDRES PARCOURUS TROP GRAND pour DETECTION de BOUCLE

            IF( IERR .GE. 2 ) THEN

C              VITESSE AU NOEUD NONOEU AVANT SORTIE DES TETRAEDRES
               PRINT*,'f3sbop2p1: POINT:', XEBL0, YEBL0, ZEBL0,
     %                ' IERR=',IERR

               IF( IERR .EQ. 4 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT 10000, V, VA
                  ELSE
                     PRINT 20000, V, VA
                  ENDIF
10000 FORMAT(' f3sbop2p1: CAUSE sortie des TETRAEDRES: VITESSE CALCULEE=
     %'            ,3G15.6,' REMPLACEE par VITESSE PRECEDENTE=', 3G15.6)
20000 FORMAT(' f3sbop2p1: Before TETRAHEDRA EXIT: COMPUTED VELOCITY='
     %                ,3G15.6,' REPLACED by PREVIOUS VELOCITY=', 3G15.6)
               ENDIF

C              RETOUR A LA VITESSE ANTERIEURE
               V(1) = VA(1)
               V(2) = VA(2)
               V(3) = VA(3)
               IERR = 0
               GOTO 80

            ENDIF

            IF( NEF0 .NE. NEF ) THEN
C              RECUPERATION DES COORDONNEES Xe(3,4) DES SOMMETS
C              DU TETRAEDRE NEF OU SE TROUVE (XEBL0, YEBL0, ZEBL0)
               DO I=1,4
                  NS = NUNOTE(NEF,I)
                  DO K=1,3
                     Xe(K,I) = DBLE( XYZNOE( K, NS ) )
                  ENDDO
               ENDDO

C              CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 de NEF
               X1 = Xe(1,1)
               Y1 = Xe(2,1)
               Z1 = Xe(3,1)
               DELTAe = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                          Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                          Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )
               IF( DELTAe .LE. 0D0 ) THEN
C                 JACOBIEN de l'EF NEF NUL ou NEGATIF
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,*) 'f3sbop2p1: Probleme EF',NEF,
     %         ' de VOLUME*6=',DELTAe,' Be NON CALCULE pour L=',L
                  ELSE
                     WRITE(IMPRIM,*) 'f3sbop2p1: Problem FE',NEF,
     %         ' of VOLUME*6=',DELTAe,' Be is NOT COMPUTED for L=',L
                  ENDIF
                  GOTO 100
               ENDIF
            ENDIF

            IF( IERR .EQ. 1 ) THEN

C              LE POINT XYZ0 EST EXTERIEUR AUX TETRAEDRES
C              IL N'Y A PAS DE TETRAEDRE derriere LA FACE NOFAMX de
C              NEF LE DERNIER TETRAEDRE DU PARCOURS DANS xyzdste
C              ----------------------------------------------------
C              CALCUL DU POINT D'INTERSECTION DU SEGMENT
C              XEBL1,YEBL1,ZEBL1-XEBL0,YEBL0,ZEBL0 ET DE LA FACE
C              FRONTIERE NOFAMX DU TETRAEDRE NEF
               CALL INARTR( XYZ1, XYZ0,
     %                      Xe( 1, NOSOFATE(1,NOFAMX) ),
     %                      Xe( 1, NOSOFATE(2,NOFAMX) ),
     %                      Xe( 1, NOSOFATE(3,NOFAMX) ),
     %                      LINTER, PTITR, CBTR )
C              LINTER : -3 SI XYZ1-XYZ0 PAS DE CALCUL DE PTITR
C                       -2 SI XYZ1-XYZ0 EST DANS  LE PLAN DU TRIANGLE
C                       -1 SI XYZ1-XYZ0 PARALLELE AU PLAN DU TRIANGLE SANS ETRE DEDANS
C                        0 SI XYZ1-XYZ0 N'INTERSECTE PAS LE TRIANGLE ENTRE SES 3 St
C                        1 SI XYZ1-XYZ0   INTERSECTE     LE TRIANGLE ENTRE SES 3 St
C                          OU SUR UNE DE SES 3 ARETES
C              PTITR  : LES 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=0
C              CBTR   : LES 3 COORDONNEES BARYCENTRIQUES DE PT DANS LE TRIANGLE

               IF( LINTER .EQ. 1 ) THEN
C
C                 IL EXISTE UN TEL POINT D'INTERSECTION SUR LA FACE NOFAMX de NEF
C                 ---------------------------------------------------------------
                  XEBL0 = PTITR(1)
                  YEBL0 = PTITR(2)
                  ZEBL0 = PTITR(3)

ccc            WRITE(IMPRIM,*)'f3sbop2p1: Point frontiere intersecte',
ccc     %                      XEBL0, YEBL0, ZEBL0,' Vitesse',V,' EF',NEF
ccc            WRITE(IMPRIM,*)
C
               ELSE
C
C                 IL N'EXISTE PAS DE POINT PROJETE SUR UNE FACE FRONTIERE
C                 => POINT EXTERIEUR AU DOMAINE
C                 -------------------------------------------------------
C                 VITESSE AU POINT AVANT SORTIE DES TETRAEDRES
                  V(1) = VA(1)
                  V(2) = VA(2)
                  V(3) = VA(3)
                  GOTO 80

               ENDIF
C
            ENDIF

C           CALCUL DES 4 COORDONNEES BARYCENTRIQUES DU
C           POINT (XEBL0,YEBL0,ZEBL0) DANS LE TETRAEDRE NEF QUI LE CONTIENT
C           CE SONT ABSCISSE,ORDONNEE,COTE DU POINT DANS L'EF DE REFERENCE

C           Lambda 1
            CB1BL = DETM33( Xe(1,2)-XEBL0, Xe(1,3)-XEBL0, Xe(1,4)-XEBL0,
     %                      Xe(2,2)-YEBL0, Xe(2,3)-YEBL0, Xe(2,4)-YEBL0,
     %                      Xe(3,2)-ZEBL0, Xe(3,3)-ZEBL0, Xe(3,4)-ZEBL0)
     %            / DELTAe
C           Lambda 2
            CB2BL = DETM33( Xe(1,3)-XEBL0, Xe(1,1)-XEBL0, Xe(1,4)-XEBL0,
     %                      Xe(2,3)-YEBL0, Xe(2,1)-YEBL0, Xe(2,4)-YEBL0,
     %                      Xe(3,3)-ZEBL0, Xe(3,1)-ZEBL0, Xe(3,4)-ZEBL0)
     %            / DELTAe
C           Lambda 3
            CB3BL = DETM33( Xe(1,4)-XEBL0, Xe(1,1)-XEBL0, Xe(1,2)-XEBL0,
     %                      Xe(2,4)-YEBL0, Xe(2,1)-YEBL0, Xe(2,2)-YEBL0,
     %                      Xe(3,4)-ZEBL0, Xe(3,1)-ZEBL0, Xe(3,2)-ZEBL0)
     %            / DELTAe
C           Lambda 4
            CB4BL = DETM33( Xe(1,1)-XEBL0, Xe(1,3)-XEBL0, Xe(1,2)-XEBL0,
     %                      Xe(2,1)-YEBL0, Xe(2,3)-YEBL0, Xe(2,2)-YEBL0,
     %                      Xe(3,1)-ZEBL0, Xe(3,3)-ZEBL0, Xe(3,2)-ZEBL0)
     %            / DELTAe

C           Verification (XEBL0,YEBL0,ZEBL0) INTERNE au TETRAEDRE NEF?
            D = ABS(CB1BL) + ABS(CB2BL) + ABS(CB3BL) + ABS(CB4BL)
            IF( D .GT. 1.01D0 ) THEN

ccc               print*,'f3sbop2p1: Pb NEF=',NEF,' L=',L,' CB1234=',D,
ccc     %                ' CB1BL=',CB1BL,' CB2BL=',CB2BL,' CB3BL=',CB3BL,
ccc     %                ' CB4BL=',CB4BL,'->  XYZBL0 NON INTERNE a NEF'
ccc               print*,'f3sbop2p1: XYZBL0=',XYZ0,
ccc     %                ' est PROJETE dans NEF.  LINTER=',LINTER

C              PROJECTION DU POINT XYZBL0 DANS LE TETRAEDRE NEF
               IF( CB1BL .LT. 0D0 ) CB1BL=0D0
               IF( CB2BL .LT. 0D0 ) CB2BL=0D0
               IF( CB3BL .LT. 0D0 ) CB3BL=0D0
               IF( CB4BL .LT. 0D0 ) CB4BL=0D0
               D = CB1BL + CB2BL + CB3BL + CB4BL
               CB1BL = CB1BL / D
               CB2BL = CB2BL / D
               CB3BL = CB3BL / D
               CB4BL = CB4BL / D

C         (XEBL0,YEBL0,ZEBL0)=Fe(CB2BL,CB3BL,CB4BL) SUR LE TETRAEDRE NEF
           XEBL0=CB1BL*Xe(1,1)+CB2BL*Xe(1,2)+CB3BL*Xe(1,3)+CB4BL*Xe(1,4)
           YEBL0=CB1BL*Xe(2,1)+CB2BL*Xe(2,2)+CB3BL*Xe(2,3)+CB4BL*Xe(2,4)
           ZEBL0=CB1BL*Xe(3,1)+CB2BL*Xe(3,2)+CB3BL*Xe(3,3)+CB4BL*Xe(3,4)

ccc               print*,'f3sbop2p1: XYZBL0=',XYZ0,' PROJETE DANS NEF'

               IF( LINTER .EQ. 1 ) GOTO 30
            ENDIF

 10      ENDDO

ccc      print*,'XYZtn+1=',XEBL1,YEBL1,ZEBL1,' XYZtn=',XEBL0,YEBL0,ZEBL0,
ccc     %       '  V=',V,' EF=',NEF,' CB=',CB1BL,CB2BL,CB3BL,CB4BL

C        Calcul de X(tn;tn+1,Fe(bl)) effectue:
C        Les Pchapeau(X(tn;tn+1,Fe(bl))) POLYNOMES LAGRANGE P2
C        DANS LE TETRAEDRE NEF contenant (XEBL0,YEBL0,ZEBL0) sont calcules
C        =================================================================
C        CALCUL DES Pchapeau(XCB3BL)   POLYNOMES LAGRANGE de DEGRE 2
 30      PRCBL0(1) = CB1BL * ( 2.D0 * CB1BL - 1.D0 )
         PRCBL0(2) = 2.D0 * CB2BL * CB2BL - CB2BL
         PRCBL0(3) = 2.D0 * CB3BL * CB3BL - CB3BL
         PRCBL0(4) = 2.D0 * CB4BL * CB4BL - CB4BL
         D         = 4.D0 * CB1BL
         PRCBL0(5) = D    * CB2BL
         PRCBL0(6) = 4.D0 * CB2BL * CB3BL
         PRCBL0(7) = D    * CB3BL
         PRCBL0(8) = D    * CB4BL
         PRCBL0(9) = 4.D0 * CB2BL * CB4BL
         PRCBL0(10)= 4.D0 * CB3BL * CB4BL

C        CALCUL DE LA VITESSE V(tn,(CB2BL,CB3BL,CB4BL)) dans l'EF NEF
         V(1) = 0D0
         V(2) = 0D0
         V(3) = 0D0
         DO I=1,10
C           NO DU DL I DANS L'EF NEF DE LA VITESSE AU TEMPS tn
C           AU POINT X(tn;tn+1,Fe(bl))
            NDL  = NDDLNO( NUNOTE(NEF,I) - 1 )
            V(1) = V(1) +  PRCBL0(I) * VXYZPNtn( NDL+1 )
            V(2) = V(2) +  PRCBL0(I) * VXYZPNtn( NDL+2 )
            V(3) = V(3) +  PRCBL0(I) * VXYZPNtn( NDL+3 )
         ENDDO

cccC        PAS D'ECRETAGE EVENTUEL SUR LA VITESSE MAX ANTERIEURE * 1.5
cccC        -----------------------------------------------------------
ccc 80      VNORM = SQRT( V(1)**2 + V(2)**2 + V(3)**2 )
ccc         D     = VNORM / VITMXtn
ccc         IF( D .GT. 1.5D0 ) THEN
ccc            print*,'f3sbop2p1: Attention ECRETAGE de VNORM=',VNORM,
ccc     %             ' en 1.5 VITMXtn=',1.5D0*VITMXtn,' V=',V
ccc            print*,'f3sbop2p1: NEF=',NEF,
ccc     %             ' CBT=',ABS(CB1BL)+ABS(CB2BL)+ABS(CB3BL)+ABS(CB4BL),
ccc     %             ' CBBL=',CB1BL,CB2BL,CB3BL,CB4BL
ccc            print*,'f3sbop2p1: PRCBL0T=',
ccc     %              PRCBL0(1)+PRCBL0(2)+PRCBL0(3)+PRCBL0(4)+PRCBL0(5)+
ccc     %              PRCBL0(6)+PRCBL0(7)+PRCBL0(8)+PRCBL0(9)+PRCBL0(10),
ccc     %             ' PRCBL0=',PRCBL0
ccc            D = 1.5D0 * VITMXtn / VNORM
ccc            V(1) = V(1) * D
ccc            V(2) = V(2) * D
ccc            V(3) = V(3) * D
ccc         ENDIF

C        CALCUL DU VECTEUR BE DE L'EF NEF00
C        ----------------------------------
C        POIDS * JACOBIEN AU POINT L D'INTEGRATION * CoefTr
 80      RDPOIDS = POTEP5(L) * DELTAe1 * CoefTr

         IF( TempouVit .EQ. 0 ) THEN

C           ICI LA TEMPERATURE AU POINT X(tn;tn+1,Fe(bl))
            TEMPERA = 0D0
            DO I=1,10
C              NO DU DL I DANS L'EF NEF DE LA TEMPERATURE AU TEMPS tn
C              AU POINT X(tn;tn+1,Fe(bl))
               NDL = NUNOTE(NEF,I)
               TEMPERA = TEMPERA + PRCBL0(I) * TEMPER(NDL)
            ENDDO

            D = RDPOIDS * TEMPERA
            DO I=1,10
C              CoefTr Poids Jacobien tP(bl) TEMPER(X(tn;tn+1,Fe(bl))
C              AU POINT X(tn;tn+1,Fe(bl))
               BE(I) = BE(I) + D * PRCBL1(I)
            ENDDO

         ELSE IF( TempouVit .EQ. 1 ) THEN

C           V(1:3) EST ICI LA VITESSE AU POINT X(tn;tn+1,Fe(bl))
C           CONSTRUCTION DU BLOC ASSOCIE A Vk DU SECOND MEMBRE ELEMENTAIRE
            IK = 0
            DO K=1,3
C           V(k) COMPOSANTE k de la VITESSE
               D = RDPOIDS * V(K)
               DO I=1,10
C                 COMPOSANTE IK du VECTEUR
C                 CoefTr Poids Jacobien tP(bl) Vk( X(tn;tn+1,Fe(bl)) )
                  BE(IK+I) = BE(IK+I) + D * PRCBL1(I)
               ENDDO
               IK = IK + 10
            ENDDO

C        ELSE IF( TempouVit .EQ. 2 ) THEN
C           CoefTr Integrale tP2 ( TEMPER(tn,x)
C                                 -u1(tn,x) d TEMPER(tn,x)/dx1
C                                 -u2(tn,x) d TEMPER(tn,x)/dx2
C                                 -u3(tn,x) d TEMPER(tn,x)/dx3 ) dx
C           N'EST PAS CALCULE ICI. A FAIRE

         ENDIF

C        FIN DE L'INTEGRATION DU POINT L
 100  ENDDO

 9999 RETURN
      END
