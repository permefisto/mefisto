      SUBROUTINE TRLLD3( XYZ1DR, XYZ2DR, NDIM, NTDL, NCAS0,NCAS1,
     %                   NTYP,   TEMPER, dptemp,
     %                   NUTYEL, NBELEM,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL,
     %                   NBCOOR, NBPOIT, XYZPOI,
     %                   MXSOMM, SOLEL,  COPOE,
     %                   MXPILE, LAPILE,
     %                   XYZSEF, XYZTEF, FBASE,
     %                   MXSGPT, NBSGPT, XYZSEG, TEMSEG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER EN 3D LES COUPLES DE POINTS D'INTERSECTION
C -----    D'UNE DROITE AVEC LES SOUS-TETRAEDRES
C          DES ELEMENTS FINIS D'UN TYPE DONNE
C          SOIT TETRAEDRE SOIT PENTAEDRE SOIT HEXAEDRE
C
C ENTREES:
C --------
C XYZ1DR : SOMMET 1 DE LA DROITE
C XYZ2DR : SOMMET 2 DE LA DROITE
C NDIM   : ESPACE DE TRAVAIL 2 OU 3
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NCAS0  : NUMERO DU PREMIER CAS A TRAITER
C NCAS1  : NUMERO DU DERNIER CAS A TRAITER
C TEMPER : LES NTDL TEMPERATURES DES NCAS0 a NCAS1 SUR LE MAILLAGE
C
C NUTYEL : NUMERO DU TYPE D'EF A TRAITER ICI
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT  FINI  DE CE TYPE
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS FINIS DE CE TYPE
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT  FINI  DE CE TYPE
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS FINIS DE CE TYPE
C
C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT (3 ou 6)
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXSOMM : NOMBRE MAXIMUM DE SOMMETS DES SOUS-TETRAEDRES DE L'EF REFERENCE
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C MXPILE : NBRE MAXIMUM DE SOUS-TETRAEDRES DANS LA PILE LAPILE
C LAPILE : PILE DES 4 SOMMETS DES SOUS-TETRAEDRES
C XYZSEF : XYZSEF(1:3,*) 3 COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES
C                        DANS L'EF DE REFERENCE
C XYZTEF : XYZTEF(1:3,*) 3 COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES
C                        DANS L'EF COURANT
C          XYZTEF(4,*) VALEUR DE LA SOLUTION EN CE MEME SOMMET
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)
C
C SORTIES:
C --------
C NBSGPT : NOMBRE TOTAL DE COUPLES XYZV INTERSECTION AVEC LES TETRAEDRES
C XYZSEG : XYZ         DES SOMMETS DES SEGMENTS DE LA DROITE D'INTERSECTION
C TEMSEG : TEMPERATURE DES SOMMETS DES SEGMENTS DE LA DROITE D'INTERSECTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS OCTOBRE 2003
C MODIFS :ALAIN PERRONNET LJLL UPMC & St Pierre du Perray   Fevrier 2013
C23456---------------------------------------------------------------012
C     SEUIL DE PRECISION
      PARAMETER    ( EPS=0.00001,UNPEPS=1.0+EPS)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NUNDEL(NBELEM,NBNOEL),
     %                  NUPTEL(NBELEM,NBPOEL),
     %                  LAPILE(4,MXPILE)
      REAL              XYZPOI(NBCOOR,NBPOIT)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp

      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1),
     %                  SOLEL(NBNOEL,NCAS0:NCAS1),
     %                  COPOE(NBPOEL,NDIM)
      DOUBLE PRECISION  XD,YD,ZD,PROSCD,
     %                  FBASE(NBNOEL,3)
      DOUBLE PRECISION  XYZTEF(4,MXSOMM)
      REAL              XYZSEF(4,MXSOMM)
      REAL              COSOTE(3,4)
      DOUBLE PRECISION  XYZ1DR(3), XYZ2DR(3)
      DOUBLE PRECISION  PT(3), COBARY(3), TEMTRI(3)
      INTEGER           NOSOTR(3)
      REAL              XYZSEG(3,2,MXSGPT)
      REAL              TEMSEG(2,NCAS0:NCAS1,MXSGPT)
      REAL              XYZ(3)
C
C     DONNEES DES SOUS-TETRAEDRES D'UN TETRAEDRE ou PYRAMIDE ou PENTAEDRE
C     ou HEXAEDRE SUBDIVISE REGULIEREMENT
      include "./incl/nostst.inc"
C
C     CALCUL DU NOMBRE DE SUBDIVISIONS DE (0,1)
C     EN FONCTION DU NOMBRE DE NOEUDS DE L ELEMENT FINI
      IF( NBNOE .LE. 8 ) THEN
         NBSOUI = 1
      ELSE
         NBSOUI = 2
      ENDIF
      NBSTTE = 2 * ( NBSOUI * NBSOUI + 1 )
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
      CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &              NBINVA, NUINVA, NUINTI, NBNDIN )
C     LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
      NOINTE = NUINTI(1)
C
CCC      IDEM POUR LE GRADIENT DE TEMPERATURE
CCC      CALL ELINT1( 'THERMIQUE', NUTYEL, NDIMF, NOINFD, NBPOF,
CCC     &              NBINVA, NUINVA, NUINTI, NBNDIN, NUNOIN )
C
C     LA DIMENSION DE L'ESPACE DES COORDONNEES  EST NDIM
      IF( NDIM .NE. NDIMF ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'INCOMPATIBILITE DIMENSION ESPACES pour F et X'
         ELSE
            KERR(1) = 'INCOMPATIBILITY DIMENSION SPACES for F and X'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NFACE .LT. 4  .OR.  NFACE .GT. 6 ) THEN
C        ELEMENT NON TETRAEDRE OU PENTAEDRE OU HEXAEDRE
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: ELEMENT FINI DIFFERENT DE'
            KERR(2) = 'TETRAEDRE ou PENTAEDRE ou HEXAEDRE'
            KERR(3) = 'PLAN-VALEURS NON CALCULABLES'
         ELSE
            KERR(1) = 'ERROR: FINITE ELEMENT DIFFERENT of'
            KERR(2) = 'TETRAHEDRON or PENTAHEDRON or HEXAHEDRON'
            KERR(3) = 'PLANE-VALUES NOT COMPUTABLE'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     FORMATION DES TABLEAUX SUR L'ELEMENT FINI DE REFERENCE
C     ======================================================
C     OUVERTURE DE LA PILE DU NO DES 4 SOMMETS DES SOUS-TETRAEDRES
C     LHPIL0 POINTE SUR LE SOMMET DE LA PILE LAPILE
      LHPIL0 = 0
C     NOSOMM NO DU DERNIER SOMMET GENERE DANS LES SOUS-TETRAEDRE
      NOSOMM = 0
C
C     GENERATION DES SOUS-TETRAEDRES A L'INTERIEUR DE L'EF DE REFERENCE
C     =================================================================
C     LE NO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C     LES COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES DANS L'EF REFERENCE
C     SONT EMPILES DANS LAPILE(4,MXPILE) ET XYZSEF(1:4,MXSOMM)
C
      IF( NFACE .EQ. 4 ) THEN
C
C        *************
C        * TETRAEDRE *
C        *************
C
C        LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C        --------------------------------------------
         CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
         IF( LHPIL0 .LT. 0 ) GOTO 9999
C
C        LES 3 COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C        -----------------------------------------------------------------
         CALL COPTTE( TETR, NBSOUI+1, 4, XYZSEF(1,NOSOMM+1) )
C
         NOSOMM = NOSOMM + NBSTTE
C
      ELSE IF( NFACE .EQ. 5 ) THEN
C
         IF( NBNSOM .LE. 5 ) THEN
C
C           *************
C           * PYRAMIDE  *
C           *************
C           BOUCLE SUR 2 SOUS-TETRAEDRES PRIMAIRES
C           ======================================
            DO 14 I = 1, 2
C
C              LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C              --------------------------------------------
               CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
               IF( LHPIL0 .LT. 0 ) GOTO 9999
C
               DO 12 J = 1, 4
                  NO = NOSY(J, I)
                  COSOTE(1, J) = PYRA(1, NO)
                  COSOTE(2, J) = PYRA(2, NO)
                  COSOTE(3, J) = PYRA(3, NO)
 12            CONTINUE
C
C              LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C              ---------------------------------------------------------------
               CALL COPTTE( COSOTE, NBSOUI+1, 4, XYZSEF(1,NOSOMM+1) )
C
               NOSOMM = NOSOMM + NBSTTE
 14         CONTINUE
C
         ELSE
C
C           *************
C           * PENTAEDRE *
C           *************
C           BOUCLE SUR 3 SOUS-TETRAEDRES PRIMAIRES
C           ======================================
            DO 18 I = 1, 3
C
C              LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C              --------------------------------------------
               CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
               IF( LHPIL0 .LT. 0 ) GOTO 9999
C
               DO 16 J = 1, 4
                  NO = NOSP(J, I)
                  COSOTE(1, J) = PENT(1, NO)
                  COSOTE(2, J) = PENT(2, NO)
                  COSOTE(3, J) = PENT(3, NO)
 16            CONTINUE
C
C              LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C              ---------------------------------------------------------------
               CALL COPTTE( COSOTE, NBSOUI+1, 4, XYZSEF(1,NOSOMM+1) )
C
               NOSOMM = NOSOMM + NBSTTE
 18         CONTINUE
         ENDIF
C
      ELSE
C
C        ************
C        * HEXAEDRE *
C        ************
C
C        BOUCLE SUR 5 SOUS-TETRAEDRES PRIMAIRES
C        ======================================
         DO 26 I = 1, 5
C
C           LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C           --------------------------------------------
            CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
            IF( LHPIL0 .LT. 0 ) GOTO 9999
C
            DO 24 J = 1, 4
               NO = NOSH(J, I)
               COSOTE(1, J) = HEXA(1, NO)
               COSOTE(2, J) = HEXA(2, NO)
               COSOTE(3, J) = HEXA(3, NO)
 24         CONTINUE
C
C           LES 3 COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C           -----------------------------------------------------------------
            CALL COPTTE( COSOTE, NBSOUI+1, 4, XYZSEF(1,NOSOMM+1) )
C
            NOSOMM = NOSOMM + NBSTTE
 26      CONTINUE
      ENDIF
C
C     ===============================
C     LA BOUCLE SUR LES EF DE CE TYPE
C     ===============================
      DO 1000 NEF=1,NBELEM
C
C        FORMATION DES TABLEAUX SUR L'EF COURANT
C        SOLEL(NBNOEL,NCAS0:NCAS1) TEMPERATURE AUX NBNOEL DE L'EF COURANT
C        COPOE(NBPOEL,NDIM) COORDONNEES DES POINTS DE L'EF COURANT
C        ================================================================
C        EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER
         DO 35 I=1,NBNOEL
            J = NUNDEL(NEF,I)
            DO NCAS = NCAS0, NCAS1
               IF( NTYP .EQ. 0 ) THEN
                  SOLEL( I, NCAS ) = TEMPER( J, NCAS )
               ELSE
                  SOLEL( I, NCAS ) = dptemp( NCAS )%dptab( J )
               ENDIF
            ENDDO
 35      CONTINUE
C
C        EXTRACTION DES COORDONNEES DES POINTS DE L'ELEMENT FINI
         DO 38 I=1,NBPOEL
            DO 37 J=1,NDIM
               COPOE(I,J) = XYZPOI( J, NUPTEL(NEF,I) )
 37         CONTINUE
 38      CONTINUE
C
C        CALCUL DE XYZ ET LA SOLUTION AUX SOMMETS DES SOUS-TETRAEDRES
C        ============================================================
C        BOUCLE SUR LES SOMMETS DES SOUS-TETRAEDRES
C        ------------------------------------------
         DO 80 I=1,NOSOMM
C
C           LES 3 COORDONNEES DU POINT I DANS L'EF DE REFERENCE
            XD = XYZSEF(1,I)
            YD = XYZSEF(2,I)
            ZD = XYZSEF(3,I)
C           LA VALEUR DES NBN FONCTIONS DE BASE EN (XD,YD,ZD)
            CALL INTERP( NOINTE, XD, YD, ZD,  NBN, FBASE )
C
C           XYZ DU SOMMET DANS L'EF COURANT
            XYZTEF(1,I) = PROSCD( FBASE, COPOE(1,1), NBN )
            XYZTEF(2,I) = PROSCD( FBASE, COPOE(1,2), NBN )
            XYZTEF(3,I) = PROSCD( FBASE, COPOE(1,3), NBN )

 80      CONTINUE
C
C        LE SOMMET DE LA PILE AU DEBUT DE CET EF
         LHPILE = LHPIL0
C
C        ***********************************************************
C        TANT QUE  LA PILE DES SOUS-TETRAEDRES N'EST PAS VIDE  FAIRE
C        ***********************************************************
 100     IF( LHPILE .LE. 0 ) GOTO 1000
C
C           LES 4 SOMMETS DU SOUS-TETRAEDRE EN SOMMET DE PILE SONT DEPILES
C           NOMBRE DE POINTS D'INTERSECTION SOUS-TETRAEDRE-DROITE
            NBPTIN = 0
            NBSGIN = NBSGPT + 1
C
C           BOUCLE SUR LES 4 TRIANGLES FACES DU SOUS-TETRAEDRE
            DO 200 NS1=1,4
C
               NS2 = MOD(NS1,4) + 1
               NS3 = MOD(NS2,4) + 1
C              NUMERO DES 3 SOMMETS DE LA FACE NS1 DANS XYZTEF
               NOSOTR(1) = LAPILE(NS1,LHPILE)
               NOSOTR(2) = LAPILE(NS2,LHPILE)
               NOSOTR(3) = LAPILE(NS3,LHPILE)
C              LE TRIANGLE NS1 DU TETRAEDRE EST IL INTERSECTE PAR LA DROITE?
               CALL INDRPL( XYZ1DR, XYZ2DR,
     %                      XYZTEF(1,NOSOTR(1)),
     %                      XYZTEF(1,NOSOTR(2)),
     %                      XYZTEF(1,NOSOTR(3)),
     %                      PT, NERR )
C              PT     : LES 3 COORDONNEES DU POINT D'INTERSECTION SI NERR=0
C              NERR   : 1 SI LA DROITE EST PARALLELE AU PLAN
C                       2 SI S1=S2
C                       0 SINON
C
               IF( NERR .EQ. 0 ) THEN
C                 LE POINT PT EST IL DANS LE TRIANGLE DE SOMMETS NOSOTR ?
                  CALL PTDSTR( PT, XYZTEF, NOSOTR, COBARY, NONOUI )
C                 COBARY : LES 3 COORDONNEES BARYCENTRIQUES
C                         (/S3 /S1 /S2) SI TRIANGLE
C                          NON DEGENERE ET NON INITIALISE SINON
C                 NONOUI : 1 SI LE POINT EST DANS OU SUR LE TRIANGLE
C                          0 SI LE POINT P EST EXTERIEUR AU TRIANGLE
C                         -1 SI LE TRIANGLE EST DEGENERE
C
                  IF( NONOUI .GT. 0 ) THEN
C                    PT EST UN POINT D'INTERSECTION
C                    EXISTE T IL DEJA?
                     XYZ(1) = REAL( PT(1) )
                     XYZ(2) = REAL( PT(2) )
                     XYZ(3) = REAL( PT(3) )
                     DO 110 N=NBPTIN,1,-1
                        CALL XYZIDE( XYZ, XYZSEG(1,N,NBSGIN), NONOUI )
                        IF( NONOUI .GT. 0 ) GOTO 200
 110                 CONTINUE
C                    NON. LE POINT D'INTERSECTION EST AJOUTE
                     NBPTIN = NBPTIN + 1
                     IF( NBPTIN .LE. 2 ) THEN

C                       LES 3 COORDONNEES DU POINT D'INTERSECTION
                        XYZSEG(1,NBPTIN,NBSGIN) = REAL( PT(1) )
                        XYZSEG(2,NBPTIN,NBSGIN) = REAL( PT(2) )
                        XYZSEG(3,NBPTIN,NBSGIN) = REAL( PT(3) )

C                       LA VALEUR DES NCAS0:NCAS1 TEMPER EN CE POINT
                        DO K=1,3
C                          LES 3 COORDONNEES DU POINT NOSOTR(K) DANS L'EF DE REFERENCE
                           N = NOSOTR(K)
                           XD = XYZSEF(1,N)
                           YD = XYZSEF(2,N)
                           ZD = XYZSEF(3,N)
C                          LA VALEUR DES NBN FONCTIONS DE BASE EN (XD,YD,ZD)
                           CALL INTERP( NOINTE, XD, YD, ZD,
     %                                  NBN, FBASE(1,K) )
                        ENDDO

                        DO NCAS = NCAS0, NCAS1
                           DO K=1,3
C                             LA VALEUR DE LA SOLUTION AU SOMMET K DE NOSOTR(1:3)
C                            (MEME VALEUR DE LA SOLUTION SUR L'EF REFERENCE OU COURANT)
                              TEMTRI(K) = PROSCD( FBASE(1,K),
     %                                            SOLEL(1,NCAS), NBN )
                           ENDDO
                           TEMSEG( NBPTIN, NCAS, NBSGIN ) = REAL(
     %                                        COBARY(1) * TEMTRI(3)
     %                                      + COBARY(2) * TEMTRI(1)
     %                                      + COBARY(3) * TEMTRI(2) )
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
 200        CONTINUE
C
C           BILAN SUR LE NOMBRE DE POINTS D'INTERSECTION
            IF( NBPTIN .EQ. 2 ) NBSGPT = NBSGIN
            IF( NBSGPT .GE. MXSGPT ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'TRLLD3: TROP de POINTS d''INTERSECTION'
               ELSE
                  KERR(1) = 'TRLLD3: TOO MANY POINTS of INTERSECTION'
               ENDIF
               CALL LEREUR
               GOTO 9999
            ENDIF
C
C           REMONTEE POUR TRAITER LE SOUS-TETRAEDRE SUIVANT
C           ===============================================
            LHPILE = LHPILE - 1
            GOTO 100
C
C        *************************************************
C        FIN DU TRAITEMENT DES SOUS-TETRAEDRES DE L'EF NEF
C        *************************************************
 1000 CONTINUE
C
 9999 RETURN
      END
