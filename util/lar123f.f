      SUBROUTINE LAR123F( NBSOMM,    XYZSOM,
     %                    N1SOEF,    NBEFOB,    NOSOEF,
     %                    L1ARET,    L2ARET,    LARETE,
     %                    MXAR1F,    NBAR1F,    MNAR1F,
     %                    MXAR2F,    NBAR2F,    MNAR2F,
     %                    MXAR3F,    NBAR3F,    MNAR3F,
     %                    MXA2SF,    NBA2SF,    MNA2SF,
     %                    MXAR2FCOL, NBAR2FCOL, MNAR2FCOL, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : A PARTIR DU TABLEAU LARETE DES ARETES CONSTRUCTION DES LISTES
C ----- D'ARETES APPARTENANT A 1 FACE  (TRIANGLE ou QUADRANGLE)
C       D'ARETES APPARTENANT A 2 FACES (TRIANGLE ou QUADRANGLE)
C       D'ARETES APPARTENANT A 3 FACES (TRIANGLE ou QUADRANGLE)
C       D'ARETES APPARTENANT A 2 FACES DONT LE RAPPORT DES SURFACES
C                                      EST <0.0667 TROP FAIBLE
C       D'ARETES APPARTENANT A 2 FACES COLLEES C-A-D
C       LES 2 FACES SONT COPLANAIRES DU MEME COTE DE LEUR ARETE COMMUNE

C ENTREES:
C --------
C NBSOMM : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : XYZ DES NBSOMM SOMMETS DU MAILLAGE
C N1SOEF : NOMBRE DE SOMMETS DECLARABLES DANS NOSOEF PAR EF
C NBEFOB : NOMBRE D'EF ACTIFS DANS LE TABLEAU NOSOEF
C NOSOEF : NUMERO XYZSOM DES N1SOEF SOMMETS DES NBEFOB EF

C L1ARET : NOMBRE DE MOTS PAR ARET DU TABLEAU LARETE
C L2ARET : NOMBRE DE FACES DU TABLEAU LARETE
C LARETE : LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          LARETE(4:3+MXFAAR,I)= NO NOSOFA DE LA FACE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE MXFAAR FACES, 
C          LE NUMERO DE FACE MXFAAR EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES FACES EST INCOMPLETE

C MODIFIES:
C --------
C MXAR1F : NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC AR1F
C NBAR1F : NOMBRE         D'ARETES DE NO LARETE DANS LE TMC AR1F
C MNAR1F : >0 ADRESSE MCN DU TMC NAR1F NO LARETE DE L'ARETE DES
C             ARETES APPARTENANT A 1 FACE
C          =0 SI NBAR1F=0

C MXAR2F : NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC AR2F
C NBAR2F : NOMBRE         D'ARETES DE NO LARETE DANS LE TMC AR2F
C MNAR2F : >0 ADRESSE MCN DU TMC NAR2F NO LARETE DE L'ARETE DES
C             ARETES APPARTENANT A 2 FACES
C          =0 SI NBAR2F=0

C MXAR3F : NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC AR3F
C NBAR3F : NOMBRE         D'ARETES DE NO LARETE DANS LE TMC AR3F
C MNAR3F : >0 ADRESSE MCN DU TMC NAR3F NO LARETE DE L'ARETE DES
C             ARETES APPARTENANT A 3 FACES
C          =0 SI NBAR3F=0

C MXA2SF : NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC A2SF
C NBA2SF : NOMBRE         D'ARETES DE NO LARETE DANS LE TMC A2SF
C MNA2SF : >0 ADRESSE MCN DU TMC NA2SF NO LARETE DES ARETES DE RATIO FAIBLE
C          =0 SI NBA2SF=0

C MXAR2FCOL: NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC AR2FCOL
C NBAR2FCOL: NOMBRE         D'ARETES DE NO LARETE DANS LE TMC AR2FCOL
C MNAR2FCOL: >0 ADRESSE MCN DU TMC NAR2FCOL NO LARETE DE L'ARETE DES
C               COUPLES DE TRIANGLES ADJACENTS COLLES
C            =0 SI NBAR2FCOL=0
C IERR     : =0 PAS D'ERREUR DETECTEE
C            >0 SATURATION DE LA MEMOIRE CENTRALE MCN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY           Janvier 2021
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           LARETE(L1ARET,L2ARET), NOSOEF(N1SOEF,NBEFOB)
      REAL              XYZSOM(3,NBSOMM)

      IERR = 0

C     SI LE TABLEAU AR1F EXISTE ALORS IL EST DETRUIT
      IF( MNAR1F .GT. 0 ) CALL TNMCDS('ENTIER', MXAR1F, MNAR1F)
      NBAR1F = 0
      MXAR1F = 0

C     SI LE TABLEAU AR2F EXISTE ALORS IL EST DETRUIT
      IF( MNAR2F .GT. 0 ) CALL TNMCDS('ENTIER', MXAR2F, MNAR2F)
      NBAR2F = 0
      MXAR2F = 0

C     SI LE TABLEAU AR3F EXISTE ALORS IL EST DETRUIT
      IF( MNAR3F .GT. 0 ) CALL TNMCDS('ENTIER', MXAR3F, MNAR3F)
      NBAR3F = 0
      MXAR3F = 0

C     SI LE TABLEAU A2SF EXISTE ALORS IL EST DETRUIT
      IF( MNA2SF .GT. 0 ) CALL TNMCDS( 'ENTIER', MXA2SF, MNA2SF )
      NBA2SF = 0
      MXA2SF = 0

C     SI LE TABLEAU AR2FCOL EXISTE ALORS IL EST DETRUIT
      IF( MNAR2FCOL .GT. 0 ) CALL TNMCDS('ENTIER', MXAR2FCOL, MNAR2FCOL)
      NBAR2FCOL = 0
      MXAR2FCOL = 0

C     RECHERCHE DU NOMBRE DE FACES DES ARETES
C     ET CELLES DE RAPPORT FAIBLE DES SURFACES DES TRIANGLES ADJACENTS
C     OU COLLEES
      DO 100 NA = 1, L2ARET

         IF( LARETE(1,NA) .NE. 0 ) THEN
C           L'ARETE NA EST INITIALISEE

C           NOMBRE DE FACES DE CETTE ARETE NA
            NBF = 0
            DO K = 4, L1ARET
               NF1 = ABS( LARETE(K,NA) )
               IF( NF1 .GT. 0 ) THEN
C                 UNE FACE DE PLUS
                  NBF = NBF + 1
               ENDIF
            ENDDO

C           L'ARETE APPARTIENT APPARTIENT A NBF FACES
            GOTO( 100, 10, 20, 30, 30, 30 ),1+NBF

C           L'ARETE NA APPARTIENT A 1 FACE
C           ------------------------------
 10         NBAR1F = NBAR1F + 1
            IF( NBAR1F .GT. MXAR1F ) THEN
               IF( MXAR1F .EQ. 0 ) THEN
C                 LE TABLEAU AR1F EST CREE
                  MXAR1F = NBSOMM
                  CALL TNMCDC( 'ENTIER', MXAR1F, MNAR1F )
               ELSE
C                 LE TABLEAU AR1F EST RALLONGE
                  CALL TNMCAU( 'ENTIER', MXAR1F, MXAR1F*2, MXAR1F,
     %                          MNAR1F )
                  MXAR1F = MXAR1F * 2
               ENDIF
               IF( MNAR1F .LE. 0 ) THEN
                  IERR = 1
                  GOTO 9999
               ENDIF
            ENDIF

C           NUMERO LARETE DE L'ARETE SIMPLE NBAR1F
            MCN( MNAR1F-1+NBAR1F ) = NA
            GOTO 100

C           2 TRIANGLES SONT ADJACENTS A CETTE ARETE
C           ----------------------------------------

C           LES 2 FACES TRIANGULAIRES SONT ELLES COLLEES?
C           CALCUL DU COSINUS DE L'ANGLE DIEDRE DES 2 TRIANGLES
C           ADJACENTS PAR L'ARETE NA
C           ---------------------------------------------------
 20         NF1 = ABS( LARETE( 4, NA ) )
            NS1 = LARETE( 1, NA )
            NS2 = LARETE( 2, NA )

C           RECHERCHE DU 3-EME SOMMET DE NF1
            DO K = 1, N1SOEF
               NS3 = NOSOEF( K, NF1 )
               IF( NS3 .NE. NS1 .AND. NS3 .NE. NS2 ) GOTO 22
            ENDDO
            GOTO 100

C           RECHERCHE DU 3-EME SOMMET DE NF2
 22         N = 4
            N = N + 1
            NF2 = ABS( LARETE( N, NA ) )
            IF( NF2 .EQ. 0 ) GOTO 22
            DO K = 1, N1SOEF
               NS4 = NOSOEF( K, NF2 )
               IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 24
            ENDDO
            GOTO 100

C           COSINUS DE L'ANGLE DIEDRE DES 2 TRIANGLES ADJACENTS
 24         CALL COS2TR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %                   XYZSOM(1,NS2), XYZSOM(1,NS1), XYZSOM(1,NS4),
     %                   COS2T, IERR1, IERR2 )
            IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) GOTO 25

            IF( COS2T .LT. -0.99 ) THEN

C              2 TRIANGLES COLLES
C              ------------------
C              COANPL: SEUIL DU COSINUS DE L'ANGLE FORME PAR LES
C                      NORMALES AUX 2 TRIANGLES ET AU DESSUS DUQUEL
C                      ILS SONT CONSIDERES COPLANAIRES
C                 ( 0.95     => 18.2  DEGRES )
C                 ( 0.96     => 16.3  DEGRES )
C                 ( 0.97     => 14.1  DEGRES )
C                 ( 0.98     => 11.5  DEGRES )
C                 ( 0.99     =>  8.11 DEGRES )
C                 ( 0.9962   =>  5    DEGRES )
C                 ( 0.99756  =>  4    DEGRES )
C                 ( 0.99863  =>  3    DEGRES )
C                 ( 0.999    =>  2.56 DEGRES )
C                 ( 0.9999   =>  0.8  DEGRES )

               IF( NBAR2FCOL .GE. MXAR2FCOL ) THEN
                  IF( MXAR2FCOL .EQ. 0 ) THEN
C                    LE TABLEAU AR2FCOL EST CREE
                     MXAR2FCOL = NBSOMM
                     CALL TNMCDC( 'ENTIER', MXAR2FCOL, MNAR2FCOL )
                  ELSE
C                    LE TABLEAU AR2FCOL EST RALLONGE
                     CALL TNMCAU( 'ENTIER', MXAR2FCOL, MXAR2FCOL*2,
     %                             MXAR2FCOL, MNAR2FCOL )
                     MXAR2FCOL = MXAR2FCOL * 2
                  ENDIF
                  IF( MNAR2FCOL .LE. 0 ) THEN
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

C              AJOUT de l'ARETE NA du TABLEAU LARETE
               MCN( MNAR2FCOL + NBAR2FCOL ) = NA
               NBAR2FCOL = NBAR2FCOL + 1

               PRINT*
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'lar2f: Les 2 TRIANGLES',NF1,NF2,
     %                   ' SONT COLLES'
               ELSE
                  PRINT*,'lar2f: 2 TRIANGLES',NF1,NF2,' ARE GLUED'
               ENDIF

               CALL QUATRI( NOSOEF(1,NF1), XYZSOM, QF1 )
               PRINT*,'lar2f: TRIANGLE',NF1,' St:',
     %                (NOSOEF(K,NF1),K=1,3),
     %                ' QUALIT=',QF1

               CALL QUATRI( NOSOEF(1,NF2), XYZSOM, QF2 )
               PRINT*,'lar2f: TRIANGLE',NF2,' St:',
     %                (NOSOEF(K,NF2),K=1,3),
     %                ' QUALIT=',QF2
               GOTO 100

            ENDIF

C           LES 2 FACES PRESENTENT ELLES DES SURFACES DISPROPORTIONNEES?
C           COMPARAISON DES SURFACES des TRIANGLES NF1 et NF2
C           ADJACENTS PAR L'ARETE NA du TABLEAU LARETE
C           ------------------------------------------------------------
 25         NF1 = ABS( LARETE( 4, NA ) )
            SF1 = SURTRR( XYZSOM( 1, NOSOEF(1,NF1) ),
     %                    XYZSOM( 1, NOSOEF(2,NF1) ),
     %                    XYZSOM( 1, NOSOEF(3,NF1) ) )
            N = 4
 26         N = N + 1
            NF2 = ABS( LARETE( N, NA ) )
            IF( NF2 .EQ. 0 ) GOTO 26
            SF2 = SURTRR( XYZSOM( 1, NOSOEF(1,NF2) ),
     %                    XYZSOM( 1, NOSOEF(2,NF2) ),
     %                    XYZSOM( 1, NOSOEF(3,NF2) ) )

            IF( SF1 .LE. SF2 ) THEN
               RASF2T = SF1 / SF2
            ELSE
               RASF2T = SF2 / SF1
            ENDIF

            IF( RASF2T .LT. 0.0667 ) THEN

C              RAPPORT des SURFACES DES 2 TRIANGLES TROP PETIT?
C              -----------------------------------------------
               IF( NBA2SF .GE. MXA2SF ) THEN
                  IF( MXA2SF .EQ. 0 ) THEN
C                    LE TABLEAU A2SF EST CREE
                     MXA2SF = NBSOMM
                     CALL TNMCDC( 'ENTIER', MXA2SF, MNA2SF )
                  ELSE
C                    LE TABLEAU A2SF EST RALLONGE
                     CALL TNMCAU( 'ENTIER',MXA2SF,MXA2SF*2,MXA2SF,
     %                             MNA2SF )
                     MXA2SF = MXA2SF * 2
                  ENDIF
                  IF( MNA2SF .LE. 0 ) THEN
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

C              AJOUT de l'ARETE NA du TABLEAU LARETE
               MCN( MNA2SF + NBA2SF ) = NA
               NBA2SF = NBA2SF + 1

               PRINT*
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'lar2sf: Le RAPPORT',RASF2T,
     %                   '<0.0667 des SURFACES des TRIANGLES ADJACENTS',
     %                    NF1,NF2,' EST TROP FAIBLE'
               ELSE
                  PRINT*,'lar2sf: The SURFACE RATIO',RASF2T,
     %                   '<0.0667 of ADJACENT TRIANGLES',NF1,NF2,
     %                   ' is TOO SMALL'
               ENDIF

               CALL QUATRI( NOSOEF(1,NF1), XYZSOM, QF1 )
               PRINT*,'lar2sf: TRIANGLE',NF1,' St:',
     %                   (NOSOEF(K,NF1),K=1,3),
     %                   ' QUALIT=',QF1,' SURFACE=',SF1

               CALL QUATRI( NOSOEF(1,NF2), XYZSOM, QF2 )
               PRINT*,'lar2sf: TRIANGLE',NF2,' St:',
     %                   (NOSOEF(K,NF2),K=1,3),
     %                   ' QUALIT=',QF2,' SURFACE=',SF2
               GOTO 100

            ENDIF

C           LES 2 FACES SONT NON DISPROPORTIONNEES NI COLLEES
C           2 TRIANGLES SONT ADJACENTS A CETTE ARETE
C           -------------------------------------------------
            IF( NBAR2F .GE. MXAR2F ) THEN
               IF( MXAR2F .EQ. 0 ) THEN
C                 LE TABLEAU AR2F EST CREE
                  MXAR2F = 4 * NBSOMM
                  CALL TNMCDC( 'ENTIER', MXAR2F, MNAR2F )
               ELSE
C                 LE TABLEAU AR2F EST RALLONGE
                  CALL TNMCAU( 'ENTIER', MXAR2F, MXAR2F*2,
     %                          MXAR2F, MNAR2F )
                  MXAR2F = MXAR2F * 2
               ENDIF
               IF( MNAR2F .LE. 0 ) THEN
                  IERR = 1
                  GOTO 9999
               ENDIF
            ENDIF

C           AJOUT de l'ARETE NA du TABLEAU LARETE
            MCN( MNAR2F + NBAR2F ) = NA
            NBAR2F = NBAR2F + 1
            GOTO 100

C           L'ARETE NA APPARTIENT A 3 FACES ou PLUS
C           ---------------------------------------
 30         NBAR3F = NBAR3F + 1
            IF( NBAR3F .GT. MXAR3F ) THEN
               IF( MXAR3F .EQ. 0 ) THEN
C                 LE TABLEAU AR3F EST CREE
                  MXAR3F = NBSOMM
                  CALL TNMCDC( 'ENTIER', MXAR3F, MNAR3F )
               ELSE
C                 LE TABLEAU AR3F EST RALLONGE
                  CALL TNMCAU( 'ENTIER', MXAR3F, MXAR3F*2, MXAR3F,
     %                          MNAR3F )
                  MXAR3F = MXAR3F * 2
               ENDIF
               IF( MNAR3F .LE. 0 ) THEN
                  IERR = 1
                  GOTO 9999
               ENDIF
            ENDIF

C           NUMERO LARETE DE L'ARETE SIMPLE NBAR3F
            MCN( MNAR3F-1+NBAR3F ) = NA

         ENDIF

 100  ENDDO

 9999 RETURN
      END
