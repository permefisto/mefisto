      SUBROUTINE LAR2FCOL( NBSOM0,    XYZSOM,    L1ARET, L2ARET, LARETE,
     %                     N1SOEF,    NBEFOB,    NOSOEF,
     %                     MXAR2FCOL, NBAR2FCOL, MNAR2FCOL, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU TABLEAU LARETE DES ARETES CONSTRUCTION DU TABLEAU
C -----    DES COUPLES DE TRIANGLES ADJACENTS COLLES

C ENTREES:
C --------
C NBSOM0 : NOMBRE INITIAL DE SOMMETS DU MAILLAGE
C XYZSOM : XYZ DES NBSOM0 SOMMETS DU MAILLAGE
C L1ARET : NOMBRE DE MOTS PAR ARET DU TABLEAU LARETE
C L2ARET : NOMBRE DE FACES DU TABLEAU LARETE
C LARETE : LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          LARETE(4:3+MXFAAR,I)= NO NOSOFA DE LA FACE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE MXFAAR FACES, 
C          LE NUMERO DE FACE MXFAAR EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES FACES EST INCOMPLETE

C N1SOEF : NOMBRE DE SOMMETS DECLARABLES DANS NOSOEF PAR EF
C NBEFOB : NOMBRE D'EF ACTIFS DANS LE TABLEAU NOSOEF
C NOSOEF : NUMERO XYZSOM DES N1SOEF SOMMETS DES NBEFOB EF

C SORTIES:
C --------
C MXAR2FCOL: NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC AR2FCOL
C NBAR2FCOL: NOMBRE         D'ARETES DE NO LARETE DANS LE TMC AR2FCOL
C MNAR2FCOL: >0 ADRESSE MCN DU TMC NAR2FCOL NO LARETE DE L'ARETE DES
C               COUPLES DE TRIANGLES ADJACENTS COLLES
C            =0 SI NBAR2FCOL=0
C IERR     : =0 PAS D'ERREUR DETECTEE
C            >0 SATURATION DE LA MEMOIRE CENTRALE MCN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY               Mai 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           LARETE(L1ARET,L2ARET), NOSOEF(N1SOEF,NBEFOB)
      REAL              XYZSOM(3,NBSOM0)

      IERR = 0

C     SI LE TABLEAU AR2FCOL EXISTE ALORS IL EST DETRUIT
      IF( MNAR2FCOL .GT. 0 ) CALL TNMCDS('ENTIER', MXAR2FCOL, MNAR2FCOL)
      NBAR2FCOL = 0
      MXAR2FCOL = 0

C     RECHERCHE DES ARETES APPARTENANT A 2 FACES DE RAPPORT FAIBLE
C     DES SURFACES
      DO 100 NA = 1, L2ARET
C
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
            IF( NBF .EQ. 2 ) THEN

C              CALCUL DU COSINUS DE L'ANGLE DIEDRE DES 2 TRIANGLES
C              ADJACENTS PAR L'ARETE NA
C              ---------------------------------------------------
               NF1 = ABS( LARETE( 4, NA ) )
               NS1 = LARETE( 1, NA )
               NS2 = LARETE( 2, NA )

C              RECHERCHE DU 3-EME SOMMET DE NF1
               DO K = 1, N1SOEF
                  NS3 = NOSOEF( K, NF1 )
                  IF( NS3 .NE. NS1 .AND. NS3 .NE. NS2 ) GOTO 10
               ENDDO
               GOTO 100

C              RECHERCHE DU 3-EME SOMMET DE NF2
 10            N = 4
               N = N + 1
               NF2 = ABS( LARETE( N, NA ) )
               IF( NF2 .EQ. 0 ) GOTO 10
               DO K = 1, N1SOEF
                  NS4 = NOSOEF( K, NF2 )
                  IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 20
               ENDDO
               GOTO 100

C              COSINUS DE L'ANGLE DIEDRE DES 2 TRIANGLES ADJACENTS
 20            CALL COS2TR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %                      XYZSOM(1,NS2), XYZSOM(1,NS1), XYZSOM(1,NS4),
     %                      COS2T, IERR1, IERR2 )
               IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) GOTO 100

               IF( COS2T .LT. -0.99 ) THEN

C              COANPL: SEUIL DU COSINUS DE L'ANGLE FORME PAR LES NORMALES AUX
C              2 FACES ET AU DESSUS DUQUEL LES FACES SONT CONSIDEREES COPLANAIRES
C              ( 0.95     => 18.2  DEGRES )
C              ( 0.96     => 16.3  DEGRES )
C              ( 0.97     => 14.1  DEGRES )
C              ( 0.98     => 11.5  DEGRES )
C              ( 0.99     =>  8.11 DEGRES )
C              ( 0.9962   =>  5    DEGRES )
C              ( 0.99756  =>  4    DEGRES )
C              ( 0.99863  =>  3    DEGRES )
C              ( 0.999    =>  2.56 DEGRES )
C              ( 0.9999   =>  0.8  DEGRES )

                  IF( NBAR2FCOL .GE. MXAR2FCOL ) THEN
                     IF( MXAR2FCOL .EQ. 0 ) THEN
C                       LE TABLEAU AR2FCOL EST CREE
                        MXAR2FCOL = NBSOM0
                        CALL TNMCDC( 'ENTIER', MXAR2FCOL, MNAR2FCOL )
                     ELSE
C                       LE TABLEAU AR2FCOL EST RALLONGE
                        CALL TNMCAU( 'ENTIER', MXAR2FCOL, MXAR2FCOL*2,
     %                                MXAR2FCOL, MNAR2FCOL )
                        MXAR2FCOL = MXAR2FCOL * 2
                     ENDIF
                     IF( MNAR2FCOL .LE. 0 ) THEN
                        IERR = 1
                        GOTO 9999
                     ENDIF
                  ENDIF

C                 AJOUT de l'ARETE NA du TABLEAU LARETE
                  MCN( MNAR2FCOL + NBAR2FCOL ) = NA
                  NBAR2FCOL = NBAR2FCOL + 1

                  PRINT*
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT*,'lar2fcol: Les 2 TRIANGLES',NF1,NF2,
     %                      ' SONT COLLES'
                  ELSE
                     PRINT*,'lar2fcol: 2 TRIANGLES',NF1,NF2,' ARE GLUED'
                  ENDIF

                  CALL QUATRI( NOSOEF(1,NF1), XYZSOM, QF1 )
                  PRINT*,'lar2fcol: TRIANGLE',NF1,' St:',
     %                   (NOSOEF(K,NF1),K=1,3),
     %                   ' QUALIT=',QF1

                  CALL QUATRI( NOSOEF(1,NF2), XYZSOM, QF2 )
                  PRINT*,'lar2fcol: TRIANGLE',NF2,' St:',
     %                   (NOSOEF(K,NF2),K=1,3),
     %                   ' QUALIT=',QF2

               ENDIF

            ENDIF

         ENDIF

 100  ENDDO

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'lar2fcol:',NBAR2FCOL,' COUPLES de TRIANGLES COLLES'
      ELSE
         PRINT*,'lar2fcol:',NBAR2FCOL,' COUPLES of GLUED TRIANGLES'
      ENDIF

      RETURN
      END
