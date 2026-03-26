      SUBROUTINE LAR3F( NBSOM0, XYZSOM, L1ARET, L2ARET, LARETE,
     %                  MXAR3F, NBAR3F, MNAR3F, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU TABLEAU LARETE DES ARETES DES FACES D'UNE SURFACE
C -----    CONSTRUIRE LE TABLEAU AR3F DU NUMERO LARETE DES ARETES TRIPLES CAD
C          DES ARETES APPARTENANT A AU MOINS 3 FACES (QUADRANGLE ou TRIANGLE)

C ENTREES:
C --------
C NBSOM0 : NOMBRE INITIAL DE SOMMETS DU MAILLAGE
C L1ARET : NOMBRE DE MOTS PAR ARET DU TABLEAU LARETE
C L2ARET : NOMBRE DE FACES DU TABLEAU LARETE
C LARETE : LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          LARETE(4:3+MXFAAR,I)= NO NOSOFA DE LA FACE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE MXFAAR FACES, 
C          LE NUMERO DE FACE MXFAAR EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES FACES EST INCOMPLETE

C SORTIES:
C --------
C MXAR3F : NOMBRE MAXIMAL D'ARETES DANS UNE FACE RANGES DANS LE TMC AR3F
C NBAR3F : NOMBRE         D'ARETES DANS UNE FACE RANGES DANS LE TMC AR3F
C MNAR3F : >0 ADRESSE MCN DU TMC NAR3F NO LARETE DES ARETES TRIPLES
C          =0 SI NBAR3F=0
C IERR   : =0 PAS D'ERREUR DETECTEE
C          >0 SATURATION DE LA MEMOIRE CENTRALE MCN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY             Avril 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           LARETE(L1ARET,L2ARET)
      REAL              XYZSOM(3,NBSOM0)

      IERR = 0

      IF( MNAR3F .GT. 0 ) CALL TNMCDS( 'ENTIER', MXAR3F, MNAR3F )
      NBAR3F = 0
      MXAR3F = 0

C     RECHERCHE DES ARETES APPARTENANT A PLUS DE 2 FACES
      DO NA = 1, L2ARET
C
         IF( LARETE(1,NA) .NE. 0 ) THEN
C           L'ARETE EST INITIALISEE

C           NOMBRE DE FACES DE CETTE ARETE NA
            NBF = 0
            DO K=4,L1ARET
               IF( LARETE(K,NA) .NE. 0 ) THEN
C                 UNE FACE TQ DE PLUS
                  NBF = NBF + 1
               ENDIF
            ENDDO

C           L'ARETE APPARTIENT APPARTIENT A NBFFACES
            IF( NBF .GE. 3 ) THEN

               NBAR3F = NBAR3F + 1
               IF( NBAR3F .GT. MXAR3F ) THEN
                  IF( MXAR3F .EQ. 0 ) THEN
C                    LE TABLEAU AR3F EST CREE
                     MXAR3F = NBSOM0
                     CALL TNMCDC( 'ENTIER', MXAR3F, MNAR3F )
                  ELSE
C                    LE TABLEAU AR3F EST RALLONGE
                     CALL TNMCAU( 'ENTIER', MXAR3F, MXAR3F*2, MXAR3F,
     %                             MNAR3F )
                     MXAR3F = MXAR3F * 2
                  ENDIF
                  IF( MNAR3F .LE. 0 ) THEN
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

C              NUMERO LARETE DE L'ARETE TRIPLE NBAR3F
               MCN( MNAR3F-1+NBAR3F ) = NA
C
C              AFFICHAGE DE L'ARETE AU MOINS TRIPLE NBAR3F
               NS1 = LARETE(1,NA)
               NS2 = LARETE(2,NA)
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10085) NBAR3F, NS1,(XYZSOM(K,NS1),K=1,3),
     %                                        NS2,(XYZSOM(K,NS2),K=1,3)
               ELSE
                  WRITE(IMPRIM,20085) NBAR3F, NS1,(XYZSOM(K,NS1),K=1,3),
     %                                        NS2,(XYZSOM(K,NS2),K=1,3)
               ENDIF

            ENDIF

         ENDIF

      ENDDO

10085 FORMAT(' lar3f: ARETE AU MOINS TRIPLE',I10,
     %        T41,' : du SOMMET',I10,': ',3G14.7/
     %        T41,'   au SOMMET',I10,': ',3G14.7)
20085 FORMAT(' lar3f: At LEAST TRIPLE EDGE', I10,
     %        T40,' : from VERTEX',I10,': ',3G14.7/
     %        T40,'     to VERTEX',I10,': ',3G14.7)

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'lar3f :',NBAR3F,' ARETES APPARTENANT A 3 FACES'
      ELSE
         PRINT*,'lar3f :',NBAR3F,' EDGES BELONGING to 3 FACES'
      ENDIF

      RETURN
      END
