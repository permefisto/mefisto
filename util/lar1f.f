      SUBROUTINE LAR1F( NBSOM0, XYZSOM, L1ARET, L2ARET, LARETE,
     %                  MXAR1F, NBAR1F, MNAR1F, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU TABLEAU LARETE DES ARETES DES FACES D'UNE SURFACE
C -----    CONSTRUIRE LE TABLEAU AR1F DU NUMERO LARETE DES ARETES SIMPLES
C          I.E. DES ARETES APPARTENANT A UNE SEULE FACE (QUADRANGLE ou TRIANGLE)

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
C MXAR1F : NOMBRE MAXIMAL D'ARETES DANS UNE FACE RANGES DANS LE TMC AR1F
C NBAR1F : NOMBRE         D'ARETES DANS UNE FACE RANGES DANS LE TMC AR1F
C MNAR1F : >0 ADRESSE MCN DU TMC NAR1F NO LARETE DES ARETES SIMPLES
C          =0 SI NBAR1F=0
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

      IERR   = 0

      IF( MNAR1F .GT. 0 ) CALL TNMCDS( 'ENTIER', MXAR1F, MNAR1F )
      NBAR1F = 0
      MXAR1F = 0

C     RECHERCHE DES ARETES APPARTENANT A UNE SEULE FACE
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
            IF( NBF .EQ. 1 ) THEN

               NBAR1F = NBAR1F + 1
               IF( NBAR1F .GT. MXAR1F ) THEN
                  IF( MXAR1F .EQ. 0 ) THEN
C                    LE TABLEAU AR1F EST CREE
                     MXAR1F = NBSOM0
                     CALL TNMCDC( 'ENTIER', MXAR1F, MNAR1F )
                  ELSE
C                    LE TABLEAU AR1F EST RALLONGE
                     CALL TNMCAU( 'ENTIER', MXAR1F, MXAR1F*2, MXAR1F,
     %                             MNAR1F )
                     MXAR1F = MXAR1F * 2
                  ENDIF
                  IF( MNAR1F .LE. 0 ) THEN
                     IERR = 1
                     GOTO 9999
                  ENDIF
               ENDIF

C              NUMERO LARETE DE L'ARETE SIMPLE NBAR1F
               MCN( MNAR1F-1+NBAR1F ) = NA

C              AFFICHAGE DE L'ARETE SIMPLE NBAR1F
C              test bidon pour eviter l'affichage
               if( ns1.eq.0) then
               NS1 = LARETE(1,NA)
               NS2 = LARETE(2,NA)
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10085) LARETE(4,NA),
     %                  NS1,(XYZSOM(K,NS1),K=1,3),
     %                  NS2,(XYZSOM(K,NS2),K=1,3)
               ELSE
                  WRITE(IMPRIM,20085) LARETE(4,NA),
     %                  NS1,(XYZSOM(K,NS1),K=1,3),
     %                  NS2,(XYZSOM(K,NS2),K=1,3)
               ENDIF
               endif

            ENDIF

         ENDIF

      ENDDO

10085 FORMAT(' lar1f: ARETE d''UNE SEULE FACE',I10,
     %        T43,' : du SOMMET',I10,': ',3G14.7/
     %        T43,'   au SOMMET',I10,': ',3G14.7)
20085 FORMAT(' lar1f: ONLY ONE FACE EDGE', I10,
     %        T38,' : from VERTEX',I10,': ',3G14.7/
     %        T38,'     to VERTEX',I10,': ',3G14.7)

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'lar1f :',NBAR1F,' ARETES APPARTENANT A UNE SEULE FACE'
      ELSE
         PRINT*,'lar1f :',NBAR1F,' EDGES BELONGING to ONE FACE'
      ENDIF

      RETURN
      END
