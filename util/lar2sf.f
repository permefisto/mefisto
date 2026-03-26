      SUBROUTINE LAR2SF( NBSOM0, XYZSOM, L1ARET, L2ARET, LARETE,
     %                   N1SOEF, NBEFOB, NOSOEF,
     %                   MXA2SF, NBA2SF, MNA2SF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU TABLEAU LARETE DES ARETES DES FACES D'UNE SURFACE
C -----    CONSTRUIRE LE TABLEAU A2SF DU NUMERO LARETE DES ARETES 
C          DONT LE RAPPORT DES SURFACES DES 2 TRIANGLES ADJACENTS
C          EST <0.0667 TROP FAIBLE

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

C N1SOEF : NOMBRE DE SOMMETS DECLARBLES DANS NOSOEF PAR EF
C NBEFOB : NOMBRE D'EF ACTIFS DANS LE TABLEAU NOSOEF
C NOSOEF : NUMERO XYZSOM DES N1SOEF SOMMETS DES NBEFOB EF

C SORTIES:
C --------
C MXA2SF : NOMBRE MAXIMAL D'ARETES DECLARABLES  DANS LE TMC A2SF
C NBA2SF : NOMBRE         D'ARETES DE NO LARETE DANS LE TMC A2SF
C MNA2SF : >0 ADRESSE MCN DU TMC NA2SF NO LARETE DES ARETES DE RATIO FAIBLE
C          =0 SI NBA2SF=0
C IERR   : =0 PAS D'ERREUR DETECTEE
C          >0 SATURATION DE LA MEMOIRE CENTRALE MCN
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

C     SI LE TABLEAU A2SF EXISTE ALORS IL EST DETRUIT
      IF( MNA2SF .GT. 0 ) CALL TNMCDS( 'ENTIER', MXA2SF, MNA2SF )
      NBA2SF = 0
      MXA2SF = 0

C     RECHERCHE DES ARETES APPARTENANT A 2 FACES DE RAPPORT FAIBLE
C     DES SURFACES
      DO 100 NA = 1, L2ARET
C
         IF( LARETE(1,NA) .NE. 0 ) THEN
C           L'ARETE EST INITIALISEE

C           NOMBRE DE FACES DE CETTE ARETE NA
            NBF = 0
            DO K=4,L1ARET
               NF1 = ABS( LARETE(K,NA) )
               IF( NF1 .GT. 0 ) THEN
C                 UNE FACE DE PLUS
                  NBF = NBF + 1
               ENDIF
            ENDDO

C           L'ARETE APPARTIENT APPARTIENT A NBF FACES
            IF( NBF .EQ. 2 ) THEN

C              COMPARAISON DES SURFACES des TRIANGLES NF1 et NF2
C              ADJACENTS PAR L'ARETE NA du TABLEAU LARETE
C              --------------------------------------------------
               NF1 = ABS( LARETE( 4, NA ) )
               SF1 = SURTRR( XYZSOM( 1, NOSOEF(1,NF1) ),
     %                       XYZSOM( 1, NOSOEF(2,NF1) ),
     %                       XYZSOM( 1, NOSOEF(3,NF1) ) )
               N = 4
 10            N = N + 1
               NF2 = ABS( LARETE( N, NA ) )
               IF( NF2 .EQ. 0 ) GOTO 10
               SF2 = SURTRR( XYZSOM( 1, NOSOEF(1,NF2) ),
     %                       XYZSOM( 1, NOSOEF(2,NF2) ),
     %                       XYZSOM( 1, NOSOEF(3,NF2) ) )

               IF( SF1 .LE. SF2 ) THEN
                  RASF2T = SF1 / SF2
               ELSE
                  RASF2T = SF2 / SF1
               ENDIF

               IF( RASF2T .LT. 0.0667 ) THEN

C                 RAPPORT des 2 SURFACES TROP PETIT: A SIGNALER
                  IF( NBA2SF .GE. MXA2SF ) THEN
                     IF( MXA2SF .EQ. 0 ) THEN
C                       LE TABLEAU A2SF EST CREE
                        MXA2SF = NBSOM0
                        CALL TNMCDC( 'ENTIER', MXA2SF, MNA2SF )
                     ELSE
C                       LE TABLEAU A2SF EST RALLONGE
                        CALL TNMCAU( 'ENTIER',MXA2SF,MXA2SF*2,MXA2SF,
     %                                MNA2SF )
                        MXA2SF = MXA2SF * 2
                     ENDIF
                     IF( MNA2SF .LE. 0 ) THEN
                        IERR = 1
                        GOTO 9999
                     ENDIF
                  ENDIF

C                 AJOUT de l'ARETE NA du TABLEAU LARETE
                  MCN( MNA2SF + NBA2SF ) = NA
                  NBA2SF = NBA2SF + 1

                  PRINT*
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT*,'lar2sf: Le RAPPORT',RASF2T,
     %                 '<0.0667 des SURFACES des TRIANGLES ADJACENTS',
     %                  NF1,NF2,' EST TROP FAIBLE'
                  ELSE
                     PRINT*,'lar2sf: The SURFACE RATIO',RASF2T,
     %                      '<0.0667 of ADJACENT TRIANGLES',NF1,NF2,
     %                      ' is TOO SMALL'
                  ENDIF

                  CALL QUATRI( NOSOEF(1,NF1), XYZSOM, QF1 )
                  PRINT*,'lar2sf: TRIANGLE',NF1,' St:',
     %                   (NOSOEF(K,NF1),K=1,3),
     %                   ' QUALIT=',QF1,' SURFACE=',SF1

                  CALL QUATRI( NOSOEF(1,NF2), XYZSOM, QF2 )
                  PRINT*,'lar2sf: TRIANGLE',NF2,' St:',
     %                   (NOSOEF(K,NF2),K=1,3),
     %                   ' QUALIT=',QF2,' SURFACE=',SF2

               ENDIF

            ENDIF

         ENDIF

 100  ENDDO

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'lar2sf:',NBA2SF,' ARETES de RAPPORT<0.0667 des SURFACES
     % des 2 TRIANGLES ADJACENTS'
      ELSE
         PRINT*,'lar2sf:',NBA2SF,' EDGES with a RATIO<0.0667 of 2 ADJACE
     %NT TRIANGLE SURFACES'
      ENDIF

      RETURN
      END
