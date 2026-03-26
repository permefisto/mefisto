      SUBROUTINE TUERTA1F( NBEFOB, NOSOEF, L1ARET, L2ARET, LARETE,
     %                     NBFSUP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER TOUTE FACE (QT-ANGLES) dont TOUTES LES ARETES
C -----    APPARTIENNENT SEULEMENT A ELLE-MEME I.E. LA FACE N'A PAS DE
C          CONTACT AVEC LES AUTRES FACES. ELLE EST ISOLEE DONC SUPPRIMEE

C ENTREES:
C --------
C NBEFOB : NOMBRE DE FACES (Q-TRIANGLES) DU TABLEAU NOSOEF
C L1ARET : NOMBRE DE MOTS PAR ARET DU TABLEAU LARETE
C L2ARET : NOMBRE MAXIMAL DE FACES DU TABLEAU LARETE
C LARETE : LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          LARETE(4:L1ARET,I)= NO NOSOEF DU FACE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE L1ARET-3 FACES, 
C          LE DERNIER NUMERO DE FACE EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES FACES EST INCOMPLETE

C MODIFIE :
C ---------
C NOSOEF : NUMERO DES 4 SOMMETS DE CHACUN DES NBEFOB FACES
C          UNE FACE DONT TOUTES LES ARETES N'ONT PAS DE FACE OPPOSEE
C          A TOUS SES NUMEROS DE SOMMETS REMIS A ZERO

C SORTIE :
C --------
C NBFSUP : NOMBRE DE FACES SUPPRIMEES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C2345X7..............................................................012
      INTEGER    NOSOEF(4,NBEFOB), LARETE(L1ARET,L2ARET),
     %           NOSOAR(2), NOFA1F(4)

C     RECHERCHE DES ARETES SIMPLES AR1F DE CHAQUE FACE
C     ------------------------------------------------
      NBFSUP = 0
      DO 100 NEF = 1, NBEFOB

         IF( NOSOEF(1,NEF) .EQ. 0 ) GOTO 100

C        NBS NOMBRE DE SOMMETS ou ARETES DE LA FACE
         IF( NOSOEF(4,NEF) .EQ. 0 ) THEN
            NBS = 3
         ELSE
            NBS = 4
         ENDIF

C        NOMBRE D'ARETES APPARTENANT A UNE SEULE FACE
         NBA1F = 0
         DO N = 1, NBS

C           NO DES 2 SOMMETS DE L'ARETE N DE NEF
            IF( N .EQ. NBS ) THEN
               N1 = 1
            ELSE
               N1 = N+1
            ENDIF
            NS1 = NOSOEF(N ,NEF)
            NS2 = NOSOEF(N1,NEF)
            IF( NS1 .LE. NS2 ) THEN
               NOSOAR(1) = NS1
               NOSOAR(2) = NS2
            ELSE
               NOSOAR(1) = NS2
               NOSOAR(2) = NS1
            ENDIF

C           NAR LE NUMERO DE L'ARETE NS1 NS2 DANS LE TABLEAU LARETE
            CALL HACHAR( 2, NOSOAR, L1ARET, L2ARET, LARETE, 3, NAR )
C           NAR =0 SI LE TABLEAU NOSOAR N'A PAS ETE RETROUVE
C               >0 SI LE TABLEAU NOSOAR   A     ETE RETROUVE

            IF( NAR .LE. 0 ) THEN
C              PAS D'ARETE NOSOAR DANS LARETE. EF DEFECTUEUX
               GOTO 100
            ENDIF

C           QUEL EST LE NOMBRE DE FACES CETTE ARETE?
            NBFA = 0
            DO K = 4, L1ARET
               NF = LARETE(K,NAR)
               IF( NF .NE. 0 ) THEN
                  NBFA = NBFA + 1
                  NF1  = NF
               ENDIF
            ENDDO

            IF( NBFA .NE. 1 ) GOTO 100

C           NOMBRE D'ARETES APPARTENANT A UNE SEULE FACE
            NBA1F = NBA1F + 1
C           LE NUMERO DE LA FACE EST NF1
            NOFA1F( NBA1F ) = NF1

         ENDDO

C        BILAN DES ARETES A 1 SEULE FACE DE NEF
         IF( NBA1F .EQ. NBS ) THEN

C           TOUTES LES ARETES APPARTIENNENT ELLES A 1 SEULE et MEME FACE?
            NF1 = NOFA1F( 1 )
            DO K = 2, NBS
               IF( NOFA1F( K ) .NE. NF1 ) GOTO 100
            ENDDO

C           OUI: LA FACE NEF EST DETRUITE
            NBFSUP = NBFSUP + 1
            PRINT*,'tuerta1f: SUPPRESSION de la FACE',NF1,
     %             ' St:',(NOSOEF(K,NEF),K=1,4)
            DO K = 1, 4
               NOSOEF( K, NEF ) = 0
            ENDDO
         ENDIF

 100  ENDDO

      PRINT*,'tuerta1f: SUPPRESSION de',NBFSUP,' FACES ISOLEES'

      RETURN
      END
