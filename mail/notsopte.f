      SUBROUTINE NOTSOPTE( NT,  NO1TSO, NOTESO, NBSOTE, NSTETR,
     %                     NFR, NTOP,   NSOP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LE NUMERO DES TETRAEDRES ET SOMMETS OPPOSES
C -----   AUX 4 FACES DU TETRAEDRE NT
C
C ENTREES:
C --------
C NT     : NUMERO DANS NSTETR DU TETRAEDRE A TRAITER
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DE CHAQUE SOMMET
C NOTESO : NUMERO DU TETRAEDRE ET SUIVANT DES SOMMETS
C NBSOTE : NOMBRE DE SOMMETS DECLARABLES PAR TETRAEDRE DANS NSTETR(>3)
C NSTETR : LISTE DES 4 SOMMETS DES TETRAEDRES
C
C SORTIES:
C --------
C NFR    : NOMBRE DE TETRAEDRES OPPOSES RETROUVES ( 0 A 4 )
C NTOP   : NUMERO DU TETRAEDRE OPPOSE A CHACUNE DES 4 FACES 123 234 341 412
C          0 SI PAS DE TETRAEDRE OPPOSE
C NSOP   : NUMERO DU SOMMET OPPOSE A CHACUNE DES 4 FACES
C          0 SI PAS DE TETRAEDRE OPPOSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1992
C2345X...............................................................012
      INTEGER           NSTETR(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*),
     %                  NTOP(4),NSOP(4),
     %                  NSF(3)
C
C     RECHERCHE D'UN TETRAEDRE ADJACENT PAR UNE FACE
C     ET CONTENANT LE 4 EME SOMMET DE NT
C     ----------------------------------------------
      NFR = 0
      DO 6 I=1,4
C
C        MISE A ZERO DE NTOP ET NSOP
         NTOP(I) = 0
         NSOP(I) = 0
C
C        LE NUMERO DES 3 SOMMETS DE LA FACE
         NSF(1) = NSTETR(I,NT)

         IF( I .EQ. 4 ) THEN
            NDT = 1
         ELSE
            NDT = I + 1
         ENDIF
         NSF(2) = NSTETR(NDT,NT)

         IF( NDT .EQ. 4 ) THEN
            NDT = 1
         ELSE
            NDT = NDT + 1
         ENDIF
         NSF(3) = NSTETR(NDT,NT)
C
C        TRI CROISSANT DES SOMMETS DE LA FACE
         CALL TRI3NO( NSF, NSF )
C
C        POSITION DANS NOTESO DU 1-ER TETRAEDRE DE SOMMET NS
         NDT = NO1TSO( NSTETR(I,NT) )
C
C        TANT QU'IL EXISTE UN TETRAEDRE DE SOMMET NS FAIRE
 3       IF( NDT .GT. 0 ) THEN

C           LE NUMERO DU TETRAEDRE DANS NSTETR
            NT1 = NOTESO(1,NDT)
            IF( NT .EQ. NT1 ) GOTO 5

C           LE NUMERO DE FACE NF DANS NT1 DE CETTE FACE NSF (NO CROISSANTS)
            CALL NO1F1T( NSF, NSTETR(1,NT1), NF )

            IF( NF .GT. 0 ) THEN
C              FACE RETROUVEE
               NFR = NFR + 1
C              LE TETRAEDRE OPPOSE A LA FACE
               NTOP( I ) = NT1

C              LE 4-EME SOMMET NS5 OPPOSE A LA FACE NSF
               IF( NF .EQ. 1 ) THEN
                  NS5 = 4
               ELSE
                  NS5 = NF - 1
               ENDIF
               NSOP( I ) = NSTETR( NS5, NT1 )
               GOTO 6

            ENDIF

C           LE TETRAEDRE SUIVANT
 5          NDT = NOTESO(2,NDT)
            GOTO 3

         ENDIF

 6    ENDDO

      RETURN
      END
