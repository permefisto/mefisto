      SUBROUTINE TETR1ST( NBEF, NUSOTE, NOST, MXTE1S, NBTE1S, NOTE1S )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LISTER DANS NOTE1S LE NUMERO NOTETR DES NBTE1S TETRAEDRES
C -----    DE SOMMET NOST    (TROP BESTIAL... A AMELIORER)

C ENTREES:
C --------
C NOST   : NUMERO DU SOMMET COMMUN A TOUS LES TETRAEDRES
C NUSOTE : LISTE DU NO DES SOMMETS DES NBEF TETRAEDRES ( POUR UN OBJET )
C          NUSOTE(NTE,NST) NO XYZPOI DU SOMMET NST (1a4) DU TETRAEDRE NTE
C MXTE1S : NOMBRE D'ENTIERS DU TABLEAU NOTE1S

C SORTIES:
C --------
C NBTE1S : NOMBRE DE TETRAEDRES DE SOMMET NOST SI IERR=0
C NOTE1S : NOTE1S(I) >0 NUMERO NOTETR DU TETRAEDRE I DE SOMMET NOST
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     AUTEUR : ALAIN PERRONNET Saint Pierre du Perray       Octobre 2020
C2345X7..............................................................012
      INTEGER  NBEF, NUSOTE(NBEF,4), NOTE1S(MXTE1S)

      NBTE1S = 0

      DO 10 NTE = 1, NBEF

         DO K = 1, 4

            IF( NUSOTE( NTE, K ) .EQ. NOST ) THEN

               IF( NBTE1S .GE. MXTE1S ) THEN
              PRINT*,'tetr1st: TABLEAU NOTE1S SATURE. AUGMENTER MXTE1S='
     %              ,MXTE1S,' TABLEAU INCOMPLET en SORTIE'
              GOTO 9999
               ENDIF

C              AJOUT DANS NOTE1S DU TETRAEDRE NTE DE SOMMET NOST
               NBTE1S = NBTE1S + 1
               NOTE1S( NBTE1S ) = NTE
               GOTO 10

            ENDIF

         ENDDO

 10   ENDDO

 9999 RETURN
      END

