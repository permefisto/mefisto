      SUBROUTINE EMPITE( NB, NOSOM1, MXSOMM, MXPIL3, IASPI3, NPILE3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EMPILER DANS LA PILE NPILE3 LE NO DES 4 SOMMETS DES
C -----    SOUS-TETRAEDRES
C
C ENTREES:
C --------
C NB     : NOMBRE D INTERVALLES DE SUBDIVISION DE CHAQUE ARETE
C NOSOM1 : NO - 1 DU 1-ER SOMMET A GENERER
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS   DANS LA PILE DES SOMMETS
C MXPIL3 : NOMBRE MAXIMAL DE SOUS-TETRAEDRES DANS LA PILE NPILE3
C
C MODIFIES:
C ---------
C IASPI3 : POINTEUR SUR LE SOMMET DE LA PILE NPILE3
C          -1 EN SORTIE EN CAS DE SATURATION DE LA PILE
C NPILE3 : PILE ( 4 , MXPIL3 ) LA PILE DES 4 NO DES SOMMETS DE
C          SOUS-TETRAEDRES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NPILE3(4,MXPIL3)
      INTEGER           NOS(4)
      INTEGER           NPILOC(4, 8)
      DATA      NPILOC/ 1, 5, 7,  8,
     6                  9, 8, 5,  7,
     6                  9, 5, 6,  7,
     6                  5, 2, 9,  6,
     6                  7, 6, 9, 10,
     6                  6, 3, 7, 10,
     6                  8,10, 9,  7,
     6                  9,10, 8,  4 /
C
10100 FORMAT('ERREUR EMPITE:PILE VALST(4,',I8,') INSUFFISANTE')
C
C     GENERATION DU NO DES SOMMETS DES SOUS-TETRAEDRES
C     ================================================
      IF( NOSOM1 + 2 * (NB**2 +1)  .GT. MXSOMM ) THEN
C        PILE SATUREE
         WRITE(IMPRIM,10100) MXSOMM
         IASPI3 = -1
         RETURN
      ENDIF
C
      K1 = NOSOM1
      IF( NB .EQ. 1 ) THEN
         DO 10 I = 1, 4
            NOS( I ) = K1 + I
 10      CONTINUE
         CALL EMPILN ( IASPI3 , MXPIL3 , NPILE3 , 4 , NOS )
C
      ELSE IF( NB .EQ. 2 ) THEN
C
C        BOUCLE SUR 8 SOUS-TETRAEDRES
C        ----------------------------
         DO 20 J =1, 8
C
            DO 15 I =1, 4
               NOS( I ) = K1 + NPILOC( I, J )
 15         CONTINUE
            CALL EMPILN ( IASPI3 , MXPIL3 , NPILE3 , 4 , NOS )
C
 20      CONTINUE
C
      ELSE
C
C        ERREUR OU CAS RESTANT A PROGRAMMER
         WRITE(IMPRIM,*) 'ERREUR EMPITE: NB INCORRECT ',NB
C
      ENDIF
      END
