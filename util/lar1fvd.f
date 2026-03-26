      SUBROUTINE LAR1FVD( XYZSOM, NBSA1F, NUSA1F,  NUNVST, DISVOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU TABLEAU NUSA1F DES SOMMETS NUSA1F DES ARETES
C -----    SIMPLES D'UNE SURFACE, CONSTRUIRE LA TABLEAU DES SOMMETS
C          DES ARETES SIMPLES LES PLUS PROCHES ET LEUR DISTANCE

C ENTREES:
C --------
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE
C NBAR1F : NOMBRE         D'ARETES DANS UNE FACE RANGES DANS LE TMC AR1F
C NBSA1F : NOMBRE         DE SOMMETS DES ARETES SIMPLES DANS SA1F
C NUSA1F : NO XYZSOM DES SOMMETS DES ARETES SIMPLES

C SORTIES:
C --------
C NUNVST : INDICE DANS NUSA1F DU SOMMET VOISIN
C DISVOI :  DISTANCE AU SOMMET VOISIN LE PLUS PROCHE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2015
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER            NUSA1F(NBSA1F), NUNVST(NBSA1F)
      REAL               XYZSOM(3,*), DISVOI(NBSA1F)

C     MISE A ZERO DE NUNVST
      CALL AZEROI( NBSA1F, NUNVST )

      DO 100 J=1,NBSA1F

         IMIN = 0
         DMIN = 1E28

C        NUMERO XYZSOM DU SOMMET D'AU MOINS UNE ARETE SIMPLE
         NS1 = NUSA1F( J )

         DO I=1,NBSA1F
            IF( I .NE. J ) THEN
C              CALCUL DU CARRE DE LA DISTANCE NS1-NS2
               NS2  = NUSA1F( I )
               D    = ( XYZSOM(1,NS2) - XYZSOM(1,NS1) ) ** 2
     %              + ( XYZSOM(2,NS2) - XYZSOM(2,NS1) ) ** 2
     %              + ( XYZSOM(3,NS2) - XYZSOM(3,NS1) ) ** 2
               IF( D .LT. DMIN ) THEN
                  IMIN = I
                  DMIN = D
               ENDIF
            ENDIF
         ENDDO

C        IMIN DANS NUSA1F EST LE SOMMET LE PLUS PROCHE DU SOMMET J
         IF( IMIN .GT. 0 ) THEN

C           NUMERO XYZSOM DU SOMMET LE PLUS PROCHE DE NS1
            NS2 = NUSA1F( IMIN )

C           NUMERO NUSA1F DU SOMMET NS2 LE PLUS PROCHE DE NS1
            NUNVST( J ) = IMIN

C           DISTANCE ENTRE CES 2 SOMMETS DE SA1F
            DMIN = SQRT( DMIN )
            DISVOI( IMIN ) = DMIN
            DISVOI( J    ) = DMIN

            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10100) DMIN, NS1, NS2
               WRITE(IMPRIM,10101) NS1,(XYZSOM(K,NS1),K=1,3)
               WRITE(IMPRIM,10101) NS2,(XYZSOM(K,NS2),K=1,3)
            ELSE
               WRITE(IMPRIM,20100) DMIN, NS1, NS2
               WRITE(IMPRIM,20101) NS1,(XYZSOM(K,NS1),K=1,3)
               WRITE(IMPRIM,20101) NS2,(XYZSOM(K,NS2),K=1,3)
            ENDIF

C           CES 2 SOMMETS SONT IL IDENTIFIABLES?
            CALL XYZIDE( XYZSOM(1,NS1), XYZSOM(1,NS2), IDENTQ )
            IF( IDENTQ .EQ. 1 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10103) NS1, NS2
               ELSE
                  WRITE(IMPRIM,20103) NS1, NS2
               ENDIF
            ENDIF

         ENDIF

 100  CONTINUE

10100 FORMAT(/' lar1fvd: DISTANCE',G14.6,' ENTRE le SOMMET',I10,
     %        ' et SON PLUS PROCHE SOMMET',I10,
     %        ' des ARETES SIMPLES')
20100 FORMAT(/' lar1fvd: DISTANCE',G14.6,' BETWEEN the VERTEX',I10,
     %        ' and its NEAREST VERTEX',I10,
     %        ' of SINGLE EDGES')

10101 FORMAT(' SOMMET',I10,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
20101 FORMAT(' VERTEX',I10,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)

10103 FORMAT(' lar1fvd: Le SOMMET',I10,' EST IDENTIQUE a SON PLUS PROCHE
     % SOMMET',I10,' et DEVRAIENT ETRE IDENTIFIES' )
20103 FORMAT(' lar1fvd: the VERTEX ',I10,' IS IDENTICAL to its NEAREST V
     %ERTEX',I10,' and SHOULD BE IDENTIFIED' )

      RETURN
      END
