      SUBROUTINE TRIAMAXQ( NBTR1, NBTR2, NUSTRA, XYZST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECHANGER LES DIAGONALES DE 2 TRIANGLES ADJACENTS PAR UNE ARETE
C -----    SI LA QUALITE MINIMALE DU COUPLE DE TRIANGLES EST AUGMENTEE
C
C ENTREES :
C ---------
C NBTR1  : NUMERO DU PREMIER TRIANGLE A AMELIORER LA QUALITE
C NBTR2  : NUMERO DU DERNIER TRIANGLE A AMELIORER LA QUALITE
C XYZST  : 3 XYZ  DES SOMMETS DES TRIANGLES
C
C MODIFIES:
C ---------
C NUSTRA : (3,NBTR1:NBTR2) NO DANS XYZST DES 3 SOMMETS DES NBTRA TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray  Octobre 2011
C2345X7..............................................................012
      DOUBLE PRECISION  XYZST(3,*)
      INTEGER           NUSTRA(3,*), NOSOTR(3)
C
      DOUBLE PRECISION  SURTRD, QUALIT1, QUALIT2,
     %                  AIRE1, AIRE2, AIRE3, AIRE4
C
C     NOMBRE D'ECHANGES DE DIAGONALE
 10   NBECHA = 0
C
      IF( NBTR1 .EQ. NBTR2 ) RETURN
C
      DO 100 NT1 = NBTR1, NBTR2
C
C        BOUCLE SUR LES 3 ARETES DU TRIANGLE NT
         DO 90 K1 = 1, 3
C
C           CALCUL DE LA QUALITE DU TRIANGLE NT1
            CALL QUATRID( NUSTRA(1,NT1), XYZST, QUALIT1 )
            IF( K1 .LT. 3 ) THEN
               K2 = K1+1
            ELSE
               K2 = 1
            ENDIF
            IF( K1 .GT. 1 ) THEN
               K0 = K1-1
            ELSE
               K0 = 3
            ENDIF
C
C           LE NO DES 2 SOMMETS DE L'ARETE K1 DU TRIANGLE NT1
            NS1 = NUSTRA(K1,NT1)
            NS2 = NUSTRA(K2,NT1)
C
C           ET DU SOMMET OPPOSE
            NS0 = NUSTRA(K0,NT1)
C
C           CALCUL DE L'AIRE1 DU TRIANGLE NT1
            AIRE1 = SURTRD( XYZST(1,NS1), XYZST(1,NS2), XYZST(1,NS0) )
C
C           EXISTE T IL UN TRIANGLE OPPOSE A CETTE ARETE K1 DE NT1?
            DO NT2 = NBTR1, NBTR2
C
               IF( NT2 .NE. NT1 ) THEN
C
                  DO M1 = 1, 3
                     IF( M1 .LT. 3 ) THEN
                        M2 = M1+1
                     ELSE
                        M2 = 1
                     ENDIF
C
C                    LE NO DES 2 SOMMETS DE L'ARETE M1 DU TRIANGLE NT2
                     NS3 = NUSTRA(M1,NT2)
                     NS4 = NUSTRA(M2,NT2)
C
C                    CETTE ARETE M1 DE NT2 EST ELLE L'ARETE K1 DE NT1?
                     IF( (NS1 .EQ. NS4 .AND. NS2 .EQ. NS3) .OR.
     %                   (NS2 .EQ. NS4 .AND. NS1 .EQ. NS3) ) THEN
C                       LES TRIANGLES NT1 ET NT2 SONT ADJACENTS PAR CETTE ARETE
C                       CALCUL DE LA QUALITE DU TRIANGLE NT2
                        CALL QUATRID( NUSTRA(1,NT2), XYZST, QUALIT2 )
C
C                       LA QUALITE MINIMALE DES 2 TRIANGLES NT1 NT2
                        QUALIT1 = MIN( QUALIT1, QUALIT2 ) * 1.0001D0
C                       AUGMENTEE POUR EVITER UN ECHANGE DE DIAGONALES CYCLIQUE
C                       ECHANGE SEULEMENT SI LA QUALITE S'ELEVE VRAIMENT
C                       ET PAS EN CAS D'EGALITE DES QUALITES
C
C                       LA QUALITE DES 2 TRIANGLES DE L'AUTRE SUBDIVISION
C                       DU QUADRANGLE NT1-NT2
                        IF( M1 .GT. 1 ) THEN
                           M0 = M1-1
                        ELSE
                           M0 = 3
                        ENDIF
                        NS5 = NUSTRA(M0,NT2)
C
C                       LE TRIANGLE NS0-NS5-NS2
                        NOSOTR(1) = NS0
                        NOSOTR(2) = NS5
                        NOSOTR(3) = NS2
C                       CALCUL DE LA QUALITE DU TRIANGLE NS0-NS5-NS2
                        CALL QUATRID( NOSOTR, XYZST, QUALIT2 )
                        IF( QUALIT2 .LE. QUALIT1 ) GOTO 90
C
C                       LE TRIANGLE NS0-NS1-NS5
                        NOSOTR(1) = NS0
                        NOSOTR(2) = NS1
                        NOSOTR(3) = NS5
C                       CALCUL DE LA QUALITE DU TRIANGLE NS0-NS1-NS5
                        CALL QUATRID( NOSOTR, XYZST, QUALIT2 )
                        IF( QUALIT2 .LE. QUALIT1 ) GOTO 90
C
C                       CALCUL DE L'AIRE2 DU TRIANGLE NT2
                        AIRE2 = SURTRD( XYZST(1,NS3), XYZST(1,NS4),
     %                                  XYZST(1,NS5) )
C
C                       CALCUL DE L'AIRE3 DU TRIANGLE NS0-NS5-NS2
                        AIRE3 = SURTRD( XYZST(1,NS0), XYZST(1,NS5),
     %                                  XYZST(1,NS2) )
C
C                       CALCUL DE L'AIRE4 DU TRIANGLE NS0-NS1-NS5
                        AIRE4 = SURTRD( XYZST(1,NS0), XYZST(1,NS1),
     %                                  XYZST(1,NS5) )
C
C                       L'UNION DES 2 TRIANGLES EST ELLE UN QUADRANGLE CONVEXE?
                        AIRE2 = AIRE1 + AIRE2
                        AIRE3 = AIRE3 + AIRE4
                        IF( ABS( AIRE2 - AIRE3 ) .LE. AIRE2*1D-8 ) THEN
C
c                          QUADRANGLE CONVEXE
C                          LES TRIANGLES NS0-NS5-NS2 NS0-NS1-NS5
C                          ONT UNE MEILLEURE MINIMALE QUALITE
C                          => ECHANGE DE NT1+NT2 EN NS0-NS5-NS2 + NS0-NS1-NS5
                           NUSTRA(1,NT1)=NS0
                           NUSTRA(2,NT1)=NS5
                           NUSTRA(3,NT1)=NS2
C
                           NUSTRA(1,NT2)=NS0
                           NUSTRA(2,NT2)=NS1
                           NUSTRA(3,NT2)=NS5
C
                           NBECHA = NBECHA + 1
C
C                          PASSAGE AU TRIANGLE NT1 SUIVANT
                           GOTO 100
C
                        ENDIF
C
                     ENDIF
C
                  ENDDO
C
               ENDIF
C
            ENDDO
C
 90      CONTINUE
C
 100  CONTINUE
C
C     ITERATION JUSQU'A NE PLUS AVOIR D'ECHANGE DES DIAGONALES
      IF( NBECHA .GT. 0 ) GOTO 10
C
      RETURN
      END
