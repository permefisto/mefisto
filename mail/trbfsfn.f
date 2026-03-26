      SUBROUTINE TRBFSFN( NBFACE, MOFACE, MXFACE, LFACES, XYZPOI,
     %                    NOFACE, DIFACE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE BARYCENTRE DE CHAQUE FACE DES SURFACES NOMMEES
C -----    CALCULER LA COTE AXONOMETRIQUE DE CES BARYCENTRES
C
C ENTREES:
C --------
C NBFACE : NOMBRE DE FACE A TRIER
C MOFACE : NOMBRE DE MOTS PAR FACE DU TABLEAU LFACES
C          =5 SI INTERP=1, 9 SI INTERP>=2
C MXFACE : NOMBRE DE FACES DU TABLEAU NFACE
C LFACES : TABLEAU NUMERO DES FACES DES SURFACES NOMMEES
C          LFACES(1,I)= NO DE LA SURFACE NOMMEE DE LA FACE I
C          LFACES(2...MOFACE,I)= NO DU POINT DE LA FACE I
C          LFACES(MOFACE,I)= 0 SI QUADRANGLE
C XYZPOI : 3 COORDONNEES DES POINTS
C
C SORTIES:
C --------
C NOFACE : NUMERO DE LA FACE DANS LFACES POUR LE TRI
C DIFACE : COTE AXONOMETRIQUE DU BARYCENTRE DES FACES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY SEPTEMBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      INTEGER           LFACES(1:MOFACE,1:MXFACE)
      REAL              XYZPOI(1:3,1:*)
      INTEGER           NOFACE(1:NBFACE)
      REAL              DIFACE(1:NBFACE),
     %                  XYZ(3)
C
C     BOUCLE SUR LES FACES A TRIER
      DO N = 1, NBFACE
         NOFACE(N) = N
C
C        LES COORDONNEES DU BARYCENTRE DE LA FACE
         XYZ(1) = 0.0
         XYZ(2) = 0.0
         XYZ(3) = 0.0
         IF( MOFACE .EQ. 4 ) THEN
            IF( LFACES(MOFACE,N) .EQ. 0 ) THEN
C              TRIANGLE P1
               NBN = 3
            ELSE
C              QUADRANGLE Q1
               NBN = 4
            ENDIF
         ELSE
            IF( LFACES(MOFACE,N) .EQ. 0 ) THEN
C              TRIANGLE P2
               NBN = 6
            ELSE
C              QUADRANGLE Q2 SANS LE BARYCENTRE
               NBN = 8
            ENDIF
         ENDIF
C
         DO 20 J=1,NBN
C           LE NUMERO DU POINT
            NS = LFACES( J, N )
            DO 10 K=1,3
               XYZ(K) = XYZ(K) + XYZPOI(K,NS)
 10         CONTINUE
 20      CONTINUE
         XYZ(1) = XYZ(1) / NBN
         XYZ(2) = XYZ(2) / NBN
         XYZ(3) = XYZ(3) / NBN
C
C        COTE Z AXONOMETRIQUE DANS LA DIRECTION DE VISEE
         CALL XYZAXO( XYZ, XYZ )
         DIFACE(N) = XYZ(3)
C
      ENDDO
C
      RETURN
      END
