      SUBROUTINE MATAXO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CONSTRUCTION DE LA MATRICE DE LA PROJECTION AXONOMETRIQUE
C -----  A PARTIR DES 3 COORDONNEES DU POINT VISE PTV
C                 DES 3 COORDONNEES POSITION DE L'OEIL
C
C        VERSION xvue
C
C SORTIE :
C --------
C AXOMAT(3,3) DANS LE COMMON TRVAR4 DE ./incl/trvari.inc"
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C2345X789............................................................012
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
      IF( NOTYVI .EQ. 0 ) THEN
C
C        PAS D'AXONOMETRIE
C        =================
C        LA MATRICE IDENTITE
         CALL MATAXI
C
      ELSE IF( NOTYVI .EQ. 11 ) THEN
C
C        AXONOMETRIE PAR PTV, OEIL
C        =========================
C        LA MATRICE DE PASSAGE XYZ => AXOXYZ
         CALL XYZIDE( AXOPTV, AXOEIL, I )
         IF( I .NE. 0 ) THEN
C
C           PTV = OEIL
C           ----------
C           LONGITUDE  30 DEGRES  LATITUDE 20 DEGRES  IMPOSEES
            CALL LONLAT( 30.0, 20.0 )
C
         ENDIF
C
C        PTV DIFFERENT DE OEIL
C        ---------------------
         AXOMAT(1,3) = AXOEIL(1) - AXOPTV(1)
         AXOMAT(2,3) = AXOEIL(2) - AXOPTV(2)
         AXOMAT(3,3) = AXOEIL(3) - AXOPTV(3)
         AXODIS = SQRT(AXOMAT(1,3)**2+AXOMAT(2,3)**2+AXOMAT(3,3)**2)
         AXOMAT(1,3) = AXOMAT(1,3) / AXODIS
         AXOMAT(2,3) = AXOMAT(2,3) / AXODIS
         AXOMAT(3,3) = AXOMAT(3,3) / AXODIS
C
         AXOMAT(1,1) =-AXOMAT(2,3)
         AXOMAT(2,1) = AXOMAT(1,3)
         AXOMAT(3,1) = 0
         R = SQRT(AXOMAT(1,1)**2 + AXOMAT(2,1)**2 )
         IF( R .GE. 1E-3 ) THEN
C
C           z N'EST PAS COLINEAIRE A PTV - OEIL
C           -----------------------------------
            AXOMAT(1,1) = AXOMAT(1,1) / R
            AXOMAT(2,1) = AXOMAT(2,1) / R
C
            AXOMAT(1,2) = -AXOMAT(2,1) * AXOMAT(3,3)
            AXOMAT(2,2) =  AXOMAT(1,1) * AXOMAT(3,3)
            AXOMAT(3,2) =  AXOMAT(2,1) * AXOMAT(1,3)
     %                   - AXOMAT(1,1) * AXOMAT(2,3)
C
         ELSE IF( AXOEIL(3) .GE. AXOPTV(3) ) THEN
C
C           z EST COLINEAIRE A PTV - OEIL   OEIL AU DESSUS DE PTV
C           -----------------------------------------------------
C           e1=x
C           e2=y
C           e3=z
            CALL MATAXI
         ELSE
C
C           z EST COLINEAIRE A PTV - OEIL   OEIL AU DESSOUS DE PTV
C           ------------------------------------------------------
C           e1=-x
C           e2= y
C           e3=-z
            CALL MATAXI
            AXOMAT(1,1) = -1.0
            AXOMAT(3,3) = -1.0
         ENDIF
      ENDIF
C
CCC      WRITE(IMPRIM,*) 'MATRICE AXOMAT EN SORTIE DU SP MATAXO'
CCC      DO 10 I=1,3
CCC         WRITE(IMPRIM,10010) (AXOMAT(I,J),J=1,3)
CCC 10   CONTINUE
CCC10010 FORMAT(4E20.6)
C
      RETURN
      END
