      SUBROUTINE XYZD1TE( NBEF, NBNOEF, NUNOTE, XYZPOI, XYZ,  NTE,
     %                    DELTAe, CBTE, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     LE TETRAEDRE NTE CONTIENT IL LE POINT XYZ?
C ----

C ENTREES:	
C --------
C NBEF   : NOMBRE DE TETRAEDRES DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN TETRAEDRE
C          (SEULS LES 4 SOMMETS SONT UTILES ICI et ONT POUR No 1 2 3 4)
C NUNOTE : NUMERO DES NBNOEF NOEUDS DES NBEF TETRAEDRES
C XYZPOI : 3 COORDONNEES DES SOMMETS DES TETRAEDRES
C XYZ    : 3 COORDONNEES DU POINT DANS LE TETRAEDRE NTE OU NON

C SORTIES:
C --------
C DELTAe : 6 FOIS LE VOLUME SIGNE DU TETRAEDRE NTE
C CBTE   : 4 COORDONNEES BARYCENTRIQUES DE XYZ DANS LE TETRAEDRE NTE
C NONOUI : =1 OUI XYZ   EST     INTERNE AU TETRAEDRE NTE  CBTE CALCULE
C          =0 NON XYZ N'EST PAS INTERNE AU TETRAEDRE NTE  CBTE CALCULE
C          =-1 TETRAEDRE NTE DE VOLUME<=0  CBTE NON CALCULE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2020
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      INTEGER           NBEF,   NBNOEF, NUNOTE(NBEF,NBNOEF), NTE, NONOUI
      REAL              XYZPOI(3,*)
      DOUBLE PRECISION  XYZ(3), CBTE(4)
      DOUBLE PRECISION  Xe(3,4), DELTAe, X1, Y1, Z1, DETM33
      INTEGER           I, K, NS

C     RECUPERATION DES COORDONNEES Xe(3,4) DES 4 SOMMETS DU TETRAEDRE NTE
      DO I=1,4
         NS = NUNOTE( NTE, I )
         DO K=1,3
            Xe( K, I ) = DBLE( XYZPOI( K, NS ) )
         ENDDO
      ENDDO

C     CALCUL DU JACOBIEN DE Fe POLYNOMIALE DE DEGRE 1 du TETRAEDRE NTE
C     ou 6 FOIS LE VOLUME du TETRAEDRE NTE
      X1 = Xe(1,1)
      Y1 = Xe(2,1)
      Z1 = Xe(3,1)
      DELTAe = DETM33( Xe(1,2)-X1, Xe(1,3)-X1, Xe(1,4)-X1,
     %                 Xe(2,2)-Y1, Xe(2,3)-Y1, Xe(2,4)-Y1,
     %                 Xe(3,2)-Z1, Xe(3,3)-Z1, Xe(3,4)-Z1 )

      IF( DELTAe .LE. 0D0 ) THEN
C        RENCONTRE DU TETRAEDRE NTE DE VOLUME<0
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *,'xyzd1te: ATTENTION TETRAEDRE',NTE,
     %              ' de VOLUME*6=',DELTAe
         ELSE
            PRINT *,'xyzd1te: ATTENTION TETRAHEDRON',NTE,
     %              ' of VOLUME*6=',DELTAe,' <=0'
         ENDIF
         NONOUI = -1
         GOTO 9999
      ENDIF

C     CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT XYZ DANS NTE
C     Lambda 1
      CBTE(1) = DETM33( Xe(1,2)-XYZ(1), Xe(1,3)-XYZ(1), Xe(1,4)-XYZ(1),
     %                  Xe(2,2)-XYZ(2), Xe(2,3)-XYZ(2), Xe(2,4)-XYZ(2),
     %                  Xe(3,2)-XYZ(3), Xe(3,3)-XYZ(3), Xe(3,4)-XYZ(3))
     %        / DELTAe
C     Lambda 2
      CBTE(2) = DETM33( Xe(1,3)-XYZ(1), Xe(1,1)-XYZ(1), Xe(1,4)-XYZ(1),
     %                  Xe(2,3)-XYZ(2), Xe(2,1)-XYZ(2), Xe(2,4)-XYZ(2),
     %                  Xe(3,3)-XYZ(3), Xe(3,1)-XYZ(3), Xe(3,4)-XYZ(3))
     %        / DELTAe
C     Lambda 3
      CBTE(3) = DETM33( Xe(1,4)-XYZ(1), Xe(1,1)-XYZ(1), Xe(1,2)-XYZ(1),
     %                  Xe(2,4)-XYZ(2), Xe(2,1)-XYZ(2), Xe(2,2)-XYZ(2),
     %                  Xe(3,4)-XYZ(3), Xe(3,1)-XYZ(3), Xe(3,2)-XYZ(3))
     %        / DELTAe
C     Lambda 4
      CBTE(4) = DETM33( Xe(1,1)-XYZ(1), Xe(1,3)-XYZ(1), Xe(1,2)-XYZ(1),
     %                  Xe(2,1)-XYZ(2), Xe(2,3)-XYZ(2), Xe(2,2)-XYZ(2),
     %                  Xe(3,1)-XYZ(3), Xe(3,3)-XYZ(3), Xe(3,2)-XYZ(3))
     %        / DELTAe

C     XYZ INTERNE au TETRAEDRE NTE?
      IF( CBTE(1) .GE. -0.01D0 .AND. CBTE(1) .LE. 1.01D0 .AND.
     %    CBTE(2) .GE. -0.01D0 .AND. CBTE(2) .LE. 1.01D0 .AND.
     %    CBTE(3) .GE. -0.01D0 .AND. CBTE(3) .LE. 1.01D0 .AND.
     %    CBTE(4) .GE. -0.01D0 .AND. CBTE(4) .LE. 1.01D0 ) THEN

C        OUI: XYZ est INTERNE au TETRAEDRE NTE
C        ----------------------------------------
         NONOUI = 1

      ELSE

C        NON: XYZ est EXTERNE au TETRAEDRE NTE
C        ----------------------------------------
         NONOUI = 0

      ENDIF

 9999 RETURN
      END
