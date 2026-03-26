      SUBROUTINE RAC3BGP1( NBNOMA, POL, NBDLFX, VADLFX, NODLFX, UG0,
     %                     UG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES 1 A 3 RACINES DES 2 POLYNOMES NODAUX DU 3-EME DEGRE
C -----

C ENTREES:
C --------
C NBNOMA : NOMBRE DE SOMMETS DU MAILLAGE P1
C POL    : (NBNOMA,0:3,2) LES 4 COEFFICIENTS DES 2 POLYNOMES NODAUX
C NBDLFX : NOMBRE DE  DEGRE  DE LIBERTE FIXE PAR CONDITION DE DIRICHLET
C VADLFX : VALEUR DES DEGRES DE LIBERTE FIXE PAR CONDITION DE DIRICHLET
C NODLFX : (NBNOMA,2) NODLFX(I,J)=0  SI DL LIBRE
C                                =NO DU DL FIXE DE 1 A NBDLFX
C UG0    : (NBNOMA,2) PARTIE REELLE ET IMAGINAIRE DE L'ONDE AUX NOEUDS
C          A L'ITERATION m PRECEDENTE POUR FAVORISER NEWTON

C SORTIES:
C --------
C UG     : (NBNOMA,2) PARTIE REELLE ET IMAGINAIRE DE L'ONDE AUX NOEUDS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray     Mai 2014
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  POL(NBNOMA,0:3,2), VADLFX(NBDLFX),
     %                  UG0(NBNOMA,2), UG(NBNOMA,2),
     %                  RAC2, RAC3
      INTEGER           NODLFX(NBNOMA,2)

      DO NC = 1, 2

         DO N = 1, NBNOMA

C           NO DU DL FIXE?
            NDLFX = NODLFX(N,NC)

            IF( NDLFX .GT. 0 ) THEN

C              DL FIXE => VALEUR FIXEE
               UG( N, NC ) = VADLFX( NDLFX )

            ELSE

C              DL LIBRE => UNE DES 1 A 3 RACINES DU POLYNOME AU NOEUD I
               IF(  POL(N,3,NC) .EQ. 0D0 ) THEN
                  PRINT *,'RAC3BGP1: PB COEFFICIENT A=0 NC=',NC,' N=',N
               ENDIF

               CALL RACPOL3( POL(N,3,NC), POL(N,2,NC),
     %                       POL(N,1,NC), POL(N,0,NC), UG0(N,NC),
     %                       NBRAC, UG(N,NC), RAC2, RAC3 )

C              AFFICHAGE FINAL SI 3 RACINES REELLES
               IF( NBRAC .EQ. 3 ) THEN
                  PRINT 19000, NC,N,UG0(N,NC),NBRAC,UG(N,NC),RAC2,RAC3
               ENDIF

            ENDIF

         ENDDO

      ENDDO

      RETURN
19000 FORMAT('RAC3BGP1: NC=',I1,' N=',I8,' RAC0=',G15.7,'  NBRAC=',I1,
     %       ' RAC1=',G15.7,' RAC2=',G15.7,' RAC3=',G15.7 )
      END
