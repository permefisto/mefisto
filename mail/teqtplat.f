      SUBROUTINE TEQTPLAT( PTXYZD, NOSOTE,  NTYPQT, NS1, NS2, NS3, NS4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETERMINER SI LE TETRAEDRE NOSOTE (SUPPOSE QUASI-PLAT) EST
C -----    DE TYPE TRIANGLE OU QUADRANGLE A 2 DIAGONALES NS1-NS2 NS3-NS4

C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NOSOTE : NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE A ESTIMER

C SORTIES :
C ---------
C NTYPQT  : =1 TETRAEDRE QUASI-PLAT DE TYPE QUADRANGLE AVEC 2 DIAGONALES
C              DE SOMMETS NS1-NS2 et NS3-NS4 (de 1 a 4) DANS NOSOTE
C           =2 TETRAEDRE QUASI-PLAT DE TYPE TRIANGLE AVEC UN SOMMET NS1
C              DE FACE OPPOSEE NS1+1 S'Y PROJETANT EN SON INTERIEUR
C              (NS2=NS3=NS4=0)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Janvier 2017
C2345X7..............................................................012
      DOUBLE PRECISION  PTXYZD(4,*), SURTRD, PTPROJ(3),
     %                  S123, SP23, SP31, SP12,
     %                  MILIEU(3,6), DISMIL(3), DISMIN
      INTEGER           NOSOTE(4)

C     NO DES MILIEU D'ARETE A MILIEU D'ARETE OPPOSEE
      INTEGER           NOMIMITE(2,3)
      DATA              NOMIMITE / 1,6, 2,4, 3,5 /

      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /

      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

C     LE TETRAEDRE EST UN QUADRANGLE QUASI-PLAT SI LA PROJECTION
C     DE CHAQUE SOMMET SUR LE TRIANGLE OPPOSE LUI EST EXTERIEURE
      NS4 = 4
      DO 10 NF = 1, 4

C        PTPROJ EST LE POINT DE PROJECTION DU SOMMET NF-1 DE NOSOTE SUR LA FACE NF
         NS1 = NOSOFATE( 1, NF )
         NS2 = NOSOFATE( 2, NF )
         NS3 = NOSOFATE( 3, NF )

C        VALEUR ABSOLUE DE LA SURFACE DE LA FACE NF DU TETRAEDRE NOSOTE
         S123 = SURTRD( PTXYZD( 1, NOSOTE( NS1 ) ),
     %                  PTXYZD( 1, NOSOTE( NS2 ) ),
     %                  PTXYZD( 1, NOSOTE( NS3 ) ) )
         IF( S123 .LE. 0D0 ) THEN
            PRINT*,'teqtplat: FACE',NF,'de SURFACE NULLE du TETRAEDRE',
     %              NOSOTE
            GOTO 10
         ENDIF

C        RECHERCHE DU POINT PROJECTION DU SOMMET NS4 SUR LE TRIANGLE NS1 NS2 NS3
         CALL PRPTPLD( PTXYZD( 1, NOSOTE( NS4 ) ),
     %                 PTXYZD( 1, NOSOTE( NS1 ) ),
     %                 PTXYZD( 1, NOSOTE( NS2 ) ),
     %                 PTXYZD( 1, NOSOTE( NS3 ) ),
     %                 PTPROJ, IERR )

         IF( IERR .NE. 0 ) THEN
            GOTO 9000
         ENDIF

C        VALEUR ABSOLUE DES COORDONNEES BARYCENTRIQUES DE PRPROJ
C        DANS LE TRIANGLE FACE NF
         SP23  = SURTRD( PTPROJ,
     %                   PTXYZD( 1, NOSOTE( NS2 ) ),
     %                   PTXYZD( 1, NOSOTE( NS3 ) ) )
         IF( SP23 .LE. 0D0 ) THEN
            PRINT*,'teqtplat: la PROJECTION du SOMMET',NS4,' EST SUR L''
     %ARETE',NS2,NS3,' de',NOSOTE
            GOTO 9000
         ENDIF

         SP31  = SURTRD( PTPROJ,
     %                   PTXYZD( 1, NOSOTE( NS3 ) ),
     %                   PTXYZD( 1, NOSOTE( NS1 ) ) )
         IF( SP31 .LE. 0D0 ) THEN
            PRINT*,'teqtplat: la PROJECTION du SOMMET',NS4,' EST SUR L''
     %ARETE',NS3,NS1,' de',NOSOTE
            GOTO 9000
         ENDIF

         SP12  = SURTRD( PTPROJ,
     %                   PTXYZD( 1, NOSOTE( NS1 ) ),
     %                   PTXYZD( 1, NOSOTE( NS2 ) ) )
         IF( SP12 .LE. 0D0 ) THEN
            PRINT*,'teqtplat: la PROJECTION du SOMMET',NS4,' EST SUR L''
     %ARETE',NS1,NS2,' de',NOSOTE
            GOTO 9000
         ENDIF

ccc         IF( SP23+SP31+SP12 .LE. S123*1.001D0 ) THEN
         IF( SP23+SP31+SP12 .LE. S123*1.01D0 ) THEN

C           PTPROJ EST INTERIEUR A LA FACE NF DU TETRAEDRE
C           LE TETRAEDRE NOSOTE EST DE TYPE TRIANGLE SURMONTE D'UN SOMMET NS4
            GOTO 9000

         ENDIF

         NS4 = NF

 10   ENDDO

C     TETRAEDRE NOSOTE DE TYPE QUADRANGLE
C     -----------------------------------
C     TOUS LES SOMMETS SE PROJETTENT A L'EXTERIEUR DE LEUR FACE OPPOSEE
C     IL RESTE A DEFINIR SES 2 DIAGONALES DE SOMMETS NS1-NS2 et NS3-NS4
      DO K = 1, 6
         N1 = NOSOARTE( 1, K )
         N2 = NOSOARTE( 2, K )
         DO L=1,3
            MILIEU(L,K)=(PTXYZD(L,NOSOTE(N1))+ PTXYZD(L,NOSOTE(N2)) ) /2
         ENDDO
      ENDDO

      DO K = 1 , 3
         N1 = NOMIMITE( 1, K )
         N2 = NOMIMITE( 2, K )
C        CARRE DE LA DISTANCE MILIEU ARETE N1 - MILIEU ARETE N2
         DISMIL(K) = ( MILIEU(1,N1) - MILIEU(1,N2) ) ** 2
     %             + ( MILIEU(2,N1) - MILIEU(2,N2) ) ** 2
     %             + ( MILIEU(3,N1) - MILIEU(3,N2) ) ** 2
      ENDDO

C     RECHERCHE DE LA PLUS PETITE DISTANCE ENTRE MILIEUX D'ARETES
      DISMIN = DISMIL(1)
      N1 = 1
      DO K = 2 , 3
         IF( DISMIL(K) .LT. DISMIN ) THEN
            N1 = K
            DISMIN = DISMIL(K)
         ENDIF
      ENDDO

C     DIAGONALES NS1-NS2 et NS3-NS4 FINALES
      GOTO ( 101, 102, 103 ) , N1

C     MILIEUX ARETES 1-6 => DIAGONALES NS1-NS2 et NS3-NS4
 101  NS1 = 1
      NS2 = 2
      NS3 = 3
      NS4 = 4
      GOTO 110

C     MILIEUX ARETES 2-4 => DIAGONALES NS1-NS2 et NS3-NS4
 102  NS1 = 2
      NS2 = 3
      NS3 = 1
      NS4 = 4
      GOTO 110

C     MILIEUX ARETES 3-5 => DIAGONALES NS1-NS2 et NS3-NS4
 103  NS1 = 3
      NS2 = 1
      NS3 = 2
      NS4 = 4

 110  NTYPQT = 1
      GOTO 9999

C     TETRAEDRE QUASI-PLAT DE TYPE TRIANGLE
C     -------------------------------------
C     AVEC UN SOMMET CENTRAL NS1 DE NOSOTE INTERIEUR A SA FACE OPPOSEE
C     OU SUR L'UNE DE SES 3 ARETES
 9000 NS1 = NS4
      NS2 = 0
      NS3 = 0
      NS4 = 0

      NTYPQT = 2

 9999 RETURN
      END
