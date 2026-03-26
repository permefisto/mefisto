      SUBROUTINE BASTCOPL( COSMAXPL, NBSOMM,  XYZSOM,
     %                     MXTRIA,   M1TRIA6, NBTRIA, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    BARYCENTRAGE DES SOMMETS D'UNE TRIANGULATION
C -----    SI TOUS SES TRIANGLES SONT COPLANAIRES a
C          UN ANGLE de COPLANEARTE PRES
C          REDUCTION DU NOMBRE DE TRIANGLES COPLANAIRES D'UN SOMMET
C          PAR RETRIANGULATION DE SES TRIANGLES et SA SUPPRESSION

C ENTREES:
C --------
C COSMAXPL:COSINUS DE L'ANGLE DIEDRE ENTRE LES 2 PLANS DES TRIANGLES
C          AU DESSOUS DUQUEL LES 2 TRIANGLES SONT CONSIDERES NON-COPLANAIRES
C NBSOMM : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : XYZ DES NBSOMM SOMMETS DE LA TRIANGULATION
C MXTRIA : NOMBRE MAXIMUM DE TRIANGLES DU TABLEAU NOTRIA
C NBTRIA : NOMBRE DE TRIANGLES
C NOTRIA : NUMERO DE 1 A NBSOMM DES SOMMETS DES TRIANGLES DU MAILLAGE
C          1 2 3 LE NO XYZSOM DES 3 SOMMETS DU TRIANGLE
C          4 5 6 LE NO DU TRIANGLE ADJACENT PAR L'ARETE 1 2 3

C MODIFIE:
C --------
C XYZSOM: XYZ DES NBSOMM SOMMETS + BARYCENTRES DE LA TRIANGULATION
C                                - SOMMETS SUPPRIMES PAR RETRIANGULATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Veulettes sur mer               Fevrier 2020
C23456...............................................................012
      PARAMETER       ( MXTRST=256 )
C     MXTRST=MAXIMUM DE TRIANGLES D'UN MEME SOMMET
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      INTEGER  NOTRIA(M1TRIA6,MXTRIA), NOTRST(MXTRST), NOARSI(3,MXTRST)
      REAL     XYZSOM(3,NBSOMM), XYZBAR(3)

      TRACTE0 = TRACTE
C     NOMBRE DE SUPPRESSION DE SOMMETS PAR RETRIANGULATION DE SES TRIANGLES
      NBSUST  = 0

C     1-ERE ITERATION EN NUMEROS DE SOMMET DECROISSANTS
      NSMIN = NBSOMM
      NSMAX = 1
      LEPAS =-1

      DO 200 ITER = 1, 4

      NBBAR  = 0
      NB2DIA = 0

      DO 100 NoS = NSMIN, NSMAX, LEPAS

C        NST EST UN SOMMET INTERNE DEPLACE AU BARYCENTRE DE SES VOISINS
         NST = NoS

C        INITIALISATION DES XYZ DU BARYCENTRE DES TRIANGLES DE SOMMET NST
 5       XYZBAR( 1 ) = 0.
         XYZBAR( 2 ) = 0.
         XYZBAR( 3 ) = 0.
         NBSTBA = 0
         NBTRST = 0

C        RECHERCHE DES TRIANGLES DE SOMMET NST
         DO 50 NTR = 1, NBTRIA

            DO I = 1, 3

C              NUMERO DU SOMMET I DU TRIANGLE NTR
               NSTI = NOTRIA( I, NTR )

               IF( NSTI .LE. 0 ) GOTO 50

               IF( NSTI .EQ. NST ) THEN

C                 NTR A POUR SOMMET NST
C                 LES 2 TRIANGLES ADJACENTS DE SOMMET NST
C                 SONT ILS COPLANAIRES?
                  II = 0
                  IF( I .EQ. 3 ) THEN
                     I1 = 1
                  ELSE
                     I1 = I+1
                  ENDIF

C                 NO NOTRIA DU TRIAGLE ADJACENT A L'ARETE I DE NTR
                  NTRADJ = NOTRIA( 3+I, NTR )

 10               IF( NTRADJ .LE. 0 ) GOTO 30
C                 LE TRIANGLE NTR JOUXTE LA FRONTIERE

C                 RECHERCHE DU NO NS4 DE SOMMET DE NTRADJ NON SOMMET DE NTR
                  DO K=1,3

                     NS4 = NOTRIA( K, NTRADJ )
                     IF( NS4 .NE. NOTRIA( 1, NTR ) .AND.
     %                   NS4 .NE. NOTRIA( 2, NTR ) .AND.
     %                   NS4 .NE. NOTRIA( 3, NTR ) ) THEN

C                       LE COSINUS DE L'ANGLE DES NORMALES AUX
C                       TRIANGLES NTR NTRADJ
                        CALL COS2TR( XYZSOM( 1, NOTRIA( 1, NTR ) ),
     %                               XYZSOM( 1, NOTRIA( 2, NTR ) ),
     %                               XYZSOM( 1, NOTRIA( 3, NTR ) ),
     %                               XYZSOM( 1, NOTRIA( 2, NTR ) ),
     %                               XYZSOM( 1, NOTRIA( 1, NTR ) ),
     %                               XYZSOM( 1, NS4 ),
     %                               COS2PL, IERR1, IERR2 )

                        IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) GOTO 20
C                       SI UN TRIANGLE EST SANS NORMALE
C                       ALORS LA COPLANEARITE EST FAVORISEE

                        IF( COS2PL .LT. COSMAXPL ) GOTO 100
C                       SI LES 2 TRIANGLES NE SONT PAS COPLANAIRES
C                       ALORS ABANDON DU BARYCENTRAGE DU SOMMET NST

C                       AJOUT DES SOMMETS DE NTR NON NS AU BARYCENTRE
 20                     DO L = 1, 3

C                          NUMERO DU SOMMET L DU TRIANGLE NTR
                           NSTL = NOTRIA( L, NTR )
                           IF( NSTL .GT. 0 .AND. NSTL .NE. NST ) THEN
                              DO M=1,3
                                 XYZBAR(M) = XYZBAR(M) + XYZSOM(M,NSTL)
                              ENDDO
                              NBSTBA = NBSTBA + 1
                           ENDIF

                        ENDDO

                     ENDIF

                  ENDDO

 30               IF( II .EQ. 0 ) THEN
C                    II L'ARETE PRECEDENT L'ARETE I DE NTR
C                    LES 2 ARETES ONT POUR SOMMET NST
                     IF( I .EQ. 1 ) THEN
                        II = 3
                     ELSE
                        II = I-1
                     ENDIF
C                    NO NOTRIA DU TRIAGLE ADJACENT A L'ARETE II DE NTR
                     NTRADJ = NOTRIA( 3+II, NTR )
                     GOTO 10
                  ENDIF

C                 NTR UN TRIANGLE PARTICIPANT AU BARYCENTRAGE DE NS
                  IF( NBTRST .LT. MXTRST ) THEN
                     NBTRST = NBTRST + 1
                     NOTRST( NBTRST ) = NTR
                  ENDIF

                  GOTO 50

               ENDIF
            ENDDO

C           PASSAGE AU TRIANGLE DE SOMMET NSt SUIVANT
 50      ENDDO

         IF( NBSTBA .GT. 0 ) THEN

C           NST EST PLACE AU BARYCENTRE DE SES NBSTBA VOISINS COPLANAIRES
C           ET DE LUI-MEME
            DO K=1,3
               XYZSOM(K,NST) = ( XYZSOM(K,NST) + XYZBAR(K)/NBSTBA ) / 2
            ENDDO
            NBBAR = NBBAR + 1

         ENDIF

cccC        TRACE DE NBTRST TRIANGLES DE NOTRIA, LEURS VOISINS et
cccC        LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
ccc         CALL TRTRIAN( 'bastcopl', XYZSOM, 6, MXTRST, NBTRST, NOTRST,
ccc     %                  NOTRIA )

         IF( NBTRST .GE. 8 ) THEN

C           ESSAI DE RETRIANGULATION DES NBTRST TRIANGLES COPLANAIRES
C           DE SOMMET NST EN JOIGNANT LES SOMMETS DES ARETES SIMPLES
C           DES NBTRCF TRIANGLES ET EN SUPPRIMANT LE SOMMET NST
            CALL SUTR1ST( COSMAXPL, NBSOMM, XYZSOM,
     %                    MXTRIA,   NBTRIA, NOTRIA,
     %                    NST,      MXTRST, NBTRST, NOTRST,
     %                    MXTRST,   NBARSI, NOARSI,
     %                    NBTNEW,   NB2DI,  NS1,    IERR )
            NB2DIA = NB2DIA + NB2DI
            IF( NBTNEW .GT. 0 .AND. NS1 .GT. 0 ) THEN
               NBSUST = NBSUST + 1
ccc            IF( NSt.EQ.NoS .AND. NS1.GT.0 .AND. NS1.NE.NoS ) THEN
               IF( NS1.NE.NSt .AND. NS1.GT.0 ) THEN
C                 BARYCENTRAGE DU SOMMET NS1 et RETRIANGULATION EVENTUELLE
                  NSt = NS1
                  GOTO 5
               ENDIF
            ENDIF

         ENDIF

         tracte = tracte0

 100  ENDDO

      PRINT*,'bastcopl: Iteration',ITER,': BARYCENTRAGE de',NBBAR,
     %       ' SOMMETS de TRIANGLES COPLANAIRES.',NBSUST,
     %       ' SOMMETS SUPPRIMES et',NB2DIA,
     %       ' ECHANGES de DIAGONALES 2Tr<->2Tr'

      IF( NBBAR .EQ. 0 .AND. NBSUST .EQ. 0 ) GOTO 9999

C     SECONDE ITERATION EN NUMEROS DE SOMMET CROISSANT ...
      I     = NSMIN
      NSMIN = NSMAX
      NSMAX = I
      LEPAS = -LEPAS

 200  ENDDO


 9999 TRACTE = TRACTE0

      RETURN
      END
