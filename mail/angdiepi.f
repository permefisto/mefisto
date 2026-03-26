      SUBROUTINE ANGDIEPI( COSMAXPL, NBSOMM, XYZSOM, NBTRIA, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRESSION DES ANGLES DIEDRES PROCHE DE Pi DES COUPLES DE
C -----    TRIANGLES de la TRIANGULATION PAR ECHANGE FORCE DE DIAGONALES
C          et SI LES SURFACES DES 2 TRIANGLES SONT DISPROPORTIONNEES

C ENTREES:
C --------
C COSMAXPL:COSINUS DE L'ANGLE DIEDRE ENTRE LES 2 PLANS DES TRIANGLES
C          AU DESSOUS DUQUEL LES 2 TRIANGLES SONT CONSIDERES NON-COPLANAIRES
C NBSOMM : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : XYZ DES NBSOMM SOMMETS DE LA TRIANGULATION
C NBTRIA : NOMBRE DE TRIANGLES DU TABLEAU NOTRIA

C MODIFIES:
C ---------
C NOTRIA : NUMERO DE 1 A NBSOMM DES SOMMETS DES TRIANGLES DE LA TRIANGULATION
C          1 2 3 LE NO XYZSOM DES 3 SOMMETS DU TRIANGLE
C          4 5 6 LE NO DU TRIANGLE ADJACENT PAR L'ARETE 1 2 3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint Pierre du Perray             Mars 2020
C23456...............................................................012
      PARAMETER       ( MXTRAR=256 )
C     MXTRAR=MAXIMUM DE TRIANGLES D'UN MEME SOMMET
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      INTEGER  NOTRIA(6,NBTRIA), NOTRAR(MXTRAR)
      REAL     XYZSOM(3,NBSOMM)

      TRACTE0 = TRACTE

C     NOMBRE D'ARETES MODIFIEES
      NBARMO   = 0
      RSFMIN   = 1E27
      COS2TMIN = 2.
      DO 100 NTR = 1, NBTRIA

C        NOMBRE DE TRIANGLES DE SOMMET NST
         IF( NOTRIA( 1, NTR ) .LE. 0 ) GOTO 100

C        SURFACE DU TRIANGLE NTR
         STR = SURTRR( XYZSOM( 1, NOTRIA(1,NTR) ),
     %                 XYZSOM( 1, NOTRIA(2,NTR) ),
     %                 XYZSOM( 1, NOTRIA(3,NTR) ) )

         DO 50 NA = 1, 3

C           NUMERO DU TRIANGLE ADJACENT A NTR PAR L'ARETE NA
            NTRADJ = NOTRIA( 3+NA, NTR )
            IF( NTRADJ .LE. 0 ) GOTO 50
            IF( NOTRIA(1,NTRADJ) .LE. 0 ) GOTO 50

C           SURFACE DU TRIANGLE NTRADJ
            STA = SURTRR( XYZSOM( 1, NOTRIA(1,NTRADJ) ),
     %                    XYZSOM( 1, NOTRIA(2,NTRADJ) ),
     %                    XYZSOM( 1, NOTRIA(3,NTRADJ) ) )

C           RAPPORT DES SURFACES DES TRIANGLES NTR et NTRADJ
            IF( STA .GE. STR ) THEN
               RSF = STR / STA
            ELSE
               RSF = STA / STR
            ENDIF
            RSFMIN = MIN( RSFMIN, RSF )

            IF( RSF .GE. 0.0667 ) THEN
C              LES 2 SURFACES NE SONT PAS DISPROPORTIONNEES
               GOTO 50
            ENDIF

C           TRACE DES TRIANGLES NOTRAR(1:NBTRAR) DE NOTRIA, LEURS VOISINS et
C           LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
            NBTRAR = 2
            NOTRAR( 1 ) = NTR
            NOTRAR( 2 ) = NTRADJ
            CALL TRTRIAN( 'angdiepi 0',XYZSOM,6,MXTRAR,NBTRAR,NOTRAR,
     %                     NOTRIA )

C           LES 2 SOMMETS DE L'ARETE NA DE NTR
            IF( NA .EQ. 3 ) THEN
               N1 = 1
            ELSE
               N1 = NA+1
            ENDIF

            IF( N1 .EQ. 3 ) THEN
               N2 = 1
            ELSE
               N2 = N1+1
            ENDIF

            NS1 = NOTRIA( NA, NTR )
            NS2 = NOTRIA( N1, NTR )
            NS3 = NOTRIA( N2, NTR )

C           LE 3-EME SOMMET DU TRIANGLE ADJACENT PAR L'ARETE NA
            DO N3=1,3
               NS3TRA = NOTRIA( N3, NTRADJ )
               IF( NS3TRA .NE. NS1 .AND. NS3TRA .NE. NS2 ) GOTO 10
            ENDDO

C           COSINUS DE L'ANGLE DIEDRE DES 2 TRIANGLES AUTOUR DE L'ARETE NA
C           COSINUS COMPRIS ENTRE -1. et 1. de l'ANGLE des NORMALES
C           AUX TRIANGLES 1  et 2
 10         CALL COS2TR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %                   XYZSOM(1,NS2), XYZSOM(1,NS1), XYZSOM(1,NS3TRA),
     %                   COS2T, IERR1, IERR2 )
            COS2TMIN = MIN( COS2TMIN, COS2T )

            IF( IERR1  .NE. 0 ) GOTO 20
            IF( IERR2  .NE. 0 ) GOTO 20

ccc            IF( COS2T .GE. 1.-COSMAXPL ) GOTO 50
ccc            IF( COS2T .GE. -0.5 ) GOTO 50

C           ECHANGE FORCE DE L'ARETE NS1-NS2 PAR L'ARETE NS3-NS3TRA
C           POUR LE COUPLE DE TRIANGLES NTR-NTRADJ D'ARETE COMMUNE NS1-NS2
C           --------------------------------------------------------------
C           ANGLE DES 2 NORMALES des TRIANGLES NTR et NTRADJ PROCHE DE Pi
 20         PRINT*,'angdiepi: ANGLE des 2 NORMALES',
     %              ACOS(COS2T)*45./ATAN(1.),
     %    ' PROCHE DE 180 degres -> ECHANGE DIAGONALES des 2 TRIANGLES',
     %      NTR,NTRADJ
            PRINT*,'angdiepi: TRIANGLE',NTR,   ': St:',NS1,NS2,NS3
            PRINT*,'angdiepi: TRIANGLE',NTRADJ,': St:',NS2,NS1,NS3TRA

            NBARMO = NBARMO + 1

C           TRACE DES TRIANGLES NOTRAR(1:NBTRAR) DE NOTRIA, LEURS VOISINS et
C           LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
            NBTRAR = 2
            NOTRAR( 1 ) = NTR
            NOTRAR( 2 ) = NTRADJ
            CALL TRTRIAN( 'angdiepi 1',XYZSOM,6,MXTRAR,NBTRAR,NOTRAR,
     %                     NOTRIA )

C           ECHANGE FORCE DE L'ARETE NS1-NS2 PAR L'ARETE NS3-NS3TRA
            CALL EC2DFORC( NTR, NA, NBTRIA, NOTRIA, NTRADJ )

C           TRACE DES TRIANGLES NOTRAR(1:NBTRAR) DE NOTRIA, LEURS VOISINS et
C           LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
            tracte = .true.
            CALL TRTRIAN( 'angdiepi 2',XYZSOM,6,MXTRAR,NBTRAR,NOTRAR,
     %                     NOTRIA )
            tracte = tracte0

            GOTO 100

 50      ENDDO

 100  ENDDO

      PRINT*,'angdiepi:',NBARMO,' ARETES MODIFIEES car ANGLE DIEDRE TROP
     % GRAND. MIN Rapport Surfaces=',RSFMIN,' MIN COS2T=',COS2TMIN,
     % ' COSMAXPL=',COSMAXPL

      TRACTE = TRACTE0
      RETURN
      END
