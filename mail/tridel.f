       SUBROUTINE TRIDEL( NDIM, MNXYZS, MNNSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RENDRE DE TYPE DELAUNAY LA TRIANGULATION PAR DES ECHANGES DE
C -----    DIAGONALES DES 2 TRIANGLES

C ENTREES:
C --------
C NDIM   : DIMENSION (2 ou 3) DE L'ESPACE DE LA SURFACE
C MNXYZS : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF : ADRESSE DU TABLEAU NSEF      DE LA SURFACE
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Labo J-L. LIONS  UPMC  PARIS   SEPTEMBRE 2007
C2345X7..............................................................012
      PARAMETER (ITERMX=3)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      REAL          CETRIA(3)

      IERR   = 0
      RGRAND = RINFO('GRAND')
      NBEF   = MCN( MNNSEF + WBEFOB )
      ITER   = 0
      NBECHA = 0
      NOTRMIN= 1
      NOTRMAX= NBEF

C     BOUCLE SUR LES 3 ARETES DE CHAQUE TRIANGLE
C     ==========================================
 5    NBEF   = MCN( MNNSEF + WBEFOB )
      NBECH0 = 0
      MNSOEL = MNNSEF + WUSOEF
      MNXY   = MNXYZS + WYZSOM - 3
      NOTRMN = NBEF
      NOTRMX = 0
      DO 100 NOTRIA1 = NOTRMIN, NOTRMAX-1

C        POINTEUR AVANT LES NO DES SOMMETS DU TRIANGLE 1
         MNT1 = MNSOEL + 4 * NOTRIA1 - 5
         IF( MCN(MNT1+4) .GT. 0 ) GOTO 100

C        BOUCLE SUR L'ARETE NAT1 DU TRIANGLE NOTRIA1
         DO 90 NAT1=1,3

C           SOMMET 1 DE L'ARETE NAT1
            NST11 = MCN( MNT1 + NAT1 )

C           SOMMET 2 DE L'ARETE NAT1
            IF( NAT1 .EQ. 3 ) THEN
               NAT11 = 1
            ELSE
               NAT11 = NAT1 + 1
            ENDIF
            NST12 = MCN( MNT1 + NAT11 )

C           3-EME SOMMET DU TRIANGLE 1
            IF( NAT11 .EQ. 3 ) THEN
               NAT12 = 1
            ELSE
               NAT12 = NAT11 + 1
            ENDIF
            NST13 = MCN( MNT1 + NAT12 )

C           RECHERCHE DU TRIANGLE ADJACENT PAR CETTE ARETE
            DO 80 NOTRIA2 = NOTRMIN, NOTRMAX

               IF( NOTRIA2 .EQ. NOTRIA1 ) GOTO 80

C              POINTEUR AVANT LES NO DES SOMMETS DU TRIANGLE 1
               MNT2 = MNSOEL + 4 * NOTRIA2 - 5
               IF( MCN(MNT2+4) .GT. 0 ) GOTO 80

C              BOUCLE SUR LES 3 ARETES DU TRIANGLE NOTRIA2
               DO NAT2=1,3

C                 LE NO DU SOMMET 1 DE L'ARETE NAT2 DU TRIANGLE 2
                  NST21 = MCN( MNT2 + NAT2 )

C                 LE NO DU SOMMET 2 DE L'ARETE NAT2 DU TRIANGLE 2
                  IF( NAT2 .EQ. 3 ) THEN
                     NAT21 = 1
                  ELSE
                     NAT21 = NAT2 + 1
                  ENDIF
                  NST22 = MCN( MNT2 + NAT21 )

                  IF( (NST11 .EQ. NST22 .AND. NST12 .EQ. NST21) .OR.
     %                (NST11 .EQ. NST21 .AND. NST12 .EQ. NST22) ) THEN

C                    LE QUADRANGLE EST IL CONVEXE?
                     CALL QUADCXR(RMCN(MNXY+3*NST11),RMCN(MNXY+3*NST12),
     %                            RMCN(MNXY+3*NST22),RMCN(MNXY+3*NST13),
     %                            NONOUI )
                     IF( NONOUI .EQ. 0 ) GOTO 90

C                    LE SOMMET NST23 DU TRIANGLE 2 NON SUR L'ARETE COMMUNE
                     IF( NAT21 .EQ. 3 ) THEN
                        NAT22 = 1
                     ELSE
                        NAT22 = NAT21 + 1
                     ENDIF
                     NST23 = MCN( MNT2 + NAT22 )

                     IF( NDIM .EQ. 2 ) THEN

C                       SURFACE DE R2:
C                       L'ARETE COMMUNE EST ELLE DELAUNAY? <=>
C                       LE SOMMET DU TRIANGLE NOTRIA2 EST IL EN DEHORS
C                       DU CERCLE CIRCONSCRIT AU TRIANGLE NOTRIA1?
                        CALL CENCER( RMCN(MNXY+3*NST11),
     %                               RMCN(MNXY+3*NST12),
     %                               RMCN(MNXY+3*NST13), CETRIA )
                        IF( CETRIA(3) .EQ. 1E28 ) GOTO 90
                        R2 = CETRIA(3)

C                       LE SOMMET SUIVANT DU TRIANGLE 2 NON SUR
C                       L'ARETE COMMUNE EST NST23
                        X = CETRIA(1) - RMCN(MNXY+3*NST23)
                        Y = CETRIA(2) - RMCN(MNXY+3*NST23+1)
                        D = X*X + Y*Y

                        IF( D .GE. R2 ) GOTO 90

                     ELSE

C                       SURFACE DE R3:
C                       L'ARETE COMMUNE EST ELLE DELAUNAY? <=>
C                       LA SOMME MINIMALE DES RAYONS DES CERCLES
C                       CIRCONSCRITS DES 2 CONFIGURATIONS EST REALISEE
C                       LE RAYON DE LA BOULE CERCLE CIRCONSCRITE AU
C                       TRIANGLE NOT
                        CALL RACITR3D( RMCN(MNXY+3*NST11),
     %                                 RMCN(MNXY+3*NST12),
     %                                 RMCN(MNXY+3*NST13), R1, S )
                        IF( R1 .EQ. RGRAND ) GOTO 90

C                       LE RAYON DE LA BOULE CERCLE CIRCONSCRITE AU TRIANGLE NOT
                        CALL RACITR3D( RMCN(MNXY+3*NST21),
     %                                 RMCN(MNXY+3*NST22),
     %                                 RMCN(MNXY+3*NST23), R2, S )
                        IF( R2 .EQ. RGRAND ) GOTO 90

C                       LE RAYON DE LA BOULE CERCLE CIRCONSCRITE AU TRIANGLE 123
                        CALL RACITR3D( RMCN(MNXY+3*NST13),
     %                                 RMCN(MNXY+3*NST11),
     %                                 RMCN(MNXY+3*NST23), R3, S )
                        IF( R3 .EQ. RGRAND ) GOTO 90

C                       LE RAYON DE LA BOULE CERCLE CIRCONSCRITE AU TRIANGLE 134
                        CALL RACITR3D( RMCN(MNXY+3*NST13),
     %                                 RMCN(MNXY+3*NST23),
     %                                 RMCN(MNXY+3*NST12), R4, S )
                        IF( R4 .EQ. RGRAND ) GOTO 90

C                       DELAUNAY: MIN(R1+R2, R3+R4)
                        IF( R3+R4 .GE. R1+R2 ) GOTO 90

                     ENDIF

C                    LE SOMMET NST23 EST DANS LE CERCLE CIRCONSCRIT AU TRIANGLE
C                    => ECHANGE DES 2 DIAGONALES FORME PAR LES 2 TRIANGLES

ccc                     CALL ECHDIA( NDIM, NST11, NST12, NST13, NST23,
ccc     %                            NOTRIA1, NOTRIA2, MNXYZS, MNNSEF,
ccc     %                            IERR )

                     MCN( MNT1 + 1 ) = NST13
                     MCN( MNT1 + 2 ) = NST11
                     MCN( MNT1 + 3 ) = NST23
                     MCN( MNT1 + 4 ) = 0

                     MCN( MNT2 + 1 ) = NST13
                     MCN( MNT2 + 2 ) = NST23
                     MCN( MNT2 + 3 ) = NST12
                     MCN( MNT2 + 4 ) = 0

ccc                     IF( IERR .EQ. 0 ) THEN
C                       ECHANGE EFFECTIF DES DIAGONALES
                        NBECH0 = NBECH0 + 1
                        NOTRMN = MIN( NOTRMN, NOTRIA1 )
                        NOTRMN = MIN( NOTRMN, NOTRIA2 )
                        NOTRMX = MAX( NOTRMX, NOTRIA1 )
                        NOTRMX = MAX( NOTRMX, NOTRIA2 )
ccc                     ENDIF

                     GOTO 90

                  ENDIF
               ENDDO
 80         ENDDO
 90      ENDDO
 100  ENDDO

      ITER    = ITER + 1
      NBECHA  = NBECHA + NBECH0
      NOTRMIN = NOTRMN
      NOTRMAX = NOTRMX

ccc      WRITE(IMPRIM,*) 'tridel: ITERATION',ITER,
ccc     %                ' Nombre ECHANGES des DIAGONALES=',NBECH0

      IF( ITER .LT. ITERMX .AND. NBECH0 .GT. 0 ) GOTO 5

      WRITE(IMPRIM,*) 'tridel: ITERATION',ITER,
     %                ' Nombre ECHANGES des DIAGONALES=',NBECHA

      RETURN
      END
