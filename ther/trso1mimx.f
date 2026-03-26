      SUBROUTINE TRSO1MIMX( NBFACE, NOFACE,
     %                      MOFACE, MXFACE, LFACES,
     %                      NBCOOR, XYZPOI,
     %                      NBCOMP, NCAS0, NCAS1,
     %                      NTYP,   SOLUTION, dptemp,
     %                      SOLMIN, SOLMAX  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE MIN ET MAX DE LA SOLUTION SUR LES FACES DES
C -----    SURFACES NOMMEES DE L'OBJET
C          DEFINIR LA VISEE POUR ENCADRER LE MIN MAX DES COORDONNEES
C          DES SURFACES NOMMEES DE NUMEROS DES FACES DANS LE TABLEAU NOFACE

C ENTREES:
C --------
C NBFACE : NOMBRE DE FACES DES SURFACES NOMMEES
C NOFACE : NUMERO DES FACES SELON LEUR DISTANCE DECROISSANTE A L'OEIL
C MOFACE : NOMBRE D'ENTIERS PAR FACE DANS LE TABLEAU LFACES
C MXFACE : MAJORATION DU NOMBRE DE FACES DU TABLEAU LFACES
C LFACES : TABLEAU NUMERO DES FACES DES SURFACES NOMMEES
C          LFACES(1...MOFACE,I)= NO DU POINT DE LA FACE I
C          LFACES(MOFACE,I)= 0 SI QUADRANGLE

C NBCOOR : NOMBRE DE COORDONNEES DES POINTS (3 ou 6)
C XYZPOI : COORDONNEES DES NOEUDS DU MAILLAGE DE L'OBJET

C NBCOMP : NOMBRE DE COMPOSANTES D'UN VECTEUR SOLUTION
C NCAS0  : NUMERO MINIMUM DES VECTEURS SOLUTION A TRACER
C NCAS1  : NUMERO MAXIMUM DES VECTEURS SOLUTION A TRACER
C NTYP    : =0 EMPLOI de SOLUTION
C           =1 EMPLOI de dptemp
C SOLUTION: VECTEUR"SOLUTION(NBCOMP,NCAS0:NCAS1)
C dptemp  : TABLEAU DES NCAS0:NCAS1 TABLEAU(NBCOMP) SOLUTIONS
C SOLUTION:VECTEURS SOLUTION NCAS0 a NCAS1 de NBCOMP COMPOSANTES
C SOLMIN : SOLUTION MINIMALE POUR LE CAS A TRACER
C SOLMAX : SOLUTION MAXIMALE POUR LE CAS A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A&M University at QATAR        Mars 2012
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           NOFACE(1:NBFACE)
      INTEGER           LFACES(1:MOFACE,1:MXFACE)
      DOUBLE PRECISION  SOLUTION(NBCOMP,NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp
      REAL              XYZPOI(1:NBCOOR,1:*)
      REAL              SOLMIN, SOLMAX
C
C     LE PARCOURS DES FACES
      SOLMIN = 1E28
      SOLMAX =-1E28
C
C     INITIALISATION DU MIN ET MAX DES COORDONNEES
      DO K=1,3
         COOEXT(K,1)= 1E28
         COOEXT(K,2)=-1E28
      ENDDO
C
      DO NF = 1, NBFACE
C
C        LE NUMERO DE LA FACE
         NUFA = NOFACE( NF )
C
C        LE NOMBRE DE POINTS DE LA FACE
         IF( MOFACE .EQ. 4 ) THEN
            IF( LFACES(MOFACE,NUFA) .EQ. 0 ) THEN
C              TRIANGLE P1
               NBN = 3
            ELSE
C              QUADRANGLE Q1
               NBN = 4
            ENDIF
         ELSE
            IF( LFACES(MOFACE,NUFA) .EQ. 0 ) THEN
C              TRIANGLE P2
               NBN = 6
            ELSE
C              QUADRANGLE Q2 SANS LE BARYCENTRE
               NBN = 8
            ENDIF
         ENDIF
C
C        LA SOLUTION AUX NBN POINTS DE LA FACE
         DO I=1,NBN
C
C           NUMERO DU POINT I
            NP = LFACES(I,NUFA)
C
C           MISE A JOUR DES COORDONNEES MAX OU MIN DES SURFACES
            DO K=1,3
               X = XYZPOI(K,NP)
               IF( X .LT. COOEXT(K,1) ) COOEXT(K,1)=X
               IF( X .GT. COOEXT(K,2) ) COOEXT(K,2)=X
            ENDDO
C
            DO NCAS=NCAS0,NCAS1
C
C              LA VALEUR DE LA SOLUTION
               IF( NTYP .EQ. 0 ) THEN
                  SOL = REAL( SOLUTION( NP, NCAS ) )
               ELSE
                  SOL = REAL( dptemp( NCAS )%dptab( NP ) )
               ENDIF
               IF( SOL .LT. SOLMIN ) SOLMIN = SOL
               IF( SOL .GT. SOLMAX ) SOLMAX = SOL
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
C     AFFICHAGE DES RESULTATS
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'COORDONNEES MIN et MAX des SURFACES NOMMEES:'
         WRITE(IMPRIM,*)'XMIN=',COOEXT(1,1),' XMAX=',COOEXT(1,2)
         WRITE(IMPRIM,*)'YMIN=',COOEXT(2,1),' YMAX=',COOEXT(2,2)
         WRITE(IMPRIM,*)'ZMIN=',COOEXT(3,1),' ZMAX=',COOEXT(3,2)
         WRITE(IMPRIM,*)'CAS',NCAS0,' au CAS',NCAS1,
     %   ' sur les SURFACES: Minimum=',SOLMIN,' Maximum=',SOLMAX
      ELSE
         WRITE(IMPRIM,*)'MIN MAX of NAMED SURFACES:'
         WRITE(IMPRIM,*)'XMIN=',COOEXT(1,1),' XMAX=',COOEXT(1,2)
         WRITE(IMPRIM,*)'YMIN=',COOEXT(2,1),' YMAX=',COOEXT(2,2)
         WRITE(IMPRIM,*)'ZMIN=',COOEXT(3,1),' ZMAX=',COOEXT(3,2)
         WRITE(IMPRIM,*)'CASE',NCAS0,' to CASE',NCAS1,
     %   ' on the SURFACES: Minimum=',SOLMIN,' Maximum=',SOLMAX
      ENDIF
C
      RETURN
      END
