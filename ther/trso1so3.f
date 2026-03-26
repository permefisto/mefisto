      SUBROUTINE TRSO1SO3( NBFACE, NOFACE,
     %                     MOFACE, MXFACE, LFACES,
     %                     NBCOOR, XYZPOI,
     %                     NBCOMP, NCAS0, NCAS1,
     %                     NTYP,   SOLUTION, dptemp, NCAS,
     %                     SOLMIN, SOLMAX  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA SOLUTION SUR LES FACES DES SURFACES NOMMEES
C -----    DE L'OBJET
C
C ENTREES:
C --------
C NBFACE : NOMBRE DE FACES DES SURFACES NOMMEES
C NOFACE : NUMERO DES FACES SELON LEUR DISTANCE DECROISSANTE A L'OEIL
C MOFACE : NOMBRE D'ENTIERS PAR FACE DANS LE TABLEAU LFACES
C MXFACE : MAJORATION DU NOMBRE DE FACES DU TABLEAU LFACES
C LFACES : TABLEAU NUMERO DES FACES DES SURFACES NOMMEES
C          LFACES(1...MOFACE,I)= NO DU POINT DE LA FACE I
C          LFACES(MOFACE,I)= 0 SI QUADRANGLE
C
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS (3 ou 6)
C XYZPOI : COORDONNEES DES NOEUDS DU MAILLAGE DE L'OBJET
C
C NBCOMP : NOMBRE DE COMPOSANTES D'UN VECTEUR SOLUTION
C NCAS0  : NUMERO MINIMUM DES VECTEURS SOLUTION A TRACER
C NCAS1  : NUMERO MAXIMUM DES VECTEURS SOLUTION A TRACER
C NTYP    : =0 EMPLOI de SOLUTION
C           =1 EMPLOI de dptemp
C SOLUTION: VECTEUR"SOLUTION(NBCOMP,NCAS0:NCAS1)
C dptemp  : TABLEAU DES NCAS0:NCAS1 TABLEAU(NBCOMP) SOLUTIONS
C SOLUTION:VECTEURS SOLUTION NCAS0 a NCAS1 de NBCOMP COMPOSANTES
C NCAS   : NUMERO DU VECTEUR SOLUTION A TRACER
C SOLMIN : SOLUTION MINIMALE POUR LE CAS A TRACER
C SOLMAX : SOLUTION MAXIMALE POUR LE CAS A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      PARAMETER    ( LIGCON=0 )
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      INTEGER           NOFACE(1:NBFACE)
      INTEGER           LFACES(1:MOFACE,1:MXFACE)
      REAL              XYZPOI(1:NBCOOR,1:*)
      DOUBLE PRECISION  SOLUTION(NBCOMP,NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp
      REAL              COUL(9)
C
C     LA PALETTE ARC EN CIEL
      CALL PALCDE(11)
C
C     LE NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL
C
C     LES ARETES SONT TRACEES AVEC UNE EPAISSEUR
      CALL XVEPAISSEUR( 1 )
C
C     TYPE DE TRACE DES ARETES
      CALL XVTYPETRAIT( NTLAFR )
C
C     LE TRACE DES FACES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     =======================================================
      DO 50 NF = 1, NBFACE
C
C        LE NUMERO DE LA FACE OU ARETE LA PLUS ELOIGNEE NON ENCORE TRACEE
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
C        LA COULEUR DES NBN POINTS DE LA FACE SELON LA SOLUTION
         DO I=1,NBN
C           LA COULEUR
            J = LFACES( I, NUFA )
            IF( NTYP .EQ. 0 ) THEN
               S = REAL( SOLUTION( J, NCAS ) )
            ELSE
               S = REAL( dptemp( NCAS )%dptab( J ) )
            ENDIF
            COUL(I) = ( S - SOLMIN )
     %              / ( SOLMAX - SOLMIN ) * NBCOUL + N1COUL
            IF( COUL(I) .GT. NDCOUL ) COUL(I)=NDCOUL
            IF( COUL(I) .LT. N1COUL ) COUL(I)=N1COUL
         ENDDO
C
         CALL XVTYPETRAIT( NTLAPL )
C
C        TRACE DE LA SOLUTION SUR LA FACE
         IF( NBN .EQ. 3 ) THEN
C
C           LE TRACE DE LA SOLUTION SUR LE TRIANGLE P1
            CALL TRIACOU3DBORD( XYZPOI(1,LFACES(1,NUFA)),
     %                          XYZPOI(1,LFACES(2,NUFA)),
     %                          XYZPOI(1,LFACES(3,NUFA)),
     %                          COUL(1), COUL(2), COUL(3),
     %                          NCOAFR, 1 )
         ELSE IF( NBN .EQ. 4 ) THEN
C
C           LE TRACE DE LA SOLUTION SUR LE QUADRANGLE Q1
            CALL QUADCOU3DBORD( XYZPOI(1,LFACES(1,NUFA)),
     %                          XYZPOI(1,LFACES(2,NUFA)),
     %                          XYZPOI(1,LFACES(3,NUFA)),
     %                          XYZPOI(1,LFACES(4,NUFA)),
     %                          COUL(1), COUL(2), COUL(3), COUL(4),
     %                          NCOAFR, 1 )
C
         ELSE IF( NBN .EQ. 6 ) THEN
C
C           LE TRACE DE LA SOLUTION SUR LE TRIANGLE P2
C           SOUS FORME DE 4 SOUS-TRIANGLES COLORES EN INTERPOLATION P1
C           TRIANGLE 1 4 6
            CALL TRIACOU3DBORD( XYZPOI(1,LFACES(1,NUFA)),
     %                          XYZPOI(1,LFACES(4,NUFA)),
     %                          XYZPOI(1,LFACES(6,NUFA)),
     %                          COUL(1), COUL(4), COUL(6),
     %                          NCOAFR, 1 )
C           TRIANGLE 4 5 6
            CALL TRIACOU3DBORD( XYZPOI(1,LFACES(4,NUFA)),
     %                          XYZPOI(1,LFACES(5,NUFA)),
     %                          XYZPOI(1,LFACES(6,NUFA)),
     %                          COUL(4), COUL(5), COUL(6),
     %                          NCOAFR, 1 )
C           TRIANGLE 4 2 5
            CALL TRIACOU3DBORD( XYZPOI(1,LFACES(4,NUFA)),
     %                          XYZPOI(1,LFACES(2,NUFA)),
     %                          XYZPOI(1,LFACES(5,NUFA)),
     %                          COUL(4), COUL(2), COUL(5),
     %                          NCOAFR, 1 )
C           TRIANGLE 6 5 3
            CALL TRIACOU3DBORD( XYZPOI(1,LFACES(6,NUFA)),
     %                          XYZPOI(1,LFACES(5,NUFA)),
     %                          XYZPOI(1,LFACES(3,NUFA)),
     %                          COUL(6), COUL(5), COUL(3),
     %                          NCOAFR, 1 )
C
         ELSE IF( NBN .EQ. 8 ) THEN
C
C           LE TRACE DE LA SOLUTION SUR LE QUADRANGLE Q2
C           SOUS FORME DE 4 SOUS-QUADRANGLES COLORES EN INTERPOLATION Q1
C
C           CALCUL DE LA COULEUR AU BARYCENTRE (QUADRANGLE Q2 REDUIT)
C           Cf Cours Thomas Raviart 1972  page X-22
C           P(Barycentre) = - SOM 1 a 4 P(Sommet) / 4 + Som 5 a 8 P(Milieu) / 2
            COUL(9) = (COUL(5)+COUL(6)+COUL(7)+COUL(8))*0.5
     %              - (COUL(1)+COUL(2)+COUL(3)+COUL(4))*0.25
C           QUADRANGLE 1 5 9 8
            CALL QUADCOU3DBORD( XYZPOI(1,LFACES(1,NUFA)),
     %                          XYZPOI(1,LFACES(5,NUFA)),
     %                          XYZPOI(1,LFACES(9,NUFA)),
     %                          XYZPOI(1,LFACES(8,NUFA)),
     %                          COUL(1), COUL(5), COUL(9), COUL(8),
     %                          NCOAFR, 1 )
C           QUADRANGLE 5 2 6 9
            CALL QUADCOU3DBORD( XYZPOI(1,LFACES(5,NUFA)),
     %                          XYZPOI(1,LFACES(2,NUFA)),
     %                          XYZPOI(1,LFACES(6,NUFA)),
     %                          XYZPOI(1,LFACES(9,NUFA)),
     %                          COUL(5), COUL(2), COUL(6), COUL(9),
     %                          NCOAFR, 1 )
C           QUADRANGLE 8 9 7 4
            CALL QUADCOU3DBORD( XYZPOI(1,LFACES(8,NUFA)),
     %                          XYZPOI(1,LFACES(9,NUFA)),
     %                          XYZPOI(1,LFACES(7,NUFA)),
     %                          XYZPOI(1,LFACES(4,NUFA)),
     %                          COUL(8), COUL(9), COUL(7), COUL(4),
     %                          NCOAFR, 1 )
C           QUADRANGLE 9 6 3 7
            CALL QUADCOU3DBORD( XYZPOI(1,LFACES(9,NUFA)),
     %                          XYZPOI(1,LFACES(6,NUFA)),
     %                          XYZPOI(1,LFACES(3,NUFA)),
     %                          XYZPOI(1,LFACES(7,NUFA)),
     %                          COUL(9), COUL(6), COUL(3), COUL(7),
     %                          NCOAFR, 1 )
         ENDIF
C
 50   CONTINUE
C
C     RETOUR AU TRACE CONTINU DES LIGNES
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
