      SUBROUTINE TRSOSF3( NBARFR, NBFAFR, NOFAFR,
     %                    MOFACE, MXFACE, LFACES,
     %                    MOARFR, MXARFR, LAREFR,
     %                    NBCOOR, XYZPOI,
     %                    NBCOMP, NCAS0,  NCAS1,
     %                    NTYP,   SOLUTION, dptemp, NCAS,
     %                    SOLMIN, SOLMAX  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES FRONTALIERES ET EN COULEURS
C ----     LA SOLUTION SUR LES FACES FRONTALIERES
C
C ENTREES:
C --------
C NBARFR : NOMBRE DES ARETES FRONTALIERES
C NBFAFR : NOMBRE DE FACES   FRONTALIERES
C NOFAFR : NUMERO DES FACES ET ARETES SELON LEUR DISTANCE DECROISSANTE A L'OEIL
C
C MOFACE : NOMBRE D'ENTIERS PAR FACE DANS LE TABLEAU LFACES
C MXFACE : MAJORATION DU NOMBRE DE FACES DU TABLEAU LFACES
C LFACES : TABLEAU DES FACES DU MAILLAGE
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C                       1 CUBE EST UN TETRA ou PENTA ou HEXAEDRE
C          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                       0 SI PAS DE 1-ER  CUBE
C          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                       ou CHAINAGE SUR LA FACE FRONTALIERE SUIVANTE
C          Cette valeur est remise a zero dans ce sous programme
C          (et sera regeneree dans un call macfaf.f)
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
C
C        ( LFACES(9,I)= NUMERO DU TYPE DU DERNIER EF DE CETTE FACE
C                       CETTE INFORMATION EXISTE SEULEMENT POUR UN OBJET
C                       PAS POUR UNE SURFACE OU UN VOLUME! )
C
C MOARFR : NOMBRE DE MOTS PAR ARETE FRONTALIERE DU TABLEAU LAREFR
C MXARFR : NOMBRE DE FACES DU TABLEAU LAREFR
C LAREFR : TABLEAU NUMERO DES 2 SOMMETS ET LIEN
C          LAREFR(1,I)= NO DU 1-ER  SOMMET DE L'ARETE FRONTALIERE
C          LAREFR(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LAREFR(3,I)= 0 OU NUMERO DE L'ARETE FRONTALIERE SUIVANTE
C                       DANS LE HACHAGE
C
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS (3 ou 6)
C XYZPOI : COORDONNEES DES SOMMETS DES FACES FRONTALIERES (ET INTERNES)
C
C NBCOMP : NOMBRE DE COMPOSANTES D'UN VECTEUR SOLUTION
C NCAS0  : NUMERO DU PREMIER VECTEUR SOLUTION A TRACER
C NCAS1  : NUMERO DU DERNIER VECTEUR SOLUTION A TRACER
C NTYP   : =0 EMPLOI de SOLUTION
C          =1 EMPLOI de dptemp
C SOLUTION: VECTEURS SOLUTION NCAS0 a NCAS1 de NBCOMP COMPOSANTES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NBCOMP) SOLUTIONS

C SOLUTION: VECTEURS SOLUTION NCAS0 a NCAS1 de NBCOMP COMPOSANTES
C NCAS   : NUMERO DU VECTEUR SOLUTION A TRACER
C SOLMIN : SOLUTION MINIMALE POUR LE CAS A TRACER
C SOLMAX : SOLUTION MAXIMALE POUR LE CAS A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & Saint PIERRE du PERRAY   AOUT 2010
C2345X7..............................................................012
      PARAMETER    ( LIGCON=0, LIGTIR=1 )
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      INTEGER           NOFAFR(1:NBFAFR+NBARFR)
      INTEGER           LFACES(1:MOFACE,1:MXFACE)
      INTEGER           LAREFR(1:MOARFR,1:MXARFR)
      REAL              XYZPOI(1:NBCOOR,1:*)
      DOUBLE PRECISION  SOLUTION(NBCOMP,NCAS0:NCAS1), S
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp
      REAL              COUL(4)
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
C     LE TRACE DES FACES ET ARETES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     =================================================================
      DO 50 NF = 1, NBFAFR + NBARFR
C
C        LE NUMERO DE LA FACE OU ARETE LA PLUS ELOIGNEE NON ENCORE TRACEE
         NAOF = NOFAFR( NF )
C
         IF( NAOF .GT. 0 ) THEN
C
C           TRACE DE LA SOLUTION SUR LA FACE FRONTIERE
C           ------------------------------------------
C           LE NOMBRE DE SOMMETS DE LA FACE
            IF( LFACES(4,NAOF) .EQ. 0 ) THEN
C              TRIANGLE
               NAF = 3
            ELSE
C              QUADRANGLE
               NAF = 4
            ENDIF
C
C           LA COULEUR DES NAF SOMMETS DE LA FACE FRONTALIERE SELON LA SOLUTION
            DO 10 I=1,NAF
C              LA COULEUR
               J = LFACES(I,NAOF)
               IF( NTYP .EQ. 0 ) THEN
                  S = SOLUTION( J, NCAS )
               ELSE
                  S = dptemp( NCAS )%dptab( J )
               ENDIF
               COUL(I) = REAL( ( S - SOLMIN)
     %                         / ( SOLMAX - SOLMIN ) * NBCOUL + N1COUL )
               IF( COUL(I) .GT. NDCOUL ) COUL(I)=NDCOUL
               IF( COUL(I) .LT. N1COUL ) COUL(I)=N1COUL
 10         CONTINUE
C
            CALL XVTYPETRAIT( NTLAPL )
C
C           PAS DE CALCUL DE REDUCTION DE LA FACE
            IF( NAF .EQ. 3 ) THEN
C
C              LE TRACE DE LA SOLUTION SUR LE TRIANGLE
               CALL TRIACOU3DBORD( XYZPOI(1,LFACES(1,NAOF)),
     %                             XYZPOI(1,LFACES(2,NAOF)),
     %                             XYZPOI(1,LFACES(3,NAOF)),
     %                             COUL(1), COUL(2), COUL(3),
     %                             NCOAPL, 1 )
            ELSE
C
C              LE TRACE DE LA SOLUTION SUR LE QUADRANGLE
               CALL QUADCOU3DBORD( XYZPOI(1,LFACES(1,NAOF)),
     %                             XYZPOI(1,LFACES(2,NAOF)),
     %                             XYZPOI(1,LFACES(3,NAOF)),
     %                             XYZPOI(1,LFACES(4,NAOF)),
     %                             COUL(1), COUL(2), COUL(3), COUL(4),
     %                             NCOAPL, 1 )
            ENDIF
C
         ELSE
C
C           TRACE DE L'ARETE FRONTALIERE DU MAILLAGE
C           ----------------------------------------
C           LE NUMERO DE L'ARETE DANS LE TABLEAU LAREFR
            IF( NCOAFR .GE. 0 ) THEN
               CALL XVTYPETRAIT( NTLAFR )
               NAOF = -NAOF
               CALL TRAIT3D( NCOAFR, XYZPOI(1,LAREFR(1,NAOF)),
     %                               XYZPOI(1,LAREFR(2,NAOF)) )
            ENDIF
         ENDIF
 50   CONTINUE
C
C     RETOUR AU TRACE CONTINU DES LIGNES
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
