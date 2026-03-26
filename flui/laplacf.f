      SUBROUTINE LAPLACF( NBSOM,  XYZSOM, NBEF, NONOSO,
     %                    MOARET, MXARET, LARETE,
     %                    NBNOVI, NCAS0,  NCAS1, NCAS, vitx, vity,
     %                    BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE ROT VITESSE
C ----     c-a-d Integrale de tLambda (-v nx + u ny ) dGamma
C          SUR LES ARETES DE LA FRONTIERE
C          LES VITESSES SONT ICI INTERPOLEES P1 SUR CHAQUE ARETE
C          LE DEGRE DE LIBERTE DU MILIEU D'ARETE POUR TAYLOR-HOOD
C          N'EST PAS PRIS EN COMPTE ICI
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOSO : NONOSO(I)=NO DU SOMMET DE 1 A NBSOM DU NOEUD I
C MOARET : NOMBRE DE MOTS PAR ARETE CHAINEE DU TABLEAU LARETE
C MXARET : NOMBRE MAXIMAL D'ARETES RANGEES DANS LARETE
C L1ARCH : NUMERO DE LA PREMIERE ARETE CHAINEE DANS LE TABLEAU LARETE
C LARETE : TABLEAU NUMERO DES 2 SOMMETS, LIEN, EF1, EF2, CHAINAGE
C          LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE CHAINEE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= 0 OU NUMERO DE L'ARETE CHAINEE SUIVANTE
C          LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 1-ER  TRIANGLE
C          LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 2-EME TRIANGLE
C          LARETE(6,I)= NUMERO DANS LARETE DE L'ARETE SUIVANTE
C            SOIT DANS LE CHAINAGE D'UNE LIGNE J ENTRE NUMILF ET NUMXLF
C            SOIT DANS LE CHAINAGE DES ARETES CHAINEES
C            0 SI C'EST LA DERNIERE
C NBNOVI : NOMBRE DE NOEUDS VITESSE DU MAILLAGE
C NCAS0:NCAS1 : LES VECTEURS VITESSE STOCKES
C NCAS   : NUMERO DU CAS A TRAITE
C vitx   : COMPOSANTE X des VECTEURS VITESSE NCAS0:NCAS1
C vitx   : COMPOSANTE Y des VECTEURS VITESSE NCAS0:NCAS1

C MODIFIES:
C ---------
C BG     : COEFFICIENTS DU SECOND MEMBRE GLOBAL ASSEMBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Janvier 2011
C MODIFS : ALAIN PERRONNET             St PIERRE du PERRAY  Mai     2021
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NONOSO(NBNOVI)
      INTEGER           LARETE(1:MOARET,1:MXARET)

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(NCAS0:NCAS1) :: vitx, vity

      DOUBLE PRECISION  BG(NBSOM), TX, TY, S, SENS
      INTEGER           NBSOM, NBEF, MOARET, MXARET, NBNOVI,
     %                  NCAS0, NCAS1, NCAS, NS1, NS2, NEF, N
ccc      INTEGER           NBARFR
C
ccc      NBARFR = 0
      DO 100 N = 1, MXARET
C
         IF( LARETE(1,N) .NE. 0 ) THEN
C           L'ARETE N EXISTE
C
            IF( LARETE(5,N) .EQ. 0 ) THEN
C
C              L'ARETE N EST FRONTALIERE, APPARTENANT A UN SEUL EF
C              LA NORMALE DOIT ETRE EXTERIEURE A L'EF NEF
C              IL FAUT DONC RETROUVER LE BON SENS
               NEF = LARETE(4,N)
               IF( NEF .EQ. 0 ) GOTO 100
               IF( NEF .LT. 0 ) THEN
C                 ARETE INDIRECTE
                  SENS = -0.25D0
               ELSE
C                 ARETE DIRECTE AVEC NORMALE EXTERIEURE A DROITE
                  SENS = 0.25D0
               ENDIF
C
C              CONTRIBUTION DE L'ARETE NS1-NS2 AU SECOND MEMBRE
ccc               NBARFR = NBARFR + 1
               NS1 = LARETE(1,N)
               NS2 = LARETE(2,N)
C
C              COMPOSANTES DE L'ARETE NS1-NS2
               TX = XYZSOM(1,NS1) - XYZSOM(1,NS2)
               TY = XYZSOM(2,NS1) - XYZSOM(2,NS2)
C
ccc               S = SENS * ( (VITESX(NS1,NCAS) + VITESX(NS2,NCAS)) * TX
ccc     %                     +(VITESY(NS1,NCAS) + VITESY(NS2,NCAS)) * TY )

               S = SENS * ( ( vitx(NCAS)%dptab(NS1)
     %                      + vitx(NCAS)%dptab(NS2) ) * TX
     %                    + ( vity(NCAS)%dptab(NS1)
     %                      + vity(NCAS)%dptab(NS2) ) * TY )
C
               IF( NBSOM+NBEF .NE. NBNOVI ) THEN
C                 EF TAYLOR-HOOD
C                 PASSAGE DU NUMERO DE NOEUD AU NUMERO DE SOMMET
                  NS1 = NONOSO(NS1)
                  NS2 = NONOSO(NS2)
               ENDIF
C
C              ASSEMBLAGE DANS LE VECTEUR GLOBAL
               BG(NS1) = BG(NS1) + S
               BG(NS2) = BG(NS2) + S
C
            ENDIF
C
         ENDIF
C
 100  CONTINUE
C
ccc      print *,NBARFR,' ARETES FRONTALIERES TRAITEES DANS LAPLACF'
      RETURN
      END
