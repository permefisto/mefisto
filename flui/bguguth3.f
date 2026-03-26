      SUBROUTINE BGUGUTH3( NBSOM,  NBNOVI, XYZNOE,
     %                     NBNOEF, NBEF,   NONOEF, NONOSO,
     %                     NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                     NCAS,   NTDLVP, NBVECT, VIXYZP, NDDL,
     %                     P2DP2,  CoTeNL,  BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL DU PROBLEME
C ----     -Coef LAPLACIEN P = Div( -Fomega + CoTeNL u. Grad u)
C          POUR LE TETRAEDRE TAYLOR-HOOD
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS VITESSE
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF (10 POUR TETRAEDRE TAYLOR-HOOD)
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF SOMMETS DES NBEF EF
C NONOSO : NONOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
C NOOBVC : NUMERO DE VOLUME DU FLUIDE
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DES VOLUMES
C NCAS   : NUMERO DU VECTEUR A TRAITER (COMPRIS ENTRE 1 et NBVECT)
C NTDLVP : NOMBRE TOTAL DE DL VITESSE+PRESSION
C NBVECT : NOMBRE DE VECTEURS VITESSE+PRESSION
C VIXYZP : 3 COMPOSANTES DE LA VITESSE + PRESSION AUX SOMMETS ET
C          3 COMPOSANTES DE LA VITESSE AU MILIEU DES ARETES
C          POUR CHACUN DES NBVECT (ou TEMPS) VECTEURS
C NDDL   : NDDL(I)=NO DU DERNIER DL DU NOEUD VITESSE I NDDL(0)=0
C P2DP2  : P2DP2(i,k,j) = Integrale P2i DP2j/Dxk dx
C CoTeNL : Coefficient du TERME NON LINEAIRE DE NAVIER-STOKES (TRANSPORT)
C
C SORTIE :
C --------
C BG     : VALEUR DES COEFFICIENTS DU VECTEUR GLOBAL SECOND MEMBRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     Juin 2011
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/donflu.inc"
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOVI, NBNOEF, NBEF,
     %                  NONOEF(NBEF,NBNOEF), NONOSO(NBNOVI),
     %                  NOOBVC, NUMIVO, NUMAVO
      DOUBLE PRECISION  CoTeNL, VIXYZP(NTDLVP,NBVECT), BG(NBSOM)
C
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO), NDDL(0:NBNOVI)
      DOUBLE PRECISION  DELTAe, DF(3,3), DFM1(3,3), DFM1DLa(3,4),
     %                  S, FL, X, Y, Z
      INTEGER           NOSOTE(4), NS1, NS2, NS3, NS4
      EQUIVALENCE      (NOSOTE(1),NS1),(NOSOTE(2),NS2),(NOSOTE(3),NS3),
     %                 (NOSOTE(4),NS4)
      DOUBLE PRECISION  FORCE(3,10)
      INTEGER           I, J, JJ, K, L, N, NSJ, NSJJ, NEF, NDJ,
     %                  NCAS, NBVECT, NTDLVP
      DOUBLE PRECISION  P2DP2(10,3,10), VJ, VJJ
C
C     INTEGRALE P2i dx SUR LE TETRAEDRE UNITE
      DOUBLE PRECISION  INTP2(10)
      DATA              INTP2/
     %  -0.8333333333333333333D-2,  -0.8333333333333333333D-2,
     %  -0.8333333333333333333D-2,  -0.8333333333333333333D-2,
     %   0.3333333333333333333D-1,   0.3333333333333333333D-1,
     %   0.3333333333333333333D-1,   0.3333333333333333333D-1,
     %   0.3333333333333333333D-1,   0.3333333333333333333D-1 /
C
C     MISE A ZERO DU SECOND MEMBRE GLOBAL
      DO I = 1, NBSOM
         BG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NUMERO DES 4 SOMMETS DU TETRAEDRE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
         NS4 = NONOEF(NEF,4)
C
C        CONSTRUCTION DE LA MATRICE DF
         X = XYZNOE(1,NS1)
         DF(1,1) = XYZNOE(1,NS2) - X
         DF(2,1) = XYZNOE(1,NS3) - X
         DF(3,1) = XYZNOE(1,NS4) - X
C
         Y = XYZNOE(2,NS1)
         DF(1,2) = XYZNOE(2,NS2) - Y
         DF(2,2) = XYZNOE(2,NS3) - Y
         DF(3,2) = XYZNOE(2,NS4) - Y
C
         Z = XYZNOE(3,NS1)
         DF(1,3) = XYZNOE(3,NS2) - Z
         DF(2,3) = XYZNOE(3,NS3) - Z
         DF(3,3) = XYZNOE(3,NS4) - Z
C
C        LE DETERMINANT DE DF
         DELTAe = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %          + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %          + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
C
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
         DFM1(1,1) = DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)
         DFM1(2,1) = DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1)
         DFM1(3,1) = DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2)
C
         DFM1(1,2) = DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3)
         DFM1(2,2) = DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1)
         DFM1(3,2) = DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2)
C
         DFM1(1,3) = DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)
         DFM1(2,3) = DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1)
         DFM1(3,3) = DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2)
C
C        [DFM1] [DLa]
         DFM1DLa(1,1) = -DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
         DFM1DLa(1,2) =  DFM1(1,1)
         DFM1DLa(1,3) =  DFM1(1,2)
         DFM1DLa(1,4) =  DFM1(1,3)
C
         DFM1DLa(2,1) = -DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
         DFM1DLa(2,2) =  DFM1(2,1)
         DFM1DLa(2,3) =  DFM1(2,2)
         DFM1DLa(2,4) =  DFM1(2,3)
C
         DFM1DLa(3,1) = -DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
         DFM1DLa(3,2) =  DFM1(3,1)
         DFM1DLa(3,3) =  DFM1(3,2)
         DFM1DLa(3,4) =  DFM1(3,3)
C
C        CONTRIBUTION DES EFFORTS VOLUMIQUES INTERPOLES TAYLOR-HOOD
C        ----------------------------------------------------------
         IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
C
C           VALEUR DES EFFORTS VOLUMIQUES AUX 10 NOEUDS DU TETRAEDRE
            DO J=1,10
               NSJ = NONOEF(NEF,J)
               X = XYZNOE(1,NSJ)
               Y = XYZNOE(2,NSJ)
               Z = XYZNOE(3,NSJ)
               CALL REFORC( 4,NOOBVC, 3, X,Y,Z,  0D0,0D0,0D0,
     %                      LTDEVO(LPFORC,NOOBVC), FORCE(1,J) )
            ENDDO
C
         ENDIF
C
         DO I=1,4
C
C           I-eme COEFFICIENT DU VECTEUR ELEMENTAIRE BE
            S = 0D0
C
            DO L=1,3
C
C              COEFFICIENTS DES EFFORTS VOLUMIQUES L
               FL = 0D0
               IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
                  DO J=1,10
                     FL = FL + INTP2(J) * FORCE(L,J)
                  ENDDO
               ENDIF
C
C              - Integrale  Grad P1  CoTeNL ( u. Grad u) dx
               DO K=1,3
                  DO J=1,10
C
C                    NO GLOBAL DU NOEUD J DU TETRAEDRE TH
                     NSJ = NONOEF(NEF,J)
C                    NO DU DERNIER DL DU NOEUD PRECEDENT NSJ
                     NDJ = NDDL( NSJ-1 )
C                    COMPOSANTE K DE LA VITESSE AU NOEUD NSJ
                     VJ  = VIXYZP(NDJ+K,NCAS)
C
                     DO N=1,3
                        DO JJ=1,10
C
C                          NO GLOBAL DU NOEUD JJ DU TETRAEDRE TH
                           NSJJ = NONOEF(NEF,JJ)
C                          NO DU DERNIER DL DU NOEUD PRECEDENT NSJJ
                           NDJ = NDDL( NSJJ-1 )
C                          COMPOSANTE L DE LA VITESSE AU NOEUD NSJJ
                           VJJ = VIXYZP(NDJ+L,NCAS)
C                          P2DP2(i,k,j) = integrale P2j dP2jj/dxn dx dy dz
                           FL = FL
     %                        - CoTeNL * P2DP2(J,N,JJ) * VJ
     %                                 * DFM1(K,N) * VJJ / DELTAe
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               S = S + DFM1DLa(L,I) * FL
            ENDDO
C
C           ASSEMBLAGE DE BE(I) DANS BG( NONOSO( NOSOTE(I) ) )
            NSJ = NONOSO( NOSOTE(I) )
            BG(NSJ) = BG(NSJ) + S
C
         ENDDO
C
 100  CONTINUE
C
      RETURN
      END
