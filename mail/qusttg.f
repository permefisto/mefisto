      SUBROUTINE QUSTTG( NBS1,   NBS2,   NBA1,   NBA2,
     S                   MNXYT1, MNNTG1, MNXYT2, MNNTG2,
     S                   MNXYT3, MNNTG3, MNXYT4, MNNTG4,
     S                   COSO,   COTG,   STCARR, TGPURE  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES DES SOMMETS INTERNES DU
C -----     QUADRANGLE A MAILLER PAR LA METHODE DE COOK-GORDON-HALL
C           EN PRESENCE DE TANGENTES
C
C ENTREES :
C ---------
C NBS1    : NOMBRE DE POINTS SUR LES LIGNES "HORIZONTALES"
C NBS2    : NOMBRE DE POINTS SUR LES LIGNES "VERTICALES"
C MNXYT1234: ADRESSE MCN DE LA PREMIERE COMPOSANTE DES TANGENTES DU COTE 1234
C MNNTG1234: ADRESSE MCN DES 2 NUMEROS DES TANGENTES DES ARETES  DU COTE 1234
C STCARR  : COORDONNEES DES NBS1 x NBS2 SOMMETS DU CARRE UNITE
C COSO    : COORDONNEES DES SOMMETS DES 4 COTES (BORD) DU MAILLAGE
C
C SORTIES :
C ---------
C COSO    : 3 COORDONNEES DES NBS1 x NBS2 SOMMETS DU MAILLAGE
C COTG    : 3 COMPOSANTES DES 2 TANGENTES AUX 4 SOMMETS
C           DES NBA1 x NBA2 QUADRANGLES DU QUADRANGLE A MAILLER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1996
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              STCARR(2,NBS1,NBS2)
      REAL              COSO(3,NBS1,NBS2)
      REAL              COTG(3,8,NBA1,NBA2)
      REAL              TGPURE(3,2,NBS1,NBS2)
      REAL              TG2AR(3,2), FB(4), DFB(4), F(3,4), DF(3,4)
      REAL              ST4C(2,0:1,0:1)
      DATA              ST4C / 0.,0., 1.,0., 0.,1., 1.,1. /
C
C     CALCUL DES Ds/Du ET Ds/Dv AUX SOMMETS INTERNES DU CARRE UNITE
C     CALCUL DES 3 COORDONNEES DES SOMMETS INTERNES DU QUADRANGLE A MAILLER
C     =====================================================================
C     BOUCLE SUR LES POINTS INTERNES DU DOMAINE
      DO 100 J = 2, NBS2-1
C
         DO 90 I = 2, NBS1-1
C
C           LES 2 COORDONNEES DU SOMMET I,J DU CARRE UNITE
            XM = STCARR(1,I,J)
            YM = STCARR(2,I,J)
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 1
C           ET DE LA DERIVEE TANGENTIELLE EN CE POINT
C           -----------------------------------------
            N2=1
            DO 10 I1=2,NBS1
               IF( XM .LE. STCARR(1,I1,1)+EPSXYZ ) THEN
                  N2=I1
                  GO TO 11
               END IF
 10         CONTINUE
C           LES NUMEROS DES SOMMETS DE L'ARETE CONTENANT LE POINT PROJETE
 11         N1 = N2-1
            S1 = (    XM           - STCARR(1,N1,1) )
     S          / ( STCARR(1,N2,1) - STCARR(1,N1,1) )
C           RECHERCHE DES 2 TANGENTES DE L'ARETE N1
C           ARETE P3 HERMITE => C1
            CALL TG2ARL( N1, COSO(1,N1,1), COSO(1,N2,1),
     %                   MNXYT1, MNNTG1, TG2AR )
C           LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
            CALL VFBP3H( S1, FB )
C           LA VALEUR DE LA DERIVEE DESS POLYNOMES DE BASE P3 HERMITE AU POINT S
            CALL DFBP3H( S1, DFB )
            DO 13 K=1,3
C              LA COORDONNEE K DU POINT SUR LE COTE 1
               F(K,1) = FB(1) * COSO(K,N1,1)
     %                + FB(2) * COSO(K,N2,1)
     %                + FB(3) * TG2AR(K,1)
     %                + FB(4) * TG2AR(K,2)
C              LA COMPOSANTE K DE LA DERIVEE TANGENTIELLE AU POINT SUR L'ARETE 1
               DF(K,1) = ( DFB(1) * COSO(K,N1,1)
     %                   + DFB(2) * COSO(K,N2,1)
     %                   + DFB(3) * TG2AR(K,1)
     %                   + DFB(4) * TG2AR(K,2) )
     %                 / ( STCARR(1,N2,1) - STCARR(1,N1,1) )
 13         CONTINUE
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 2
C           ET DE LA DERIVEE TANGENTIELLE EN CE POINT
C           -----------------------------------------
            N2=1
            DO 20 I1=2,NBS2
               IF( YM .LE. STCARR(2,NBS1,I1)+EPSXYZ ) THEN
                  N2=I1
                  GO TO 21
               END IF
 20         CONTINUE
C           LES NUMEROS DES SOMMETS VOISINS
 21         N1 = N2-1
            S1 = (   YM               - STCARR(2,NBS1,N1) )
     S          / ( STCARR(2,NBS1,N2) - STCARR(2,NBS1,N1) )
C           RECHERCHE DES 2 TANGENTES DE L'ARETE IM2
C           ARETE P3 HERMITE => C1
            CALL TG2ARL( N1, COSO(1,NBS1,N1), COSO(1,NBS1,N2),
     %                   MNXYT2, MNNTG2, TG2AR )
C           LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
            CALL VFBP3H( S1, FB )
C           LA VALEUR DE LA DERIVEE DESS POLYNOMES DE BASE P3 HERMITE AU POINT S
            CALL DFBP3H( S1, DFB )
            DO 23 K=1,3
C              LA COORDONNEE K DU POINT SUR LE COTE 2
               F(K,2) = FB(1) * COSO(K,NBS1,N1)
     %                + FB(2) * COSO(K,NBS1,N2)
     %                + FB(3) * TG2AR(K,1)
     %                + FB(4) * TG2AR(K,2)
C              LA COMPOSANTE K DE LA DERIVEE TANGENTIELLE AU POINT SUR L'ARETE 1
               DF(K,2) =( DFB(1) * COSO(K,NBS1,N1)
     %                  + DFB(2) * COSO(K,NBS1,N2)
     %                  + DFB(3) * TG2AR(K,1)
     %                  + DFB(4) * TG2AR(K,2) )
     %                / ( STCARR(2,NBS1,N2) - STCARR(2,NBS1,N1) )
 23         CONTINUE
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 3
C           ET DE LA DERIVEE TANGENTIELLE EN CE POINT
C           -----------------------------------------
            N2 = 1
            DO 31 I1=2,NBS1
               IF( XM .LE. STCARR(1,I1,NBS2)+EPSXYZ ) THEN
                  N2 = I1
                  GO TO 32
               END IF
 31         CONTINUE
C           LES NUMEROS DES SOMMETS VOISINS
 32         N1 = N2-1
            S1 = (   XM               - STCARR(1,N1,NBS2) )
     S          / ( STCARR(1,N2,NBS2) - STCARR(1,N1,NBS2) )
C           RECHERCHE DES 2 TANGENTES DE L'ARETE N1
C           ARETE P3 HERMITE => C1
            CALL TG2ARL( N1, COSO(1,N1,NBS2), COSO(1,N2,NBS2),
     %                   MNXYT3, MNNTG3, TG2AR )
C           LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
            CALL VFBP3H( S1, FB )
C           LA VALEUR DE LA DERIVEE DESS POLYNOMES DE BASE P3 HERMITE AU POINT S
            CALL DFBP3H( S1, DFB )
            DO 33 K=1,3
C              LA COORDONNEE K DU POINT SUR LE COTE 3
               F(K,3) = FB(1) * COSO(K,N1,NBS2)
     %                + FB(2) * COSO(K,N2,NBS2)
     %                + FB(3) * TG2AR(K,1)
     %                + FB(4) * TG2AR(K,2)
C              LA COMPOSANTE K DE LA DERIVEE TANGENTIELLE AU POINT SUR L'ARETE 1
               DF(K,3) =( DFB(1) * COSO(K,N1,NBS2)
     %                  + DFB(2) * COSO(K,N2,NBS2)
     %                  + DFB(3) * TG2AR(K,1)
     %                  + DFB(4) * TG2AR(K,2) )
     %                / ( STCARR(1,N2,NBS2) - STCARR(1,N1,NBS2) )
 33         CONTINUE
C
C           RECHERCHE DU POINT PROJETE SUR LA LIGNE 4
C           ET DE LA DERIVEE TANGENTIELLE EN CE POINT
C           -----------------------------------------
            N2 = 1
            DO 41 I1=2,NBS2
               IF( YM .LE. STCARR(2,1,I1)+EPSXYZ ) THEN
                  N2=I1
                  GO TO 42
               END IF
 41         CONTINUE
C           LES NUMEROS DES SOMMETS
 42         N1 = N2 - 1
            S1 = (   YM           - STCARR(2,1,N1) )
     S         / ( STCARR(2,1,N2) - STCARR(2,1,N1) )
C           RECHERCHE DES 2 TANGENTES DE L'ARETE N1
C           ARETE P3 HERMITE => C1
            CALL TG2ARL( N1, COSO(1,1,N1), COSO(1,1,N2),
     %                   MNXYT4, MNNTG4, TG2AR )
C           LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
            CALL VFBP3H( S1, FB )
C           LA VALEUR DE LA DERIVEE DESS POLYNOMES DE BASE P3 HERMITE AU POINT S
            CALL DFBP3H( S1, DFB )
            DO 43 K=1,3
C              LA COORDONNEE K DU POINT SUR LE COTE 4
               F(K,4) = FB(1) * COSO(K,1,N1)
     %                + FB(2) * COSO(K,1,N2)
     %                + FB(3) * TG2AR(K,1)
     %                + FB(4) * TG2AR(K,2)
C              LA COMPOSANTE K DE LA DERIVEE TANGENTIELLE AU POINT SUR L'ARETE 1
               DF(K,4) = ( DFB(1) * COSO(K,1,N1)
     %                   + DFB(2) * COSO(K,1,N2)
     %                   + DFB(3) * TG2AR(K,1)
     %                   + DFB(4) * TG2AR(K,2) )
     %                 / ( STCARR(2,1,N2) - STCARR(2,1,N1) )
 43         CONTINUE
C
C           LE POINT (I,J) DU MAILLAGE ET SES 2 DERIVEES
C           D'APRES LA FORMULE DE COOK / GORDON ET HALL
C           --------------------------------------------
            XM1 = 1 - XM
            YM1 = 1 - YM
            DO 50 K=1,3
C              s(u,v) SUR LE CARRE UNITE (u,v)
               COSO(K,I,J) = YM1 * ( F(K,1) - XM1 * COSO(K,   1,   1) )
     S                     + XM  * ( F(K,2) - YM1 * COSO(K,NBS1,   1) )
     S                     + YM  * ( F(K,3) - XM  * COSO(K,NBS1,NBS2) )
     S                     + XM1 * ( F(K,4) - YM  * COSO(K,   1,NBS2) )
 50         CONTINUE
C
            DO 60 K=1,3
C              DERIVEE s/ Du SUR LE CARRE UNITE (u,v)
               TGPURE(K,1,I,J) = F(K,2) - F(K,4)
     S          + YM  * (DF(K,3)-COSO(K,NBS1,NBS2)+COSO(K,1,NBS2))
     S          + YM1 * (DF(K,1)-COSO(K,NBS1,   1)+COSO(K,1,   1))
C              DERIVEE s/ Dv SUR LE CARRE UNITE (u,v)
               TGPURE(K,2,I,J) = F(K,3) - F(K,1)
     S          + XM  * (DF(K,2)-COSO(K,NBS1,NBS2)+COSO(K,NBS1,1))
     S          + XM1 * (DF(K,4)-COSO(K,   1,NBS2)+COSO(K,   1,1))
 60         CONTINUE
 90      CONTINUE
 100  CONTINUE
C
C
C     LA RECHERCHE DES 2 TANGENTES DES 4 SOMMETS DU CARRE UNITE
C     ---------------------------------------------------------
C     SOMMET 1 : DS/DU ARETE P3 HERMITE => C1
      CALL TG2ARL( 1, COSO(1,1,1), COSO(1,2,1),
     S             MNXYT1, MNNTG1, TG2AR )
C     LA LONGUEUR DE L'ARETE 1 DU COTE 1 DU CARRE
      S1 = STCARR(1,2,1) - STCARR(1,1,1)
      TGPURE(1,1,1,1) = TG2AR(1,1) / S1
      TGPURE(2,1,1,1) = TG2AR(2,1) / S1
      TGPURE(3,1,1,1) = TG2AR(3,1) / S1
C
C     SOMMET 1 : DS/DV ARETE P3 HERMITE => C1
      CALL TG2ARL( 1, COSO(1,1,1), COSO(1,1,2),
     S             MNXYT4, MNNTG4, TG2AR )
C     LA LONGUEUR DE L'ARETE 1 DU COTE 4 DU CARRE
      S1 = STCARR(2,1,2) - STCARR(2,1,1)
      TGPURE(1,2,1,1) = TG2AR(1,1) / S1
      TGPURE(2,2,1,1) = TG2AR(2,1) / S1
      TGPURE(3,2,1,1) = TG2AR(3,1) / S1
C
C     SOMMET 2 : DS/DU  ARETE P3 HERMITE => C1
      CALL TG2ARL( NBA1, COSO(1,NBA1,1), COSO(1,NBS1,1),
     S             MNXYT1, MNNTG1, TG2AR )
C     LA LONGUEUR DE L'ARETE NBA1 DU COTE 1 DU CARRE
      S1 = STCARR(1,NBS1,1) - STCARR(1,NBA1,1)
      TGPURE(1,1,NBS1,1) = -TG2AR(1,2) / S1
      TGPURE(2,1,NBS1,1) = -TG2AR(2,2) / S1
      TGPURE(3,1,NBS1,1) = -TG2AR(3,2) / S1
C
C     SOMMET 2 : DS/DV  ARETE P3 HERMITE => C1
      CALL TG2ARL( 1, COSO(1,NBS1,1), COSO(1,NBS1,2),
     S             MNXYT2, MNNTG2, TG2AR )
C     LA LONGUEUR DE L'ARETE 1 DU COTE 2 DU CARRE
      S1 = STCARR(2,NBS1,2) - STCARR(2,NBS1,1)
      TGPURE(1,2,NBS1,1) = TG2AR(1,1) / S1
      TGPURE(2,2,NBS1,1) = TG2AR(2,1) / S1
      TGPURE(3,2,NBS1,1) = TG2AR(3,1) / S1
C
C     SOMMET 3: DS/DU ARETE P3 HERMITE => C1
      CALL TG2ARL( NBA1, COSO(1,NBA1,NBS2), COSO(1,NBS1,NBS2),
     S             MNXYT3, MNNTG3, TG2AR )
C     LA LONGUEUR DE L'ARETE NBA1 DU COTE 3 DU CARRE
      S1 = STCARR(1,NBS1,NBS2) - STCARR(1,NBA1,NBS2)
      TGPURE(1,1,NBS1,NBS2) = -TG2AR(1,2) / S1
      TGPURE(2,1,NBS1,NBS2) = -TG2AR(2,2) / S1
      TGPURE(3,1,NBS1,NBS2) = -TG2AR(3,2) / S1
C
C     SOMMET 3 : DS/DV ARETE P3 HERMITE => C1
      CALL TG2ARL( NBA2, COSO(1,NBS1,NBA2), COSO(1,NBS1,NBS2),
     S             MNXYT2, MNNTG2, TG2AR )
C     LA LONGUEUR DE L'ARETE DU COTE 2 DU CARRE
      S1 = STCARR(2,NBS1,NBS2) - STCARR(2,NBS1,NBA2)
      TGPURE(1,2,NBS1,NBS2) = -TG2AR(1,2) / S1
      TGPURE(2,2,NBS1,NBS2) = -TG2AR(2,2) / S1
      TGPURE(3,2,NBS1,NBS2) = -TG2AR(3,2) / S1
C
C     SOMMET 4 : DS/DU  ARETE P3 HERMITE => C1
      CALL TG2ARL( 1, COSO(1,1,NBS2), COSO(1,2,NBS2),
     S             MNXYT3, MNNTG3, TG2AR )
C     LA LONGUEUR DE L'ARETE 1 DU COTE 3 DU CARRE
      S1 = STCARR(1,2,NBS2) - STCARR(1,1,NBS2)
      TGPURE(1,1,1,NBS2) = TG2AR(1,1) / S1
      TGPURE(2,1,1,NBS2) = TG2AR(2,1) / S1
      TGPURE(3,1,1,NBS2) = TG2AR(3,1) / S1
C
C     SOMMET 4 : DS/DV  ARETE P3 HERMITE => C1
      CALL TG2ARL( NBA2, COSO(1,1,NBA2), COSO(1,1,NBS2),
     S             MNXYT4, MNNTG4, TG2AR )
C     LA LONGUEUR DE L'ARETE NBA2 DU COTE 4 DU CARRE
      S1 = STCARR(2,1,NBS2) - STCARR(2,1,NBA2)
      TGPURE(1,2,1,NBS2) = -TG2AR(1,2) / S1
      TGPURE(2,2,1,NBS2) = -TG2AR(2,2) / S1
      TGPURE(3,2,1,NBS2) = -TG2AR(3,2) / S1
C
C     LA RECHERCHE DES 2 TANGENTES DES SOMMETS DU COTE 1 DU CARRE UNITE
C     -----------------------------------------------------------------
      DO 110 I=2,NBA1
C
C        Ds/Du AU POINT I DU COTE 1 DU CARRE
         CALL TG2ARL( I, COSO(1,I,1), COSO(1,I+1,1),
     S                MNXYT1, MNNTG1, TG2AR )
C        LA LONGUEUR DE L'ARETE I DU COTE 1 DU CARRE
         S1 = STCARR(1,I+1,1) - STCARR(1,I,1)
         TGPURE(1,1,I,1) = TG2AR(1,1) / S1
         TGPURE(2,1,I,1) = TG2AR(2,1) / S1
         TGPURE(3,1,I,1) = TG2AR(3,1) / S1
C
C        Ds/Dv AU POINT I DU COTE 1 DU CARRE PAR LA FORMULE DE GORDON HALL
C        L'INTERVALLE DE XM SUR LE COTE 3 DU CARRE
         XM  = STCARR(1,I,1)
         XM1 = 1.0 - XM
         N2 = 1
         DO 103 I1=2,NBS1
            IF( XM .LE. STCARR(1,I1,NBS2)+EPSXYZ ) THEN
               N2 = I1
               GO TO 104
            END IF
 103     CONTINUE
C        LES NUMEROS DES SOMMETS VOISINS
 104     N1 = N2-1
C        LES 2 TANGENTES DE L'ARETE N1 DU COTE 3
         CALL TG2ARL( N1, COSO(1,N1,NBS2), COSO(1,N2,NBS2),
     S                MNXYT3, MNNTG3, TG2AR )
C        LE PARAMETRE ENTRE 0 ET 1
         S1 = (   XM               - STCARR(1,N1,NBS2) )
     S       / ( STCARR(1,N2,NBS2) - STCARR(1,N1,NBS2) )
C        LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
         CALL VFBP3H( S1, FB )
         DO 105 K=1,3
C           LA COORDONNEE K DU POINT SUR L'ARETE 3
            F(K,3) = FB(1) * COSO(K,N1,NBS2)
     %             + FB(2) * COSO(K,N2,NBS2)
     %             + FB(3) * TG2AR(K,1)
     %             + FB(4) * TG2AR(K,2)
 105     CONTINUE
         DO 107 K=1,3
            TGPURE(K,2,I,1) = F(K,3) - COSO(K,I,1)
     S    + XM *(TGPURE(K,2,NBS1,1)-COSO(K,NBS1,NBS2)+COSO(K,NBS1,1))
     S    + XM1*(TGPURE(K,2,1,1)-COSO(K,1,NBS2)+COSO(K,1,1))
 107     CONTINUE
 110  CONTINUE
C
C     LA RECHERCHE DES 2 TANGENTES DES SOMMETS DU COTE 2 DU CARRE UNITE
C     -----------------------------------------------------------------
      DO 120 I=2,NBA2
C
C        Ds/Dv AU POINT I DU COTE 2 DU CARRE
         CALL TG2ARL( I, COSO(1,NBS1,I), COSO(1,NBS1,I+1),
     S                MNXYT2, MNNTG2, TG2AR )
C        LA LONGUEUR DE L'ARETE I DU COTE 2 DU CARRE
         S1 = STCARR(2,NBS1,I+1) - STCARR(2,NBS1,I)
         TGPURE(1,2,NBS1,I) = TG2AR(1,1) / S1
         TGPURE(2,2,NBS1,I) = TG2AR(2,1) / S1
         TGPURE(3,2,NBS1,I) = TG2AR(3,1) / S1
C
C        Ds/Du AU POINT I DU COTE 2 DU CARRE PAR LA FORMULE DE GORDON HALL
C        L'INTERVALLE DE YM SUR LE COTE 4 DU CARRE
         YM  = STCARR(2,NBS1,I)
         YM1 = 1.0 - YM
         N2  = 1
C        INTERVALLE SUR LE COTE 4 CONTENANT YM
         DO 113 I1=2,NBS2
            IF( YM .LE. STCARR(2,1,I1)+EPSXYZ ) THEN
               N2 = I1
               GO TO 114
            END IF
 113     CONTINUE
C        LES NUMEROS DES SOMMETS VOISINS
 114     N1 = N2-1
C        LES 2 TANGENTES DE L'ARETE N1 DU COTE 4
         CALL TG2ARL( N1, COSO(1,1,N1), COSO(1,1,N2),
     S                MNXYT4, MNNTG4, TG2AR )
C        LE PARAMETRE ENTRE 0 ET 1
         S1  = (   YM            - STCARR(2,1,N1) )
     S        / ( STCARR(2,1,N2) - STCARR(2,1,N1) )
C        LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
         CALL VFBP3H( S1, FB )
         DO 115 K=1,3
C           LA COORDONNEE K DU POINT SUR L'ARETE 4
            F(K,4) = FB(1) * COSO(K,1,N1)
     %             + FB(2) * COSO(K,1,N2)
     %             + FB(3) * TG2AR(K,1)
     %             + FB(4) * TG2AR(K,2)
 115     CONTINUE
         DO 117 K=1,3
            TGPURE(K,1,NBS1,I) = COSO(K,NBS1,I) - F(K,4)
     S    + YM *(TGPURE(K,1,NBS1,NBS2)-COSO(K,NBS1,NBS2)+COSO(K,1,NBS2))
     S    + YM1*(TGPURE(K,1,NBS1,1)-COSO(K,NBS1,1)+COSO(K,1,1))
 117     CONTINUE
 120  CONTINUE
C
C     LA RECHERCHE DES 2 TANGENTES DES SOMMETS DU COTE 3 DU CARRE UNITE
C     -----------------------------------------------------------------
      DO 130 I=2,NBA1
C
C        Ds/Du AU POINT I DU COTE 3 DU CARRE
         CALL TG2ARL( I, COSO(1,I,NBS2), COSO(1,I+1,NBS2),
     S                MNXYT3, MNNTG3, TG2AR )
C        LA LONGUEUR DE L'ARETE I DU COTE 3 DU CARRE
         S1 = STCARR(1,I+1,NBS2) - STCARR(1,I,NBS2)
         TGPURE(1,1,I,NBS2) = TG2AR(1,1) / S1
         TGPURE(2,1,I,NBS2) = TG2AR(2,1) / S1
         TGPURE(3,1,I,NBS2) = TG2AR(3,1) / S1
C
C        Ds/Dv AU POINT I DU COTE 3 DU CARRE PAR LA FORMULE DE GORDON HALL
C        L'INTERVALLE DE XM SUR LE COTE 1 DU CARRE
         XM  = STCARR(1,I,NBS2)
         XM1 = 1.0 - XM
         N2  = 1
         DO 123 I1=2,NBS1
            IF( XM .LE. STCARR(1,I1,1)+EPSXYZ ) THEN
               N2 = I1
               GO TO 124
            END IF
 123     CONTINUE
C        LES NUMEROS DES SOMMETS VOISINS
 124     N1 = N2-1
C        LES 2 TANGENTES DE L'ARETE N1 DU COTE 1
         CALL TG2ARL( N1, COSO(1,N1,1), COSO(1,N2,1),
     S                MNXYT1, MNNTG1, TG2AR )
C        LE PARAMETRE ENTRE 0 ET 1
         S1 = (   XM            - STCARR(1,N1,1) )
     S       / ( STCARR(1,N2,1) - STCARR(1,N1,1) )
C        LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
         CALL VFBP3H( S1, FB )
         DO 125 K=1,3
C           LA COORDONNEE K DU POINT SUR LE COTE 1 DU CARRE
            F(K,1) = FB(1) * COSO(K,N1,1)
     %             + FB(2) * COSO(K,N2,1)
     %             + FB(3) * TG2AR(K,1)
     %             + FB(4) * TG2AR(K,2)
 125     CONTINUE
         DO 127 K=1,3
            TGPURE(K,2,I,NBS2) = COSO(K,I,NBS2) - F(K,1)
     S    + XM *(TGPURE(K,2,NBS1,NBS2)-COSO(K,NBS1,NBS2)+COSO(K,NBS1,1))
     S    + XM1*(TGPURE(K,2,1,NBS2)-COSO(K,1,NBS2)+COSO(K,1,1))
 127     CONTINUE
 130  CONTINUE
C
C     LA RECHERCHE DES 2 TANGENTES DES SOMMETS DU COTE 4 DU CARRE UNITE
C     -----------------------------------------------------------------
      DO 140 I=2,NBA2
C
C        Ds/Dv AU POINT I DU COTE 4 DU CARRE
         CALL TG2ARL( I, COSO(1,1,I), COSO(1,1,I+1),
     S                MNXYT4, MNNTG4, TG2AR )
C        LA LONGUEUR DE L'ARETE I DU COTE 4 DU CARRE
         S1 = STCARR(2,1,I+1) - STCARR(2,1,I)
         TGPURE(1,2,1,I) = TG2AR(1,1) / S1
         TGPURE(2,2,1,I) = TG2AR(2,1) / S1
         TGPURE(3,2,1,I) = TG2AR(3,1) / S1
C
C        Ds/Du AU POINT I DU COTE 4 DU CARRE PAR LA FORMULE DE GORDON HALL
C        L'INTERVALLE DE YM SUR LE COTE 2 DU CARRE
         YM  = STCARR(2,1,I)
         YM1 = 1.0 - YM
         N2  = 1
C        INTERVALLE SUR LE COTE 2 CONTENANT YM
         DO 133 I1=2,NBS2
            IF( YM .LE. STCARR(2,NBS1,I1)+EPSXYZ ) THEN
               N2 = I1
               GO TO 134
            END IF
 133     CONTINUE
C        LES NUMEROS DES SOMMETS VOISINS
 134     N1 = N2-1
C        LES 2 TANGENTES DE L'ARETE N1 DU COTE 2
         CALL TG2ARL( N1, COSO(1,NBS1,N1), COSO(1,NBS1,N2),
     S                MNXYT2, MNNTG2, TG2AR )
C        LE PARAMETRE ENTRE 0 ET 1
         S1  = (   YM               - STCARR(2,NBS1,N1) )
     S        / ( STCARR(2,NBS1,N2) - STCARR(2,NBS1,N1) )
C        LA VALEUR DES POLYNOMES DE BASE P3 HERMITE AU POINT S1
         CALL VFBP3H( S1, FB )
         DO 135 K=1,3
C           LA COORDONNEE K DU POINT SUR LE COTE 2 DU CARRE
            F(K,2) = FB(1) * COSO(K,NBS1,N1)
     %             + FB(2) * COSO(K,NBS1,N2)
     %             + FB(3) * TG2AR(K,1)
     %             + FB(4) * TG2AR(K,2)
 135     CONTINUE
         DO 137 K=1,3
            TGPURE(K,1,1,I) = F(K,2) - COSO(K,1,I)
     S    + YM  * (TGPURE(K,1,1,NBS2)-COSO(K,NBS1,NBS2)+COSO(K,1,NBS2))
     S    + YM1 * (TGPURE(K,1,1,1)-COSO(K,NBS1,1)+COSO(K,1,1))
 137     CONTINUE
 140  CONTINUE
C
C     LE TRAITEMENT DES TANGENTES DE CHACUN DES QUADRANGLES
C     =====================================================
      DO 250 J=1,NBA2
         DO 240 I=1,NBA1
            NTG = 0
            DO 235 KJ=0,1
               DO 230 KI=0,1
C
C                 AU SOMMET (I+KI,J+KJ) DU QUADRANGLE I,J DU CARRE UNITE
C                 LES 2 COORDONNEES DU SOMMET I,J DU CARRE UNITE
                  IKI = I + KI
                  JKJ = J + KJ
                  XM = ST4C(1,KI,KJ)
                  XM1 = 1 - XM
                  YM = ST4C(2,KI,KJ)
                  YM1 = 1 - YM
C
C                 CALCUL DE LA MATRICE JACOBIENNE DE F:CARRE UNITE->QUADRANGLE D
                  DF(1,1)=YM1 * ( STCARR(1,I+1,J)  -STCARR(1,I,J)   )
     %                   +YM  * ( STCARR(1,I+1,J+1)-STCARR(1,I,J+1) )
                  DF(1,2)=YM1 * ( STCARR(2,I+1,J)  -STCARR(2,I,J)   )
     %                   +YM  * ( STCARR(2,I+1,J+1)-STCARR(2,I,J+1) )
                  DF(2,1)=XM1 * ( STCARR(1,I,J+1)  -STCARR(1,I,J)   )
     %                   +XM  * ( STCARR(1,I+1,J+1)-STCARR(1,I+1,J) )
                  DF(2,2)=XM1 * ( STCARR(2,I,J+1)  -STCARR(2,I,J)   )
     %                   +XM  * ( STCARR(2,I+1,J+1)-STCARR(2,I+1,J) )
C
C                | dS/dU |   |dF1/dU  dF2/dU| |ds/du|
C                |       | = |              | |     |
C                | dS/dV |   |dF1/dV  dF2/dV| |ds/dv|
C
                  NTG = NTG + 2
                  DO 220 K=1,3
                     S0 = TGPURE(K,1,IKI,JKJ)
                     S1 = TGPURE(K,2,IKI,JKJ)
                     COTG(K,NTG-1,I,J) = DF(1,1) * S0 + DF(1,2) * S1
                     COTG(K,NTG  ,I,J) = DF(2,1) * S0 + DF(2,2) * S1
 220              CONTINUE
 230           CONTINUE
 235        CONTINUE
C
C           RECTIFICATION DES SENS DES TANGENTES SUR LE QUADRANGLE I,J
            DO 238 K=1,3
                S0 = COTG(K,4,I,J)
C               LA TG CALCULEE 3 EST EN FAIT LA -TG4 DU QUADRANGLE I,J
                COTG(K,4,I,J) = -COTG(K,3,I,J)
C               LA TG CALCULEE 4 EST EN FAIT LA +TG3 DU QUADRANGLE I,J
                COTG(K,3,I,J) = S0
C
C               LA TG CALCULEE 5 EST EN FAIT LA  TG8 DU QUADRANGLE I,J
                S0 = COTG(K,8,I,J)
                COTG(K,8,I,J) = COTG(K,5,I,J)
C               LA TG CALCULEE 8 EST EN FAIT LA -TG6 DU QUADRANGLE I,J
                S1 = COTG(K,6,I,J)
                COTG(K,6,I,J) = -S0
C               LA TG CALCULEE 6 EST EN FAIT LA -TG7 DU QUADRANGLE I,J
                S0 = COTG(K,7,I,J)
                COTG(K,7,I,J) = -S1
C               LA TG CALCULEE 7 EST EN FAIT LA -TG5 DU QUADRANGLE I,J
                COTG(K,5,I,J) = -S0
 238        CONTINUE
C
 240     CONTINUE
 250  CONTINUE
      RETURN
      END
