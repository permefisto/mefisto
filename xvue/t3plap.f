      SUBROUTINE T3PLAP
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER DES 3 ARETES MANQUANTES ET VUES DU TRIEDRE
C -----  CELA TERMINE UN CALL T3PLAV
C        LES PARAMETRES DE TRACE SONT DANS $MEFISTO/incl/mecoit.inc
C        ET LE COMMON / T3PLAN /
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS         DECEMBRE 1997
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / T3PLAN / XYZTRI(3,2), H(3), NDIV(3), NOFA3P(3), NOET3P
      REAL              XYZ1(3), XYZ2(3)
C
C     PAS DE TRACE SI NOET3P N'EST PAS >0 CAR NON PASSAGE PAR T3PLAV
      IF( NOET3P .LE. 0 ) RETURN
C
C     PROTECTION SI NOET3P A ETE REGENERE ET LES TABLEAUX SONT INCORRECTS
      IF( NDIV(1) .LE. 0 ) RETURN
C
C     PAS DE TRACE SI NTY3PL EST NEGATIF OU NUL
      IF( NTY3PL .LE. 0 ) RETURN
C
C     L'EPAISSEUR DES TRAITS
      IF( NEP3PL .GE. 0 ) CALL XVEPAISSEUR( NEP3PL )
C
C     RECHERCHE DES 3 ARETES A TRACER
C     CE SONT LES 3 ARETES QUI N'APPRTIENNENT PAS AUX FACES TRACEES
C
C     LE NUMERO DE LA FACE ARRIERE OU DEVANT TRACEE
C     NOFA3P(1) = L   =>1 OU 2
C     L'ABSCISSE FIXE
C     X = XYZTRI(1,L)
C
C     LE NUMERO DE LA FACE GAUCHE OU DROITE TRACEE
C     NOFA3P(2) = L   =>1 OU 2
C     L'ORDONNEE FIXE
C     Y = XYZTRI(2,L)
C
C     LE NUMERO DE LA FACE BASSE OU HAUTE TRACEE
C     NOFA3P(3) = L   =>1 OU 2
C     LA COTE FIXE
C     Z = XYZTRI(3,L)
C
      IF( NOFA3P(1) .EQ. 1 ) THEN
C
         IF( NOFA3P(2) .EQ. 1 ) THEN
C
            IF( NOFA3P(3) .EQ. 1 ) THEN
C
C              SOMMET COMMUN 7 A JOINDRE AUX SOMMETS 3 6 8 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,2)
               XYZ1(2) = XYZTRI(2,2)
               XYZ1(3) = XYZTRI(3,2)
C
C              LE SOMMET 3
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 6
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 8
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ELSE
C
C              SOMMET COMMUN 3 A JOINDRE AUX SOMMETS 2 4 7 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,2)
               XYZ1(2) = XYZTRI(2,2)
               XYZ1(3) = XYZTRI(3,1)
C
C              LE SOMMET 2
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 4
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 7
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ENDIF
C
         ELSE
C
            IF( NOFA3P(3) .EQ. 1 ) THEN
C
C              SOMMET COMMUN 6 A JOINDRE AUX SOMMETS 2 5 7 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,2)
               XYZ1(2) = XYZTRI(2,1)
               XYZ1(3) = XYZTRI(3,2)
C
C              LE SOMMET 2
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 5
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 7
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ELSE
C
C              SOMMET COMMUN 2 A JOINDRE AUX SOMMETS 1 3 6 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,2)
               XYZ1(2) = XYZTRI(2,1)
               XYZ1(3) = XYZTRI(3,1)
C
C              LE SOMMET 1
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 3
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 6
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ENDIF
C
         ENDIF
C
      ELSE
C
         IF( NOFA3P(2) .EQ. 1 ) THEN
C
            IF( NOFA3P(3) .EQ. 1 ) THEN
C
C              SOMMET COMMUN 8 A JOINDRE AUX SOMMETS 4 5 7 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,1)
               XYZ1(2) = XYZTRI(2,2)
               XYZ1(3) = XYZTRI(3,2)
C
C              LE SOMMET 4
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 5
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 7
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ELSE
C
C              SOMMET COMMUN 4 A JOINDRE AUX SOMMETS 1 3 8 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,1)
               XYZ1(2) = XYZTRI(2,2)
               XYZ1(3) = XYZTRI(3,1)
C
C              LE SOMMET 1
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 3
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 8
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ENDIF
C
         ELSE
C
            IF( NOFA3P(3) .EQ. 1 ) THEN
C
C              SOMMET COMMUN 5 A JOINDRE AUX SOMMETS 1 6 8 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,1)
               XYZ1(2) = XYZTRI(2,1)
               XYZ1(3) = XYZTRI(3,2)
C
C              LE SOMMET 1
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 6
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 8
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ELSE
C
C              SOMMET COMMUN 1 A JOINDRE AUX SOMMETS 2 4 5 DE L'HEXAEDRE
C              =========================================================
               XYZ1(1) = XYZTRI(1,1)
               XYZ1(2) = XYZTRI(2,1)
               XYZ1(3) = XYZTRI(3,1)
C
C              LE SOMMET 2
               XYZ2(1) = XYZTRI(1,2)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 4
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,2)
               XYZ2(3) = XYZTRI(3,1)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
C              LE SOMMET 5
               XYZ2(1) = XYZTRI(1,1)
               XYZ2(2) = XYZTRI(2,1)
               XYZ2(3) = XYZTRI(3,2)
C              LE TRACE DU TRAIT DE L'ARETE VUE
               CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
            ENDIF
C
         ENDIF
C
      ENDIF
C
C     ETAT DE SORTIE CORRECT LES TABLEAUX DU COMMON/T3PLAN/ ONT ETE TRAITES
      NOET3P = 0
      RETURN
      END
