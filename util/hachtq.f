      SUBROUTINE HACHTQ( NVALEU, VALEUR, NU1TG, L1, L2, LISTE, LIEN,
     %                   NOCOLO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHER LES NVALEU ENTIERS DU TABLEAU VALEUR PARMI
C ----- LES NVALEU-ERES VALEURS DES COLONNES DU TABLEAU LISTE
C       LA METHODE EMPLOYEE ICI EST CELLE DU HACHAGE
C       ADRESSAGE SUR LA SOMME DES VALEURS MODULO L2
C
C       CAS PARTICULIER DU TRIANGLE OU QUADRANGLE OU IL FAUT CHERCHER
C       DANS UN SENS ET DANS L'AUTRE LES NUMEROS DES SOMMETS
C
C ENTREES:
C --------
C NVALEU : NOMBRE DES VALEURS A IDENTIFIER OU AJOUTER
C          = 3 OU 4  SINON => SORTIE AVEC NOCOLO=0
C VALEUR : TABLEAU DE NVALEU ENTIERS DES VALEURS A RETROUVER OU AJOUTER
C NU1TG  : NUMERO DANS VALEUR DE LA PREMIERE TANGENTE
C          SI =0 PAS DE TANGENTES, SINON 2 TGS PAR ARETE
C L1     : NOMBRE DE LIGNES DU TABLEAU LISTE (>NVALEU)
C L2     : NOMBRE DE COLONNES DU TABLEAU LISTE
C LISTE  : TABLEAU CONTENANT LES VALEURS DEJA AJOUTEES
C LIEN   : NUMERO DE LA LIGNE CONTENANT LE CHAINAGE SUR LE SUIVANT
C          ET SANS SUIVANT OU LIGNE NON UTILISEE  : LA VALEUR 0
C
C SORTIES:
C --------
C VALEUR : TABLEAU DE NVALEU ENTIERS DES VALEURS EVENTUELLEMENT PERMUTE
C NOCOLO : NO DE LA COLONNE DU TABLEAU VALEUR RETROUVEE
C          > 0 SI LE TABLEAU VALEUR A ETE RETROUVE
C          =<0 SI LE TABLEAU VALEUR N'A PAS ETE RETROUVE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1994
C ......................................................................
      INTEGER  VALEUR(1:*) , LISTE(L1,L2)
C
C     LA FONCTION D ADRESSAGE
C     =======================
      IF( NVALEU .EQ. 3 ) THEN
C
C         LE TRIANGLE :
C         =============
C         TRI CROISSANT DES 3 NUMEROS
          CALL HACREN( NVALEU, NVALEU, NU1TG, VALEUR )
C
          NOCOLO = VALEUR(1) + VALEUR(2) + VALEUR(3)
          NOCOLO = MOD( ABS(NOCOLO) , L2 )
          IF( NOCOLO .EQ. 0 ) NOCOLO = L2
C
C         LA RECHERCHE DU TABLEAU
 10       IF( VALEUR(1) .EQ. LISTE(1,NOCOLO)  ) THEN
             IF( VALEUR(2) .EQ. LISTE(2,NOCOLO) ) THEN
                IF( VALEUR(3) .EQ. LISTE(3,NOCOLO) ) THEN
C                  TRIANGLE RETROUVE SENS DIRECT
                   RETURN
                ENDIF
             ELSE IF( VALEUR(2) .EQ. LISTE(3,NOCOLO) ) THEN
                IF( VALEUR(3) .EQ. LISTE(2,NOCOLO) ) THEN
C                  TRIANGLE RETROUVE SENS RETROGRADE
                   IF( NU1TG .GT. 0 ) THEN
                      NU0TG = NU1TG - 1
C                     PERMUTATION DES TANGENTES
                      L               = VALEUR(NU0TG+1)
                      VALEUR(NU0TG+1) = VALEUR(NU0TG+2)
                      VALEUR(NU0TG+2) = L
                      L               = VALEUR(NU0TG+3)
                      VALEUR(NU0TG+3) = VALEUR(NU0TG+6)
                      VALEUR(NU0TG+6) = L
                      L               = VALEUR(NU0TG+4)
                      VALEUR(NU0TG+4) = VALEUR(NU0TG+5)
                      VALEUR(NU0TG+5) = L
                   ENDIF
                   RETURN
                ENDIF
             ENDIF
          ENDIF
C
C         TABLEAU VALEUR RECHERCHE PARMI LES COLONNES CHAINEES A NOCOLO
          NOCOLO = LISTE( LIEN , NOCOLO )
          IF( NOCOLO .GT. 0 ) GOTO 10
C         VALEUR NON RETROUVEE
          RETURN
C
      ELSE IF( NVALEU .EQ. 4 ) THEN
C
C         LE QUADRANGLE :
C         ===============
C         TRI CROISSANT DES 4 NUMEROS
          CALL HACREN( NVALEU, NVALEU, NU1TG, VALEUR )
C
          NOCOLO = VALEUR(1) + VALEUR(2) + VALEUR(3) + VALEUR(4)
          NOCOLO = MOD( ABS(NOCOLO) , L2 )
          IF( NOCOLO .EQ. 0 ) NOCOLO = L2
C
C         LA RECHERCHE DU TABLEAU
 20       IF( VALEUR(1) .EQ. LISTE(1,NOCOLO)  ) THEN
             IF( VALEUR(2) .EQ. LISTE(2,NOCOLO) ) THEN
                IF( VALEUR(3) .EQ. LISTE(3,NOCOLO) ) THEN
                   IF( VALEUR(4) .EQ. LISTE(4,NOCOLO) ) THEN
C                     QUADRANGLE RETROUVE SENS DIRECT
                      RETURN
                   ENDIF
                ENDIF
             ELSE IF( VALEUR(2) .EQ. LISTE(4,NOCOLO) ) THEN
                IF( VALEUR(3) .EQ. LISTE(3,NOCOLO) ) THEN
                   IF( VALEUR(4) .EQ. LISTE(2,NOCOLO) ) THEN
C                     QUADRANGLE RETROUVE SENS RETROGRADE
                      IF( NU1TG .GT. 0 ) THEN
                         NU0TG = NU1TG - 1
C                        PERMUTATION DES TANGENTES
                         L               = VALEUR(NU0TG+1)
                         VALEUR(NU0TG+1) = VALEUR(NU0TG+2)
                         VALEUR(NU0TG+2) = L
                         L               = VALEUR(NU0TG+3)
                         VALEUR(NU0TG+3) = VALEUR(NU0TG+8)
                         VALEUR(NU0TG+8) = L
                         L               = VALEUR(NU0TG+4)
                         VALEUR(NU0TG+4) = VALEUR(NU0TG+7)
                         VALEUR(NU0TG+7) = L
                         L               = VALEUR(NU0TG+5)
                         VALEUR(NU0TG+5) = VALEUR(NU0TG+6)
                         VALEUR(NU0TG+6) = L
                      ENDIF
                      RETURN
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
C
C         TABLEAU VALEUR RECHERCHE PARMI LES COLONNES CHAINEES A NOCOLO
          NOCOLO = LISTE( LIEN , NOCOLO )
          IF( NOCOLO .GT. 0 ) GOTO 20
C         VALEUR NON RETROUVEE
          RETURN
C
      ELSE
C
C        CAS NON TRAITE
C        ==============
         NOCOLO = 0
      ENDIF

      RETURN
      END
