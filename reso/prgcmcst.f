      SUBROUTINE PRGCMCST( NBSOM,  NONOSO, MNTOPO, MNNPEF, MNXYZN,
     %                     NBDLNO, NCODSA, MNLPLI, MNLPCO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE STOCKAGE CONDENSE DE LA MATRICE DU MAILLAGE POUR
C -----    NBDLNO DEGRES DE LIBERTE ( CONSTANT EN TOUT SOMMET )
C          A PARTIR D'UN MAILLAGE P2 POUR LE MAILLAGE DES SEULS SOMMETS
C          C'EST A DIRE POUR UNE INTERPOLATION P1 AUX SOMMETS
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NONOSO : NONOSO(I)=NO SOMMET SI SOMMET, 0 SI MILIEU D'ARETE
C MNTOPO : ADRESSE  MCN DU TABLEAU TOPOLOGIE   DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"XXXX DE L'OBJET
C MNXYZN : ADRESSE  MCN DU TABLEAU XYZNOEUD    DE L'OBJET
C NBDLNO : NOMBRE CONSTANT DE DEGRES DE LIBERTE PAR NOEUD
C NCODSA : CODE DE STOCKAGE DE LA MATRICE PROFIL
C           0 : MATRICE DIAGONALE
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C
C SORTIES:
C --------
C MNLPLI : ADRESSE MCN DU POINTEUR SUR LES COEFFICIENTS DIAGONAUX
C          DE LA MATRICE DU SYSTEME LINEAIRE POUR LES NBSOM SOMMETS
C MNLPCO : ADRESSE MCN DE LA LISTE DES NUMEROS DE COLONNES
C          DES COEFFICIENTS NON NULS DE LA MATRICE POUR LES NBSOM SOMMETS
C IERR   : 0 SI PAS D'ERREUR
C          1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     Avril 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/gsmenu.inc"
C
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
      INTEGER         MNNPEF(1:*), NONOSO(1:*)
C
      IERR = 0
      IF( NCODSA .EQ. 0 ) RETURN
C
C     FORMATION DE LA LISTE DES VOISINS
C     =================================
      MNLIVO = 0
      MNLPVO = 0
      CALL LISVOI( MNTOPO, MNNPEF, MNXYZN,
     %             MNLIVO, MNLPVO, MXVOIS, IERR  )
C     MXVOIS NOMBRE MAXIMAL DE VOISINS D'UN NOEUD DU MAILLAGE
      IF( IERR .NE. 0 ) GOTO 9000
C
C     LE NOMBRE TOTAL DE NOEUDS ET DE DEGRES DE LIBERTE DE L'OBJET
      NBNOE  = MCN( MNXYZN + WNBNOE )
      LOLIVO = MCN( MNLPVO + NBNOE )
C
C     FILTRAGE NONOSO DU NUMERO DES NOEUDS VOISINS
C     LIMITES AUX SEULS SOMMETS DES EF
C     ============================================
      NUS = 0
      DO I=1,NBNOE
C
C        NUMERO SOMMET DU NOEUD I
         NOSOI = NONOSO(I)
         IF( NOSOI .GT. 0 ) THEN
C
C           LE NOEUD I EST LE SOMMET NOSOI
C           DEBUT DES VOISINS DU NOEUD I
            MN1 = MNLIVO + MCN(MNLPVO+I-1)
C           FIN DES VOISINS DU NOEUD I
            MN2 = MNLIVO-1 + MCN(MNLPVO+I)
C
C           RECHERCHE DES VOISINS SOMMETS
            DO MNV=MN1,MN2
C              NUMERO DU NOEUD VOISIN
               NUMV = MCN( MNV )
C              NUMERO DU SOMMET VOISIN
               NUSV = NONOSO( NUMV )
               IF( NUSV .GT. 0 ) THEN
C                 AJOUT DU NO DE SOMMET VOISIN
                  MCN(MNLIVO+NUS) = NUSV
                  NUS = NUS + 1
               ENDIF
            ENDDO
C
C           MISE A JOUR DU POINTEUR SUR LE DERNIER
C           SOMMET VOISIN DU SOMMET NOSOI
            MCN(MNLPVO+NOSOI) = NUS
C
         ENDIF
C
      ENDDO
C
C     LE POINTEUR SUR LES LIGNES SOMMETS
C     ==================================
      NTDL = NBSOM * NBDLNO
      CALL TNMCDC( 'ENTIER', NTDL+1, MNLPLI )
      IF( MNLPLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PRGCMCST: PAS ASSEZ DE MEMOIRE MCN'
         ELSE
            KERR(1) = 'PRGCMCST: NOT ENOUGH MEMORY MCN'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF
C
C     MAJORATION DU NOMBRE DE COLONNES A STOCKER
C     EGAL AU NOMBRE DE COEFFICIENTS NON NUL DE LA MATRICE
      NUMV = MCN(MNLPVO+NBSOM)
      IF( NCODSA .GT. 0 ) THEN
C        MATRICE SYMETRIQUE
         LOLPCI = NBDLNO**2 * ( NUMV + NBSOM*2 + MXVOIS ) / 2
      ELSE
C        MATRICE NON SYMETRIQUE
         LOLPCI = NBDLNO**2 * ( NUMV + NBSOM )
      ENDIF
C
C     LE POINTEUR SUR LES COLONNES
      CALL TNMCDC( 'ENTIER', LOLPCI, MNLPCO )
      IF( MNLPCO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SP PRGCMCST: PAS ASSEZ DE MEMOIRE MCN'
         ELSE
            KERR(1) = 'SP PRGCMCST: NOT ENOUGH MEMORY MCN'
         ENDIF
         CALL LEREUR
         IERR = 1
         IF( MNLPLI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLI )
         GOTO 9000
      ENDIF
C
C     FORMATION DES TABLEAUX LPLIGN ET LPCOLO
C     =======================================
      MCN(MNLPLI) = 0
      IAN = 0
      IA  = 0
      IF( NCODSA .GT. 0 ) THEN
C
C        MATRICE SYMETRIQUE
C        ------------------
         MN1 = MNLIVO + MCN(MNLPVO)
         DO I=1,NBSOM
            MN2 = MNLIVO + MCN(MNLPVO+I) - 1
            IANVOI = 0
C           CALCUL DU NOMBRE DE VOISINS DU NOEUD I
            NBV = 0
            DO MNV=MN1,MN2
C              NUMERO DU NOEUD VOISIN
               NUMV = MCN(MNV)
               IF( NUMV .LT. I ) THEN
                  NBV = NBV + 1
               ENDIF
            ENDDO
C           PARCOURS DES NBV NOEUDS VOISINS DU NOEUD I
            DO MNV=MN1,MN1+NBV-1
C              NUMERO DU NOEUD VOISIN
               NUMV = MCN(MNV)
C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD NUMV
               NUMIV = ( NUMV - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,NBDLNO
                     IANV = IA + JN + ( ( IN - 1 ) * IN ) / 2
     &                    + NBV * NBDLNO * ( IN - 1 )
     &                    + IANVOI
                     MCN(MNLPCO+IANV-1) = NUMIV + JN
                  ENDDO
               ENDDO
               IANVOI = IANVOI + NBDLNO
            ENDDO
C           LES NBDLNO DEGRES DE LIBERTE DU NOEUD I
            NUMI = ( I - 1 ) * NBDLNO
            DO IN=1,NBDLNO
               DO JN=1,IN
                  IAN = IA + JN + ( ( IN - 1 ) * IN ) / 2
     &                + NBV * NBDLNO *  IN
                  MCN(MNLPCO+IAN-1) = NUMI + JN
               ENDDO
               MCN(MNLPLI+NUMI+IN) = IAN
            ENDDO
            IA  = IAN
            MN1 = MN2 + 1
         ENDDO
C
      ELSE
C
C        MATRICE NON SYMETRIQUE
C        ----------------------
         MN1 = MNLIVO + MCN(MNLPVO)
         DO I=1,NBSOM
            MN2 = MNLIVO + MCN(MNLPVO+I) - 1
            IANVOI = 0
            NBV = MN2 - MN1 + 1
            DO MNV=MN1,MN2
               NUMV = MCN(MNV)
C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD K
               NUMIV = ( NUMV - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,NBDLNO
                     IANV = IA + JN + IANVOI
     &               + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     MCN(MNLPCO+IANV-1) = NUMIV + JN
                  ENDDO
               ENDDO
               IANVOI = IANVOI + NBDLNO
            ENDDO
C           LES NBDLNO DEGRES DE LIBERTE DU NOEUD I
            IAN  = IA
            NUMI = ( I - 1 ) * NBDLNO
            DO IN=1,NBDLNO
               KN = 0
               DO JN=1,NBDLNO
                  IF (JN.NE.IN) THEN
                     KN = KN + 1
                     IAN = IA + KN + MCN(MNLPLI+I) * NBDLNO
     &                + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     MCN(MNLPCO+IAN-1) = NUMI + JN
                  ENDIF
               ENDDO
               IAN = IA + ( NBV + 1 ) * NBDLNO *  IN
               MCN(MNLPCO+IAN-1) = NUMI + IN
               MCN(MNLPLI+NUMI+IN) = IAN
            ENDDO
            IA  = IAN
            MN1 = MN2 + 1
         ENDDO
      ENDIF
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
C     ------------------------------------
      IF( MNLIVO .LE. 0 .OR. MNLPVO .LE. 0 ) GOTO 9000
      CALL TNMCDS( 'ENTIER', LOLIVO , MNLIVO )
      CALL TNMCDS( 'ENTIER', NBNOE+1, MNLPVO )
C
C     NOMBRE DE COEFFICIENTS DE LA MATRICE CONDENSEE
      LOLPCO = MCN(MNLPLI+NTDL)
      IF( LOLPCI .LT. LOLPCO ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PRGCMCST: VALEUR DE LOLPCI A AUGMENTER'
         ELSE
            KERR(1) = 'PRGCMCST: AUGMENT THE LOLPCI VALUE'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF
C
C     REDUCTION DE LA TAILLE DU TABLEAU DES NUMEROS DE COLONNES
      CALL TNMCRA( 'ENTIER', LOLPCI, LOLPCO, MNLPCO )
C
 9000 RETURN
      END
