      SUBROUTINE NLSESCH( TPPDIAG, Unp1,   Un,     DTALFA, Ulapl,
     %                    DTBETA,  Vn,     Wn,     Rn,     DT,
     %                    BG,
     %                    NBNOFX,  NONOFX, VANOFX, NOPARTUnp1,
     %                    NBSOM,   XYZSOM, NBTRIA, NOSOTR,
     %                    NOOBSF,  NUMISU, NUMASU, LTDESU, NOPARTF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER UG SOLUTION DE
C ----- tPP Unp1 = tPP Un + dt alfa tDPDP Ulapl + tPU2P (Vn**2+Wn**2) Rn
C                + dt FOmega(t,X,Vn,Wn)
C       SUR UN MAILLAGE DE TRIANGLES 2P1D
C
C ENTREES:
C --------
C TPPDIAG: MATRICE DIAGONALE tP1P1 DIAGONALE
C Un, Ulapl, Vn, Wn, Rn: VECTEURS DECRITS DANS L'EQUATION CI-DESSUS
C BG     : SECOND MEMBRE GLOBAL NECESSAIRE A L'ASSEMBLAGE
C NBNOFX : NOMBRE DE NOEUDS FIXES (CONDITION DIRICHLET HOMOGENE)
C NONOFX : NBNOFX NUMEROS DES NOEUDS FIXES
C VANOFX : VALEUR FIXEE AU NOEUDS FIXES
C NOPARTUnp1: 1 PARTIE REELLE     de Unp1=V TRAITEE EN CONDITION DIRICHLET
C             2 PARTIE IMAGINAIRE de Unp1=W TRAITEE EN CONDITION DIRICHLET
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NBSOM TRIANGLES
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE
C NOSOTR : NO DES 3 SOMMETS DES NBTRIA TRIANGLES
C NOOBSF : NO DE LA SURFACE
C NUMISU : NO MINIMAL DES SURFACES DE L'OBJET
C NUMISU : NO MAXIMAL DES SURFACES DE L'OBJET
C LTDESU : ADRESSE DES TABLEAUX DES DONNEES SURFACIQUES DE L'ONDE NLSE
C NOPARTF: 1 PARTIE REELLE     de FOmega(1) a UTILISER
C          2 PARTIE IMAGINAIRE de FOmega(2) a UTILISER
C
C SORTIES:
C --------
C Unp1   : VECTEUR GLOBAL A CALCULER SOLUTION DE L'EQUATION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NOSOTR(NBTRIA,3), NONOFX(NBNOFX)
      INTEGER           LTDESU(1:MXDOTH, NUMISU:NUMASU)
      DOUBLE PRECISION  TPPDIAG(NBSOM), Unp1(NBSOM), Un(NBSOM), DT,
     %                  DTALFA,Ulapl(NBSOM), DTBETA,Vn(NBSOM),Wn(NBSOM),
     %                  Rn(NBSOM), BG(NBSOM), VANOFX(NBNOFX)
      DOUBLE PRECISION  tPP(3), TDPDP(3,3), tPU2P(3), FOmega(3,2)
C
C     MISE A ZERO DU SECOND MEMBRE
C     ============================
      DO NS = 1, NBSOM
         BG( NS ) = 0D0
      ENDDO
C
C     LA BOUCLE SUR LES TRIANGLES 2P1D
C     ================================
      DO NUELEM = 1, NBTRIA
C
C        LES INTEGRALES DES POLYNOMES DE BASE DU TRIANGLE 2P1D
         CALL NLSE2P1D( DT,     DTALFA, DTBETA, NBSOM, XYZSOM,
     %                  NBTRIA, NUELEM, NOSOTR,
     %                  NOOBSF, NUMISU, NUMASU, LTDESU,
     %                  Vn,     Wn,
     %                  tPP,    tDPDP,  tPU2P, FOmega )
C
C        tPP Unp1 = tPP Un + dt Alfa tDPDP Ulapl
C                 + dt Beta tPU2P (Vn**2+Wn**2) Rn + dt FOmega
         DO K = 1, 3
C
C           NUMERO DU SOMMET K DE L'EF
            NS = NOSOTR(NUELEM,K)
C
C           ASSEMBLAGE DU SECOND MEMBRE DANS LE SECOND MEMBRE GLOBAL
C           LE SIGNE - DEVANT tDPDP EST DU AU TERME +dt Alfa Laplacien
C           c-a-d UN SIGNE + DEVANT LE LAPLACIEN
            BG( NS ) = BG( NS )
     %               + tPP(K)     * Un(NS)
     %               - tDPDP(K,1) * Ulapl( NOSOTR(NUELEM,1) )
     %               - tDPDP(K,2) * Ulapl( NOSOTR(NUELEM,2) )
     %               - tDPDP(K,3) * Ulapl( NOSOTR(NUELEM,3) )
     %               + tPU2P(K)   * Rn(NS)
     %               + FOmega(K,NOPARTF)
         ENDDO
C
C        FIN DE BOUCLE SUR LES EF DE TYPE NUTYEL
      ENDDO
C
C     CONDITIONS AUX LIMITES DE DIRICHLET
C     ATTENTION: MEME VALEUR FIXEE POUR LA PARTIE REELLE ET IMAGINAIRE
C                A MODIFIER POUR AVOIR 2 VALEURS FIXEES DIFFERENTES
C                CF recont.f
      DO K=1,NBNOFX
C        NUMERO DU NOEUD K FIXE
         NS = NONOFX(K)
C        CONDITION DE DIRICHLET
         IF( NOPARTUnp1 .GE. 0 ) THEN
            BG(NS) = VANOFX(K)
         ENDIF
      ENDDO
C
C     TPPDIAG Unp1 = tPP Un + dt alfa tDPDP Ulapl
C                  + tPU2P (Vn**2+Wn**2) Rn + dt FOmega = BG
      DO NS = 1, NBSOM
         Unp1( NS ) = BG( NS ) / TPPDIAG( NS )
      ENDDO
C
      RETURN
      END
