      SUBROUTINE TS6Q1C( NBJEUX, JEU,    NBPOLY, NPI,    POLY,
     %                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                   F,      POIDEL, DP,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL  DU SECOND MEMBRE D'UN EF 6D
C -----    6CUBE  6Q1C  LAGRANGE DE DEGRE 1 ISOPARAMETRIQUE
C
C ENTREES:
C --------
C NOEF   : NUMERO DU TYPE DE L'EF
C X      : COORDONNEES X Y Z DES NBPOLY POINTS DE L'EF
C
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE SUR LE RECTANGLE
C POLY   : VALEUR DES POLYNOMES AUX POINTS D'INTEGRATION DE L'EF
C
C NOOBVO : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMIVO : NUMERO MINIMAL DES OBJETS SURFACES
C NUMAVO : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C F      : XX YY ZZ DES NPI POINTS D'INTEGRATION
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATI
C DP     : DP(3, NBPOLY, NPI) GRADIENT AUX POINTS D 'INTEGRATION DES
C          FONCTIONS DE BASE LAGRANGE ISOPARAMETRIQUES
C X      : XX YY ZZ DES NBPOLY POINTS DE L'EF
C
C SORTIE :
C --------
C BE     : BE(NBPOLY) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  SOURCE,
     %                  POLY(NBPOLY, NPI),
     %                  POIDEL(NPI),
     %                  DP(6, NBPOLY, NPI)
      DOUBLE PRECISION  F( NPI, 6 ),
     %                  BE(NBPOLY),
     %                  VITEFL(6),
     %                  PROSCD,
     %                  XYZUVW(6),
     %                  D
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
C
C     MISE A ZERO DE BE LE VECTEUR ELEMENTAIRE
C     ----------------------------------------
      CALL AZEROD( NBPOLY, BE )
C
C     ===================================
C     CONTRIBUTION DES SOURCES VOLUMIQUES
C     ===================================
      IF( LTDEVO(LPSOUR,JEU,NOOBVO) .GT. 0 ) THEN
C
         IF( TESTNL .GE. 1 ) THEN
C           LES SOURCES DEPENDENT DE LA TEMPERATURE
C           RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
            MN  = (MNTHET-1)/2
            MNT = (MNTHDL-1)/2
            DO 2 I=1,NBPOLY
               DMCN( MNT+I )=DMCN( MN+MCN(MNNODL+(I-1)) )
 2          CONTINUE
         ENDIF
C
         DO 9 L=1,NPI
C
C           LA CONTRIBUTION DES SOURCES DE CHALEUR
            IF( TESTNL .GE. 1 ) THEN
C              CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL=PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
            ENDIF
C
C           LA VALEUR DES SOURCES VOLUMIQUES EN CE POINT L
            XYZUVW(1) = F(L,1)
            XYZUVW(2) = F(L,2)
            XYZUVW(3) = F(L,3)
            XYZUVW(4) = F(L,4)
            XYZUVW(5) = F(L,5)
            XYZUVW(6) = F(L,6)
            CALL RESOUR( 4, NOOBVO, 6, XYZUVW,
     %                   LTDEVO(LPSOUR,JEU,NOOBVO), SOURCE )
C
            DO 7 I=1,NBPOLY
               BE(I) = BE(I) + POIDEL(L) * POLY(I,L) * SOURCE
  7         CONTINUE
C
  9      CONTINUE
C
      ENDIF
C
C     LA CONTRIBUTION DU TRANSPORT: - VITESSE * GRADIENT TEMPERATURE
C     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
C     -----------------------------------------------------------------
      IF( LTDEVO(LPVIFL,JEU,NOOBVO) .GT. 0 ) THEN
C
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN  = (MNTHET-1)/2
         MNT = (MNTHDL-1)/2
         DO 12 I=1,NBPOLY
            DMCN( MNT+I ) = DMCN( MN+MCN(MNNODL+I-1) )
 12      CONTINUE
C
         DO 18 L=1,NPI
C
C           LA VALEUR DE LA VITESSE DU FLUIDE AU POINT D'INTEGRATION L
            DO 16 I=1,6
               XYZUVW(I) = F(L,I)
 16         CONTINUE
            CALL REVIFL( 4,NOOBVO, 6, 6, XYZUVW,
     %                   LTDEVO(LPVIFL,JEU,NOOBVO), VITEFL )
C
C           - t[P] [V] [DP] [TEMPERATURE EF]
            DO 17 I=1,NBPOLY
               D = 0D0
               DO 15 K=1,NBPOLY
                  D = D + ( VITEFL(1)*DP(1,K,L) +
     %                      VITEFL(2)*DP(2,K,L) +
     %                      VITEFL(3)*DP(3,K,L) +
     %                      VITEFL(4)*DP(4,K,L) +
     %                      VITEFL(5)*DP(5,K,L) +
     %                      VITEFL(6)*DP(6,K,L) ) * DMCN(MNT+K)
 15            CONTINUE
                BE(I) = BE(I) - POIDEL(L) * POLY(I,L) * D
 17         CONTINUE
 18      CONTINUE
C
      ENDIF
C
      RETURN
      END
