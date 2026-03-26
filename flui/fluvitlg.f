      SUBROUTINE FLUVITLG( MNXYZP, MNELE,  NODDLN, VXYZPN,
     %                     NUMILG, NUMALG,
     %                     FLUNEG, FLUPOS, FLUVLG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE FLUX DE LA VITESSE A TRAVERS LES ARETES
C -----    FRONTALIERES DECLAREES DANS LA DEFINITION DE L'OBJET

C ENTREES:
C --------
C MNXYZP : ADRESSE MCN DU TABLEAU XYZPOINT DE L'OBJET KNOMOB
C MNELE  : ADRESSE DU TABLEAU NPEF"TYPE EF DU MAILLAGE 2D
C NODDLN : NUMERO DU DERNIER DL DE CHAQUE NOEUD
C VXYZPN : VITESSE PRESSION EN TOUS LES NOEUDS
C NUMILG : NUMERO MINIMAL DES LIGNES DANS LA DEFINITION DE L'OBJET
C NUMALG : NUMERO MAXIMAL DES LIGNES DANS LA DEFINITION DE L'OBJET

C SORTIES:
C --------
C FLUNEG : FLUX NEGATIF DE LA VITESSE A TRAVERS LES ARETES FRONTALIERES
C FLUPOS : FLUX POSITIF DE LA VITESSE A TRAVERS LES ARETES FRONTALIERES
C FLUVLG : FLUX DE LA VITESSE A TRAVERS LES LIGNES DE L'OBJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Paris &St PIERRE du PERRAY Mai 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/ponoel.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           NODDLN(0:*)
      INTEGER           NONOEF(10), NONOAK(3)
      EQUIVALENCE      (NONOAK(1),N1), (NONOAK(2),N2), (NONOAK(3),N3)
      INTEGER           NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
      DOUBLE PRECISION  VXYZPN(*)
      DOUBLE PRECISION  DN(3), FLUPOS, FLUNEG, FLUXV, VN
      DOUBLE PRECISION  FLUVLG(NUMILG:NUMALG)
      REAL              COORDP(18)
      CHARACTER*24      KNOMLG

C     INITIALISATION DES FLUX A ZERO
      FLUNEG = 0D0
      FLUPOS = 0D0
      DO N = NUMILG, NUMALG
         FLUVLG( N ) = 0D0
      ENDDO

C     NOMBRE D'ELEMENTS FINIS DU MAILLAGE
      NBELEM = MCN( MNELE + WBELEM )

C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )

C     MNPGEL ADRESSE MCN DES NUMEROS NOEUDS ET POINTS GEOMETRIQUES DES EF
      MNPGEL = MNELE + WUNDEL

C     LA BOUCLE SUR LES ELEMENTS FINIS POUR CALCULER LE FLUX A LA FRONTIERE
C     =====================================================================
      DO 100 NUELEM = 1, NBELEM

C        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C
C        NO DES POINTS LIGNES LIGNES VOLUMES DES SOMMETS ARETES FACES VOLUME
         CALL EFPLSV( MNELE , NUELEM,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )

C        UNE ARETE DU TRIANGLE EST ELLE FRONTALIERE?
         DO KL = 1, 3
            IF( NOOBLA(KL) .GT. 0 ) GOTO 10
         ENDDO
         GOTO 100

C        INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C        COORDONNEES DES POINTS=NOEUDS DE L'ELEMENT FINI
C        INTEGRATION NUMERIQUE AUX 2 SOMMETS DES ARETES DES LIGNES
 10      CALL EFXYZP( 2, MNXYZP, NBELEM, NUELEM, MNPGEL, NBNOEF,
     %                COORDP )
         NBNOE2 = 2 * NBNOEF
C
C        CALCUL DU VECTEUR NORMAL A L'ARETE K
         DO K = KL, 3

            IF( NOOBLA(K) .GT. 0 ) THEN

C              RECHERCHE DU JACOBIEN DE G
C              NUMERO LOCAL DES NOEUDS DE L'ARETE
               N1 = K
               IF( K .NE. 3 ) THEN
                  N2 = K+1
               ELSE
                  N2 = 1
               ENDIF
               N3 = N1 + 3

C              CALCUL DE LA NORMALE * DELTAK: DG/DY, -DG/DX
               DN(1) = COORDP( N2+NBNOEF )  - COORDP( N1+NBNOEF )
               DN(2) = COORDP( N1 )         - COORDP( N2 )
C
C              FLUX NORMAL DE LA VITESSE A TRAVERS L'ARETE K
               FLUXV = 0D0
               IF( NUTYEL .EQ. 13 ) THEN
C                 TRIANGLE BREZZI-FORTIN INTEGRATION AUX 2 SOMMETS
                  DO I = 1, 2
C                    NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                     NODL = NODDLN( NONOEF( NONOAK(I) ) - 1 )
C                    Vitesse Normale  DELTAK  au NOEUD I de la ARETE K
                     VN = VXYZPN(NODL+1) * DN(1)
     %                  + VXYZPN(NODL+2) * DN(2)
                     FLUXV = FLUXV + VN
                  ENDDO
                  FLUXV = FLUXV / 2D0
               ELSE
C                 TRIANGLE TAYLOR HOOD INTEGRATION AUX 3 NOEUDS DE L'ARETE
C                 POIDS 1/6 AUX 2 SOMMETS et 2/3 AU MILIEU DE L'ARETE
                  DO I = 1, 2
C                    NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                     NODL = NODDLN( NONOEF( NONOAK(I) ) - 1 )
C                    Vitesse Normale  DELTAK  au NOEUD I de la ARETE K
                     VN = VXYZPN(NODL+1) * DN(1)
     %                  + VXYZPN(NODL+2) * DN(2)
                     FLUXV = FLUXV + VN/6D0
                  ENDDO
C                 LE MILIEU DE L'ARETE
                  NODL = NODDLN( NONOEF( NONOAK(3) ) - 1 )
                  VN = VXYZPN(NODL+1) * DN(1)
     %               + VXYZPN(NODL+2) * DN(2)
                  FLUXV = FLUXV + VN*2D0/3D0
               ENDIF

C              REPARTITION EN FLUX NEGATIF ET POSITIF ET POIDS
               IF( FLUXV .LE. 0D0 ) THEN
                  FLUNEG = FLUNEG + FLUXV
               ELSE
                  FLUPOS = FLUPOS + FLUXV
               ENDIF

C              FLUX DE LA VITESSE SUR LA LIGNE NOOBLA(K)
               FLUVLG( NOOBLA(K) ) = FLUVLG( NOOBLA(K) ) + FLUXV

            ENDIF

         ENDDO
 100  CONTINUE

C     AFFICHAGE DU FLUX DE LA VITESSE SUR LES LIGNES DE L'OBJET
C     =========================================================
      DO K = NUMILG, NUMALG
C        LE FLUX DE LA VITESSE SUR LA LIGNE K A-T-ELLE ETE CALCULEE
         FLUXV = FLUVLG( K )
         IF( FLUXV .NE. 0D0 ) THEN
C            OUI: IMPRESSION DE L'INTEGRALE
C            LE NOM DE LA LIGNE
             CALL NMOBNU( 'LIGNE', K, KNOMLG )
             NN = NUDCNB( KNOMLG )
             IF( LANGAG .EQ. 0 ) THEN
                WRITE(IMPRIM,*)'FLUX de la VITESSE=',FLUXV,
     %                ' sur la LIGNE ',KNOMLG(1:NN)
             ELSE
                WRITE(IMPRIM,*)'VELOCITY FLUX=',FLUXV,
     %                ' through the LIGNE ',KNOMLG(1:NN)
             ENDIF
         ENDIF
      ENDDO

C     AFFICHAGE DU TOTAL DES FLUX POSITIFS et NEGATIFS DE LA VITESSE
C     SUR LES LIGNES DE L'OBJET
      FLUXV = FLUNEG + FLUPOS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) FLUNEG, FLUPOS, FLUXV
      ELSE
         WRITE(IMPRIM,20100) FLUNEG, FLUPOS, FLUXV
      ENDIF

10100 FORMAT(' FLUX<0 de la VITESSE=',
     %G14.6,'  FLUX>0=',G14.6,'  FLUX TOTAL=',G14.6)
20100 FORMAT(' VELOCITY NEGATIVE FLUX=',G14.6,
     %'  POSITIVE FLUX=',G14.6,'  TOTAL FLUX=',G14.6)

      RETURN
      END
