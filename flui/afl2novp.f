      SUBROUTINE AFL2NOVP( NDIM,   NBVECT,   TIMES,
     %                     VOLUME, INTL2VIT, INTL2PRE,
     %                     NOFOVI, INTL2VER, NOFOPR, INTL2PER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    AFFICHAGE DE LA NORME L2 DE LA PRESSION, DE LA |VITESSE|
C ----    ET SI FONCTION VITESSE_EXACTE(t,x,y,z,nc)
C               LA NORME L2 DES ERREURS SUR LE MODULE DE LA VITESSE
C         ET SI FONCTION PRESSION_EXACTE(t,x,y,z)
C               LA NORME L2 DES ERREURS SUR LA PRESSION
C         POUR TOUS LES PAS DE TEMPS CALCULES
C
C ENTREES:
C --------
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
C NBVECT : NOMBRE DE VECTEURS DU TABLEAU VECTEUR"L2VITESSEPRESS
C VOLUME : SURFACE EN 2D, VOLUME EN 3D du DOMAINE MAILLE
C VOLUME : VOLUME DU MAILLAGE
C          Volume= Som  Delta(e)
C                  e dans E
C                                                          (Vk1C)
C INTL2VIT:Som Jacobien Som (Vk1C,...,VknC) [int Pi Pj DX] ( ...)
C          e de E       k=1,...,d                          (VknC)
C                                                    (P1C)
C INTL2PRE:Som Jacobien (P1C,...,PmC) [int Pi Pj DX] (...)
C          e de E                                    (PmC)
C
C NOFOVI : NUMERO DE LA FONCTION VITESSE_EXACTE(t,x,y,z,nc)
C          0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C                                                                    (Vk1E-Vk1C)
C INTL2VER:Som Jacobien Som (Vk1E-Vk1C,...,VknE-VknC) [int Pi Pj DX] (    ...  )
C          e de E       k=1,...,d                                    (VknE-VknC)
C          NON CALCULE SI PRESSION_EXACTE(t,x,y,z) NON DONNEE PAR L'UTILISATEUR
C
C NOFOPR : NUMERO DE LA FONCTION PRESSION_EXACTE(t,x,y,z)
C          0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C                                                            (P1E-P1C)
C INTL2PER:Som Jacobien (P1E-P1C,...,PmE-PmC) [int Pi Pj DX] (  ...  )
C          e de E                                            (PmE-PmC)
C          NON CALCULE SI VITESSE_EXACTE(t,x,y,z,nc) NON DONNEE PAR L'USAGER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  VOLUME, SQRVOL,
     %                  INTL2PRE(NBVECT), INTL2VIT(NBVECT),
     %                  INTL2PER(NBVECT), INTL2VER(NBVECT)
      REAL              TIMES(NBVECT)
      INTRINSIC         SQRT
C
C     AFFICHAGE DU VOLUME DU DOMAINE
C     ==============================
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'VOLUME du DOMAINE=',VOLUME
      ELSE
         WRITE(IMPRIM,*) 'VOLUME of the DOMAIN=',VOLUME
      ENDIF
      SQRVOL = SQRT( VOLUME )
C
C     AFFICHAGE DE LA NORME L2 DE LA PRESSION, DE LA |VITESSE|
C     ET SI FONCTION EXACTE DE L'UTILISATEUR LA NORME L2 DES ERREURS
C     ==============================================================
ccc      DO K=1,NBVECT
cccC
cccC        RACINE CARREE POUR OBTENIR LA NORME L2 ET
cccC        DIVISION PAR VOLUME POUR OBTENIR UNE VALEUR MOYENNE
ccc         INTL2VIT(K) = SQRT( INTL2VIT(K) / VOLUME )
ccc         INTL2PRE(K) = SQRT( INTL2PRE(K) / VOLUME )
ccc         INTL2VER(K) = SQRT( INTL2VER(K) / VOLUME )
ccc         INTL2PER(K) = SQRT( INTL2PER(K) / VOLUME )
cccC
ccc      ENDDO
C
10500 FORMAT(/148('=')/'NORMES ||VITESSE||L2 sur le domaine de SURFACE='
     %,G14.6 / 148('=') )
10501 FORMAT(/148('=')/'NORMES ||VITESSE||L2 sur le domaine de VOLUME='
     %,G14.6 / 148('=') )
C
20500 FORMAT(/148('=')/'||VELOCITY||L2 NORMS on the domain of SURFACE=',
     %G14.6 / 148('=') )
20501 FORMAT(/148('=')/'||VELOCITY||L2 NORMS on the domain of VOLUME=',
     %G14.6 / 148('=') )
C
10502 FORMAT(/148('=')/'NORMES||PRESSION||L2 sur le domaine de SURFACE='
     %,G14.6 /148('=') )
10503 FORMAT(/148('=')/'NORMES ||PRESSION||L2 sur le domaine de VOLUME='
     %,G14.6 /148('=') )
C
20502 FORMAT(/148('=')/'||PRESSURE||L2 NORMS on the domain of SURFACE ='
     %,G14.6 /148('=') )
20503 FORMAT(/148('=')/'||PRESSURE||L2 NORMS on the domain of VOLUME='
     %,G14.6 /148('=') )
C
C     AFFICHAGE DES NORMES L2 VITESSE et ERREURS
      IF( NOFOPR .GT. 0 .AND. NOFOVI .GT. 0 ) THEN
C
         IF( NDIM .EQ. 2 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10500) VOLUME
            ELSE
               WRITE(IMPRIM,20500) VOLUME
            ENDIF
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10501) VOLUME
            ELSE
               WRITE(IMPRIM,20501) VOLUME
            ENDIF
         ENDIF
C
         DO K=1,NBVECT
C           NORME L2 DE LA VITESSE ET ERREURS
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10515) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2VER(K),
     %            INTL2VER(K)/INTL2VIT(K) * 100D0
            ELSE
               WRITE(IMPRIM,20515) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2VER(K),
     %            INTL2VER(K)/INTL2VIT(K) * 100D0
            ENDIF
         ENDDO
C
10515    FORMAT('Cas',I5,' Temps=',G14.6,
     % ' ||Vitesse||=',G14.6,
     % ' ||Vit/Vol||=',G14.6,
     %'  ||ErrVit||=',G14.6,
     %'  ||ErrVit||/||Vit||=',G12.4,'%')

20515    FORMAT('Case',I5,' Time=',G14.6,
     % ' ||Velocity||=',G14.6,
     % ' ||Vel/Vol||=',G14.6,
     %'  ||VelErr||=',G14.6,
     %'  ||ErrVel||/||Vel||=',G12.4,'%')
         WRITE(IMPRIM,10600)
10600    FORMAT(148('=')/)
C
C        NORME L2 DE LA PRESSION ET ERREURS
         IF( NDIM .EQ. 2 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10502) VOLUME
            ELSE
               WRITE(IMPRIM,20502) VOLUME
            ENDIF
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10503) VOLUME
            ELSE
               WRITE(IMPRIM,20503) VOLUME
            ENDIF
         ENDIF
C
         DO K=1,NBVECT
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10514) K, TIMES(K),
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL,
     %            INTL2PER(K),
     %            INTL2PER(K)/INTL2PRE(K) * 100D0
            ELSE
               WRITE(IMPRIM,20514) K, TIMES(K),
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL,
     %            INTL2PER(K),
     %            INTL2PER(K)/INTL2PRE(K) * 100D0
            ENDIF
         ENDDO
C
10514 FORMAT('Cas',I5,' Temps=',G14.6,
     % ' ||Pre-PreMin||=',G14.6,
     % ' ||Pre-PreMin/Vol||=',G14.6,
     %'  ||ErrPre||=',G14.6,
     %'  ||PreErr||/||Pre||=',G12.4,'%')

20514 FORMAT('Case',I5,' Time=',G14.6,
     % ' ||Pre-PreMin||=',G14.6,
     % ' ||Pre-PreMin/Vol||=',G14.6,
     %'  ||PreErr||=',G14.6,
     %'  ||PreErr||/||Pre||=',G12.4,'%')
C
      ELSE IF( NOFOPR .LE. 0 .AND. NOFOVI .LE. 0 ) THEN
         DO K=1,NBVECT
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10511) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL
            ELSE
               WRITE(IMPRIM,20511) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL
            ENDIF
         ENDDO
C
10511 FORMAT('Cas',i5,' Temps=',g14.6,
     % ' ||Vitesse||=',G14.6,
     % ' ||Vit/Vol||=',G14.6,
     %'  ||Pre-PreMin||=',G14.6,
     %'  ||Pre-PreMin/Vol||=',G14.6)

20511 format('Case',I5,' Time=',g14.6,
     % ' ||Velocity||=',G14.6,
     % ' ||Vel/Vol||=',G14.6,
     %'  ||Pre-PreMin||=',G14.6,
     %'  ||Pre-PreMin/Vol||=',G14.6)
C
      ELSE IF( NOFOPR .LE. 0 .AND. NOFOVI .GT. 0 ) THEN
         DO K=1,NBVECT
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10512) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2VER(K),
     %            INTL2VIT(K)/INTL2VER(K) * 100D0,
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL
            ELSE
               WRITE(IMPRIM,20512) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2VER(K),
     %            INTL2VIT(K)/INTL2VER(K) * 100D0,
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL
            ENDIF
         ENDDO
C
10512 FORMAT('Cas',I5,' Temps=',G14.6,
     % ' ||Vitesse||=',G14.6,
     % ' ||Vit/Vol||=',G14.6,
     %'  ||ErrVit||=',G14.6,
     %'  ||ErrVit||/||Vit||=',G12.4,'%',
     %'  ||Pre-PreMin||=',G14.6,
     %'  ||Pre-PreMin/Vol||=',G14.6)

20512 FORMAT('Case',I5,' Time=',G14.6,
     % ' ||Velocity||=',G14.6,
     % ' ||Vel/Vol||=',G14.6,
     %'  ||VelErr||=',G14.6,
     %'  ||ErrVel||/||Vel||=',G12.4,'%',
     %'  ||Pre-PreMin||=',G14.6,
     %'  ||Pre-PreMin/Vol||=',G14.6)
C
      ELSE IF( NOFOPR .GT. 0 .AND. NOFOVI .LE. 0 ) THEN
         DO K=1,NBVECT
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10513) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL,
     %            INTL2PER(K),
     %            INTL2PRE(K)/INTL2PER(K) * 100D0
            ELSE
               WRITE(IMPRIM,20513) K, TIMES(K),
     %            INTL2VIT(K),
     %            INTL2VIT(K)/SQRVOL,
     %            INTL2PRE(K),
     %            INTL2PRE(K)/SQRVOL,
     %            INTL2PER(K),
     %            INTL2PRE(K)/INTL2PER(K) * 100D0
            ENDIF
         ENDDO
      ENDIF
C
10513 FORMAT('Cas',I5,' Temps=',G14.6,
     % ' ||Vitesse||=',G14.6,
     % ' ||Vit/Vol||=',G14.6,
     %'  ||Pre-PreMin||=',G14.6,
     %'  ||Pre-PreMin/Vol||=',G14.6,
     %'  ||ErrPre||=',G14.6,
     %'  ||PreErr||/||Pre||=',G12.4,'%')

20513 FORMAT('Case',I5,' Time=',G14.6,
     % ' ||Velocity||=',G14.6,
     % ' ||Vel/Vol||=',G14.6,
     %'  ||Pre-PreMin||=',G14.6,
     %'  ||Pre-PreMin/Vol||=',G14.6,
     %'  ||PreErr||=',G14.6,
     %'  ||PreErr||/||Pre||=',G12.4,'%')
C
      WRITE(IMPRIM,10600)
C
      RETURN
      END
