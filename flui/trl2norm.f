      SUBROUTINE TRL2NORM( NDIM, NBVECT, TIMES,  VOLUME, INTVIT, INTPRE,
     %                     NOFOVI, INTVER, NOFOPR, INTPER, TABAUX  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES COURBES de la NORME L2 MOYENNE en FONCTION DU TEMPS de
C ----- ||Vitesse|| , PRESSION-PRESSION Min,
C       ||Vitesse Exacte-Calculee||, ||Vitesse Exacte-Calculee||/||Vitesse||
C        PRESSION Exacte-Calculee,   PRESSION Exacte-Calculee/PRESSION
C
C       LES DONNEES ont ete CALCULEES PAR fl2norm.f
C
C ENTREES:
C --------
C NDIM   : DIMENSION 2 ou 3 DE L'ESPACE
C NBVECT : NOMBRE D'ELEMENTS DE CHAQUE TABLEAU
C TIMES  : TEMPS DES NBVECT ELEMENTS
C VOLUME : VOLUME DU MAILLAGE
C          Volume= Som  Delta(e)
C                  e dans E
C                                                          (Vk1C)
C INTVIT:  Som Jacobien Som (Vk1C,...,VknC) [int Pi Pj DX] ( ...) / Volume
C          e de E       k=1,...,d                          (VknC)
C                                                    (P1C)
C INTPRE:  Som Jacobien (P1C,...,PmC) [int Pi Pj DX] (...) / Volume
C          e de E                                    (PmC)
C
C NOFOVI : NUMERO DE LA FONCTION VITESSE_EXACTE(t,x,y,z)
C          0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C                                                                    (Vk1E-Vk1C)
C INTVER:  Som Jacobien Som (Vk1E-Vk1C,...,VknE-VknC) [int Pi Pj DX] (    ...  )
C          e de E       k=1,...,d                                    (VknE-VknC)
C
C NOFOPR : NUMERO DE LA FONCTION PRESSION_EXACTE(t,x,y,z)
C          0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C                                                            (P1E-P1C)
C INTPER:  Som Jacobien (P1E-P1C,...,PmE-PmC) [int Pi Pj DX] (  ...  )
C          e de E                                            (PmE-PmC)
C
C LA NORME L2 DE L'ERREUR EST CALCULEE SEULEMENT SI LA FONCTION
C PRESSION_EXACTE(t,x,y,z) (NOFOPR>0)
C et/ou VITESSE_EXACTE(t,x,y,z,nocomp) (NOFOVI>0)
C SONT FOURNIES PAR L'UTILISATEUR
C
C TABAUX : TABLEAU AUXILIAIRE DE REELS DOUBLE POUR TRACER LES VALEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  TEXAS A & M UNIVERSITY at QATAR Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
C
      DOUBLE PRECISION  VOLUME, SQRVOL,
     %                  INTPRE(NBVECT), INTVIT(NBVECT),
     %                  INTPER(NBVECT), INTVER(NBVECT), TABAUX(NBVECT)
      REAL              TIMES(NBVECT)
      CHARACTER*7       NMVOLSUR
      INTRINSIC         SQRT
C
      IF( INTERA .LE. 0 ) RETURN
C
      IF( NDIM .EQ. 2 ) THEN
         NMVOLSUR = 'Surface'
         NB = 7
      ELSE
         NMVOLSUR = 'Volume '
         NB = 6
      ENDIF
      SQRVOL = SQRT( VOLUME )
C
C     TRACE DE ||Vitesse/Volume||L2
C     DIVISION PAR VOLUME
      DO K=1,NBVECT
         TABAUX(K) = INTVIT(K) / SQRVOL
      ENDDO
      CALL TRTABLE( NBVECT, TIMES, 1,1,TABAUX, 'L2MeanVelocity', NCNOIR,
     %              'Time',
     %              '||Velocity/'//NMVOLSUR(1:NB)//'||L2',
     %              '||Velocity/'//NMVOLSUR(1:NB)//'||L2(Time)',
     %              ' ',
     %              ' ' )
C
      IF( NOFOVI .GT. 0 ) THEN
C
C        TRACE DE || Exact|Velocity|-Computed|Velocity| ||L2 / Volume
C        DIVISION PAR VOLUME
         DO K=1,NBVECT
            TABAUX(K) = INTVER(K) / SQRVOL
         ENDDO
         CALL TRTABLE( NBVECT, TIMES, 1,1,TABAUX, 'L2VelErr',
     %                 NCROSE, 'Time',
     %   '||Exact|Velocity|-Computed|Velocity|||L2/'//NMVOLSUR(1:NB),
     %   '||Exact|Velocity|-Computed|Velocity|||L2/'//NMVOLSUR(1:NB)//
     %'(Time)',
     %      ' ',
     %      ' ' )
C
C        TRACE DE ||Exact|Velocity|-Computed|Velocity|||L2/||Computed|Velocity||
         DO K=1,NBVECT
            TABAUX(K) = INTVER(K) / INTVIT(K) * 100D0
         ENDDO
         CALL TRTABLE( NBVECT, TIMES, 1,1,TABAUX, 'L2VelErrRel',
     %                 NCJAUN, 'Time',
     %'||Exact|Vel|-Comp|Vel|||L2/||Comp|Vel|||L2 %',
     %'||Exact|Velocity|-Computed|Velocity|||L2 / ||Computed|Velocity||L
     %2 % (Time)',
     %      ' ',
     %      ' ' )
C
      ENDIF
C
C     TRACE DE ||Pression-PressionMin||L2 / volume
C     DIVISION PAR VOLUME
      DO K=1,NBVECT
         TABAUX(K) = INTPRE(K) / SQRVOL
      ENDDO
      CALL TRTABLE( NBVECT, TIMES, 1,1,TABAUX, 'L2MeanPressure',
     %              NCNOIR, 'Time',
     %             '||(Pressure-PreMin)/'//NMVOLSUR(1:NB)//'||L2',
     %             '||(Pressure-PreMin)/'//NMVOLSUR(1:NB)//'||L2(Time)',
     %             ' ',
     %             ' ' )
C
      IF( NOFOPR .GT. 0 ) THEN
C
C        TRACE DE ||ExactPressure-ComputedPressure||L2 / Volume
C        DIVISION PAR VOLUME
         DO K=1,NBVECT
            TABAUX(K) = INTPER(K) / SQRVOL
         ENDDO
         CALL TRTABLE( NBVECT, TIMES, 1,1,TABAUX, 'L2PreErr',
     %                 NCROSE, 'Time',
     %         '||ExactPressure-ComputedPressure||L2/'//NMVOLSUR(1:NB),
     %         '||ExactPressure-ComputedPressure||L2/'//NMVOLSUR(1:NB)//
     %'(Time)',
     %          ' ',
     %          ' ' )
C
C        TRACE DE ||ExactPressure-ComputedPressure||L2 / ||ComputedPressure||L2
         DO K=1,NBVECT
            TABAUX(K) = INTPER(K) / INTPRE(K) * 100D0
         ENDDO
         CALL TRTABLE( NBVECT, TIMES, 1,1,TABAUX, 'L2PreErrRel',
     %                 NCJAUN, 'Time',
     %'||ExactPres-CompPres||L2/||CompPres||L2 %',
     %'||ExactPressure-ComputedPressure||L2 / ||ComputedPressure||L2 % (
     %Time)',
     %          ' ',
     %          ' ' )
C
      ENDIF
C
      RETURN
      END
