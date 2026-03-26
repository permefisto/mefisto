      SUBROUTINE NLSEITER( MODE, TEMPS, ITERm, NBNOMA,
     %                     UC0,  UN0,   UC1,   UN1,
     %                     MODUMX, Testm )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES NORMES ABSOLUES ET RELATIVES DE U0, U1
C -----    ET LA VALEUR POUR TESTER LA FIN DES ITERATIONS m
C
C ENTREES:
C --------
C MODE   : 0 => UTILISE UC*(NBNOMA,2)
C          1 => UTILISE UN*(2,NBNOMA)
C TEMPS  : INSTANT DU CALCUL
C ITERm  : NUMERO DE L'ITERATION m
C NBNOMA : NOMBRE DE COMPOSANTES DES VECTEURS UC0, UC1(NBNOMA,2)
C          ou NOMBRE DE NOEUDS DU MAILLAGE

C UC0    : VECTEUR(NBNOMA,2) DE L'ITERATION m
C UN0    : VECTEUR(2,NBNOMA) DE L'ITERATION m

C UC1    : VECTEUR(NBNOMA,2) DE L'ITERATION m+1
C UN1    : VECTEUR(2,NBNOMA) DE L'ITERATION m+1
C
C SORTIES:
C --------
C MODUMX : Max | U1(tn+1,m+1,Noeud) | aux NOEUDS du MAILLAGE
C Testm  : Som||Um+1|(N)-|Um|(N)|/Som|Um+1|(N) aux NOEUDS N du MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Octobre 2013
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      REAL              TEMPS
      INTEGER           MODE, ITERm, NBNOMA, I, LECTEU, IMPRIM, NUNITE
      DOUBLE PRECISION  UC0(NBNOMA,2), UC1(NBNOMA,2),
     %                  UN0(2,NBNOMA), UN1(2,NBNOMA), Testm
      DOUBLE PRECISION  DUM0R, DUM0I, DUM0, DUM1R, DUM1I, DUM1, DIF,
     %                  NORMDF, NORMUM, NOMXDF, MODUMX
C
C     CALCUL DE || Un+1m+1 - Un+1m || et ||Un+1m+1|| L1 et MAX
C     =================================================================
      NORMUM = 0D0
      MODUMX = 0D0
      NORMDF = 0D0
      NOMXDF = 0D0

      DO I=1,NBNOMA

         IF( MODE .EQ. 0 ) THEN
C
C           PARTIES REELLE, IMAGINAIRE et MODULE DE U0(tn+1,m,Noeud I)
            DUM0R = UC0(I,1)
            DUM0I = UC0(I,2)
C
C           PARTIES REELLE, IMAGINAIRE et MODULE DE U1(tn+1,m+1,Noeud I)
            DUM1R = UC1(I,1)
            DUM1I = UC1(I,2)

         ELSE
C
C           PARTIES REELLE, IMAGINAIRE et MODULE DE U0(tn+1,m,Noeud I)
            DUM0R = UN0(1,I)
            DUM0I = UN0(2,I)
C
C           PARTIES REELLE, IMAGINAIRE et MODULE DE U1(tn+1,m+1,Noeud I)
            DUM1R = UN1(1,I)
            DUM1I = UN1(2,I)

         ENDIF
C
C        LES MODULES
         DUM0 = SQRT( DUM0R**2 + DUM0I**2 )
         DUM1 = SQRT( DUM1R**2 + DUM1I**2 )
C
C        Somme sur I des | U1(tn+1,m+1,Noeud I) |
         NORMUM = NORMUM + DUM1
C
C        Max sur I des | U1(tn+1,m+1,Noeud I) |
         IF( DUM1 .GT. MODUMX ) MODUMX = DUM1
C
C        | U1(tn+1,m+1,Noeud I) - U0(tn+1,m,Noeud I) | au NOEUD I
         DIF = ABS( DUM1 - DUM0 )
C
C        Somme sur I des | U1(tn+1,m+1,Noeud I) - U0(tn+1,m,Noeud I) |
         NORMDF = NORMDF + DIF
C
C        Max sur I des | U1(tn+1,m+1,Noeud I) - U0(tn+1,m,Noeud I) |
         IF( DIF .GT. NOMXDF ) NOMXDF = DIF
C
      ENDDO
C
C     NORMUM = MOYENNE de la Somme sur I des | U1(tn+1,m+1,Noeud I) |
      NORMUM = NORMUM / NBNOMA
C
C     NORMDF = MOYENNE de la Somme sur I des | |U1(tn+1,m+1,Noeud)| - |U0(tn+1,m,Noeud)|
      NORMDF = NORMDF / NBNOMA
C
C     Testm = TEST D'ARRET DE L'ITERATION m =  Som| |Um+1|-|Um| | / Som |Um+1|
      IF( NORMUM .EQ. 0D0 ) THEN
         Testm = 0.0
      ELSE
         Testm = NORMDF / NORMUM
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10130) TEMPS,  ITERm,
     %                       MODUMX, NOMXDF, NOMXDF/MODUMX,
     %                       NORMUM, NORMDF, Testm
      ELSE
         WRITE(IMPRIM,20130) TEMPS,  ITERm,
     %                       MODUMX, NOMXDF, NOMXDF/MODUMX,
     %                       NORMUM, NORMDF, Testm
      ENDIF
C
10130 FORMAT('Au TEMPS',G14.6,' ITER m=',I3,
     %T35, 'Max|Um+1(N)|     =',G14.6,
     %T69, 'Max|Um+1(N)-Um(N)|      =',G14.6,
     %T109,'Max|Um+1-Um|/Max|Um+1|    =',G14.6/
     %T35, 'Som|Um+1|/NbNoeud=',G14.6,
     %T69, 'Som||Um+1|-|Um||/NbNoeud=',G14.6,
     %T109,'Som||Um+1|-|Um||/Som|Um+1|=',G14.6)
C
20130 FORMAT('At TIME',G14.6,' ITER m=',I3,
     %T35, 'Max|Um+1(N)|    =',G14.6,
     %T69, 'Max|Um+1(N)-Um(N)|     =',G14.6,
     %T109,'Max|Um+1-Um|/Max|Um+1|    =',G14.6/
     %T35, 'Sum|Um+1|/NodeNb=',G14.6,
     %T69, 'Sum||Um+1|-|Um||/NodeNb=',G14.6,
     %T109,'Sum||Um+1|-|Um||/Sum|Um+1|=',G14.6)

      RETURN
      END
