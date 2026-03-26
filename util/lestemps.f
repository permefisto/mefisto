      SUBROUTINE LESTEMPS( KNOMOB, MNVECT, NBVECT, MNTIMES, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE DANS UN TABLEAU LES TEMPS DU CALCUL DES NBVECT VECTEURS
C -----    D'UN TMS VECTEUR"...
C          SI LES TEMPS NE SONT PAS STOCKES, ALORS LE NUMERO DU VECTEUR
C          EST STOCKE
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"...
C
C SORTIES:
C --------
C NBVECT : NOMBRE DE VECTEURS ET DE TEMPS
C MNTIMES: ADRESSE MCN DU TABLEAU DES TEMPS DE CALCUL
C IERR   : 0 SI PAS D'ERREUR, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/a___vecteur.inc"
      include"./incl/ctemps.inc"
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*(*)     KNOMOB
C
      IERR    = 0
      MNTIMES = 0
C
C     NOMBRE DE VECTEUR"...  STOCKES
      NBVECT = MCN( MNVECT + WBVECT )
      IF( NBVECT .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'AUCUN VECTEUR SOLUTION CALCULE'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'NO COMPUTED SOLUTION VECTOR'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9000
      ENDIF
C
C     NOMBRE DE COMPOSANTES DE CHACUN DES VECTEURS
      NTDL = MCN( MNVECT + WBCOVE )
C
C     ADRESSE DES TEMPS DES VECTEURS SAUVEGARDES OU NON
      CALL TNMCDC( 'REEL', NBVECT, MNTIMES )
      IF( MNTIMES .LE. 0 ) THEN
         IERR = 1
         GOTO 9000
      ENDIF
C
      IF( MCN(MNVECT+WBCPIN ) .GT. 0 ) THEN
C        LES TEMPS SONT STOCKES
         MOREE2 = MOTVAR(6)
         MNTIME = MNVECT + WECTEU + NTDL * NBVECT * MOREE2 -1
         DO 11 I=1,NBVECT
            RMCN(MNTIMES-1+I) = RMCN(MNTIME+I)
 11      CONTINUE
      ELSE
C        PAS DE TEMPS STOCKES DANS LE TMS VECTEUR
         DO 12 I=1,NBVECT
            RMCN(MNTIMES-1+I) = I
 12      CONTINUE
         MNTIME = 0
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'Il existe',NBVECT,' vecteurs de',
     %   NTDL,' composantes stockes aux TEMPS:'
      ELSE
         WRITE(IMPRIM,*)'It exists',NBVECT,' vectors of ',
     %   NTDL,' components stored at TIMES:'
      ENDIF
      WRITE(IMPRIM,10012) (k,RMCN(MNTIMES-1+k),k=1,nbvect)
10012 FORMAT(5(I5,':',G14.6,'  '))
C
 9000 RETURN
      END
