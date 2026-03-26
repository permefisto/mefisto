      SUBROUTINE AFEETH( NDSM,   NBTTEF, ERTHEF,
     %                   ESERTH, H1TEMP, EEH1TH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES ESTIMATEURS D'ERREUR SUR LE MAILLAGE TOTAL
C -----
C ENTREES :
C ---------
C NDSM   : NOMBRE DE CAS TRAITES
C NBTTEF : NOMBRE TOTAL D'EF DU MAILLAGE
C ERTHEF : ERREUR L2 A POSTERIORI SUR CHAQUE EF
C ESERTH : ESTIMATEUR D'ERREUR SUR LE MAILLAGE TOTAL POUR CHAQUE CAS
C H1TEMP : NORME H1 DE LA TEMPERATURE SUR LE MAILLAGE TOTAL POUR CHAQUE CAS
C EEH1TH : SOMME DES CARRES DES QUOTIENTS ESTIMATEUR
C          ERREUR / NORME H1 TEMPERATURE SUR LES EF + RACINE CARREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1995
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  ESERTH(1:NDSM),
     %                  H1TEMP(1:NDSM),
     %                  EEH1TH(1:NDSM), E1,
     %                  ERTHEF(1:NBTTEF,1:NDSM)
C
10000 FORMAT(1X,62('-'))
      WRITE(IMPRIM,*)
C
C     AFFICHAGE DES NORMES DES ESTIMATEURS D'ERREUR
C     =============================================
      DO 20 N=1,NDSM
C        AFFICHAGE DES ESTIMATEURS D'ERREUR POUR LE CAS N
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10010) N
         ELSE
            WRITE(IMPRIM,20010) N
         ENDIF
         WRITE(IMPRIM,10000)
         IF( ABS(H1TEMP(N)) .LE. 1D-28 ) THEN
C           POUR EVITER LA DIVISION PAR ZERO!
            E1 = 0D0
         ELSE
            E1 = ESERTH(N) / H1TEMP(N)
         ENDIF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10020) ESERTH(N), H1TEMP(N), E1, EEH1TH(N)
         ELSE
            WRITE(IMPRIM,20020) ESERTH(N), H1TEMP(N), E1, EEH1TH(N)
         ENDIF
 20   CONTINUE
C
10010 FORMAT(' CAS',I4)
20010 FORMAT(' CASE',I4)
10020 FORMAT(' | NORME L2 des ESTIMATEURS D''ERREUR    =',G16.7,T63,'|'/
     %       ' | NORME H1 de la TEMPERATURE sur les EF=', G16.7,T63,'|'/
     %       ' | NORME L2 ESTIMATEUR / NORME H1 TEMP  =', G16.7,T63,'|'/
     %       ' | MAX ESTIMATEUR 1 EF / NORME TEMP EFs =', G16.7,T63,'|')
20020 FORMAT(' | L2 NORM of the ERROR ESTIMATOR       =', G16.7,T63,'|'/
     %       ' | H1 NORM of the TEMPERATURE on FEs    =', G16.7,T63,'|'/
     %       ' | L2 NORM ESTIMATOR / H1 NORM TEMP     =', G16.7,T63,'|'/
     %       ' | MAX ESTIMATOR on 1FE / NORM TEMP FEs =', G16.7,T63,'|')
C
      WRITE(IMPRIM,10000)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10030)
      ELSE
         WRITE(IMPRIM,20030)
      ENDIF
C
10030 FORMAT(' ou NORME L2 (u) = RACINE CARREE de la SOMME des CARRES de
     % u sur tous les EF'/)
20030 FORMAT(' where L2 NORM(u) = SQUARE ROOT of SUM of SQUARE u on all
     %FE'/)
C
C
C     AFFICHAGE DE L'ESTIMATEUR D'ERREUR POUR CHAQUE EF
C     =================================================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10040)
         DO 50 I=1,MIN(10,NBTTEF)
            WRITE(IMPRIM,10050) (I,ERTHEF(I,N),N=1,NDSM)
 50      CONTINUE
      ELSE
         WRITE(IMPRIM,20040)
         DO 60 I=1,MIN(10,NBTTEF)
            WRITE(IMPRIM,20050) (I,ERTHEF(I,N),N=1,NDSM)
 60      CONTINUE
      ENDIF
C
10040 FORMAT(/' L''ESTIMATEUR d''ERREUR pour 10-ers EF du MAILLAGE :'/
     &       1X,80('='))
10050 FORMAT( ' EF',I6,' : ESTIMATEUR ERREUR =',4G16.7)
C
20040 FORMAT(/' The ERROR ESTIMATOR for the 10 first FE :'/
     &       1X,80('='))
20050 FORMAT( ' FE',I6,' : ERROR ESTIMATOR =',4G16.7)
C
      RETURN
      END
