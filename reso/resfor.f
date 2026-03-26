      SUBROUTINE RESFOR( NBNOEU, NDIM, NTDL, NDSM, FORCE, KFORCE,
     %                   RESULF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA RESULTANTE DES NDSM FORCES
C -----
C
C ENTREE :
C --------
C NBNOEU : NOMBRE DE NOEUDS DU MAILLAGE
C NDIM   : NOMBRE DE COMPOSANTES DES FORCES
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C NDSM   : NOMBRE DE CAS DE FORCES
C FORCE  : FORCE(NTDL,NDSM) LES FORCES AUX NOEUDS DU MAILLAGE
C KFORCE : NOM 'FORCE ou SOURCE' SELON LE PROBLEME TRAITE
C
C SORTIE :
C --------
C RESULF : LES NDSM RESULTANTES DES FORCES DE L'OBJET COMPLET
C          RESULF(NDIM,NDSM)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JUILLET 1989
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  FORCE(NTDL,NDSM),RESULF(NDIM,NDSM)
      CHARACTER*(*)     KFORCE
C
C     LES COMPOSANTES DES RESULTANTES SONT MISES A ZERO
      CALL AZEROD( NDSM*NDIM , RESULF )
C
      IF( NTDL .NE. NBNOEU*NDIM ) THEN
          WRITE (KERR(MXLGER-1)(1:8),'(I8)') NTDL
          WRITE (KERR(MXLGER)(1:8)  ,'(I8)') NBNOEU*NDIM
          NBLGRC(NRERR) = 3
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'NOMBRE TOTAL D''INCONNUES        :'
     %            // KERR(MXLGER-1)(1:8)
             KERR(2) = 'NOMBRE DE COMPOSANTES DES '//KFORCE//'S :'
     %            // KERR(MXLGER)(1:8)
             KERR(3) = 'ERREUR: CES 2 VALEURS DOIVENT ETRE EGALES'
          ELSE
             KERR(1) = 'DEGREES of FREEDOM TOTAL NUMBER:'
     %            // KERR(MXLGER-1)(1:8)
             KERR(2) = 'COMPONENT NUMBER of '//KFORCE//'S :'
     %            // KERR(MXLGER)(1:8)
             KERR(3) = 'ERROR: These 2 VALUES MUST BE EQUALS'
          ENDIF
          CALL LEREUR
          RETURN
      ENDIF
C
      DO 30 J=0,NBNOEU-1
         DO 20 I=1,NDIM
            NODL = J*NDIM + I
            DO 10 K=1,NDSM
               RESULF(I,K) = RESULF(I,K) + FORCE(NODL,K)
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
C
C     AFFICHAGE DES COMPOSANTES DE LA RESULTANTE
C     ==========================================
      IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,10030) ((K,KFORCE,J,RESULF(J,K),K=1,NDSM),J=1,NDIM)
      ELSE
        WRITE(IMPRIM,20030) ((K,KFORCE,J,RESULF(J,K),K=1,NDSM),J=1,NDIM)
      ENDIF
10030 FORMAT(' CAS',I3,' RESULTANTE des ',A,'S : COMPOSANTE',
     %I2,' =',G15.7)
20030 FORMAT(' CASE',I3,' RESULTANT of ',A,'S : COMPONENT',
     %I2,' =',G15.7)
      RETURN
      END
