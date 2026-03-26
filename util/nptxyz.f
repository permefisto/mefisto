      SUBROUTINE NPTXYZ( NUPOIN, MNXYZ, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETROUVER LES 3 COORDONNEES D'UN POINT A PARTIR DE SON NUMERO
C ----- DE POINT DANS LE LEXIQUE 'POINT'
C
C ENTREES :
C ---------
C NUPOIN : NUMERO DU POINT DANS LE LEXIQUE POINT
C
C SORTIES :
C ---------
C MNXYZ  : L'ADRESSE MCN DES 3 COORDONNEES DU POINT NUPOIN
C IERR   : 0 SI XYZ INITIALISE
C          1 SI LE POINT OU 2 SI LE TABLEAU XYZSOMMET N'EXISTE PAS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1993
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
C
      CALL LXNLOU( NTPOIN, NUPOIN, NTP, MNP )
      IF( NTP .LE. 0 ) THEN
         WRITE(KERR(MXLGER)(1:10),'(I10)') NUPOIN
         NBLGRC(NRERR) = 1
         KERR(1)='ERREUR: POINT INCONNU NUMERO '// KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
      CALL LXTSOU( NTP, 'XYZSOMMET', NTXYZ, MNXYZ )
      IF( NTXYZ .LE. 0 ) THEN
         CALL NMOBNU( KNM , NUOB , KERR(3) )
         NBLGRC(NRERR) = 2
         KERR(1) = 'POINT : ' // KERR(3)(1:24)
         KERR(2) = 'ERREUR: XYZ INCONNU'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LE DECALAGE POUR ARRIVER SUR LES COORDONNEES
      MNXYZ = MNXYZ + WYZSOM
      IERR  = 0
C
      RETURN
      END
