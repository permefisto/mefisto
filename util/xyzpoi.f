      SUBROUTINE XYZPOI( NUPOIN , MNDUPT , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES 3 COORDONNEES DU POINT DE NUMERO NUPOIN
C -----
C ENTREES:
C --------
C NUPOIN : NUMERO DU POINT DANS LE LEXIQUE DES POINTS
C
C SORTIES:
C --------
C MNDUPT : ADRESSE MCN DU DEBUT DES COORDONNEES DU POINT
C          DANS LE TABLEAU 'XYZSOMMET'
C IERR   : 0   PAS D'ERREUR
C          1   POINT INCONNU OU SANS TABLEAU 'SOMMET'
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  UPMC   OCTOBRE 1988
C......................................................................
      IMPLICIT          INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
C
C     OUVERTURE DU POINT DANS LE LEXIQUE DES POINTS
      CALL LXNLOU( NTPOIN , NUPOIN , NTLXPO , MNLXPO )
      IF( NTLXPO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NUPOIN
         KERR(1) = ' XYZPOI:POINT INCONNU '//KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     OUVERTURE DU TABLEAU 'XYZSOMMET' DU POINT
      CALL LXTSOU( NTLXPO , 'XYZSOMMET' , NTSOMM , MNSOMM )
      IF( NTSOMM .LE. 0 ) THEN
          NBLGRC(NRERR) = 1
          WRITE(KERR(MXLGER)(1:10),'(I10)') NUPOIN
          KERR(1) =  ' XYZPOI:POINT '//KERR(MXLGER)(1:10)
          CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     L'ADRESSE DES COORDONNEES DU POINT
      MNDUPT = MNSOMM + WYZSOM
 9999 RETURN
      END
