      SUBROUTINE RECTEF( NUR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFACER LE RECTANGLE DE NUMERO NURECT
C -----          VERSION xvue
C
C ENTREE :
C --------
C NUR : NUMERO DU RECTANGLE A EFFACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1994
C23456---------------------------------------------------------------012
      include"./incl/epombr.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
C
      if ( nur .EQ. nrmenu ) goto 9999
      IF( INTERA .LT. 3 ) RETURN
      IF( DXRECT( NUR ) .LE. 0 .OR. DYRECT( NUR ) .LE. 0 ) GOTO 9999
C
C     SI LE MODE POSTSCRIPT EST EN COURS
      IF ( LASOPS .NE. 0 ) THEN
        IF ( LASOPS .EQ. 1 ) THEN
C         MISE EN SUSPEND DU POSTSCRIPT
          LASOPS = -1
        ELSE
          IF ( LASOPS .EQ. 2 ) THEN
C           EFFACEMENT DU MENU CORRESPONDANT
            LASOPS = -3 - NUR
            LASOPS = MAX(-10,LASOPS)
          ELSE
CCC            NBLGRC(NRERR) = 2
CCC            KERR(1) = 'RECTEF: MAUVAISE VALEUR DE LASOPS'
CCC            KERR(2) = '        ARRET DU POSTSCRIPT'
CCC            CALL LEREUR
CCC            LASOPS = 0
             LASOPS = -1
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     TRACE DU REMPLISSAGE DU RECTANGLE EN COULEUR DU FOND
      CALL XVCOULEUR( NCOFON )
      CALL XVFRECTANGLE( XRECT(NUR)-NBEPO4-1,  YRECT(NUR)-NBEPO4-1,
     %                   DXRECT(NUR)+NBEPO3+2, DYRECT(NUR)+NBEPO3+3 )
C
C     SI LE MODE POSTSCRIPT ETAIT EN COURS
      IF ( LASOPS .NE. 0 ) THEN
C       RETOUR AU MODE POSTSCRIPT NORMAL
        IF ( ABS(LASOPS) .EQ. 1 ) THEN
          LASOPS = 1
        ELSE
          LASOPS = 2
        ENDIF
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     REMISE A ZERO DU RECTANGLE
 9999 MDLGRC( NUR ) = 0
      IF( NUR .NE. NRHIST ) NBLGRC( NUR ) = 0
      DXRECT( NUR ) = 0
      DYRECT( NUR ) = 0
      END
