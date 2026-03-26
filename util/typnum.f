      CHARACTER*9 FUNCTION TYPNUM( NUMTYV )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LES 9 CARACTERES DU NOM DE NUMERO DE TYPE NUMTYV
C -----
C
C ENTREE :
C --------
C NUMTYV : NUMERO DU TYPE ( ENTIER COMPRIS ENTRE 1 ET 9 )
C
C SORTIE :
C --------
C TYPNUM : CHAINE DE 9 CARACTERES
C          1=>'LOGIQUE'   2 =>'CARACTERE' 3=>'ENTIER/2'   4=>'ENTIER'
C          5=>'REEL'      6 =>'REEL2'     7=>'REEL4'      8=>'COMPLEXE'
C          9=>'COMPLEXE2' 10=>'MOTS'     11=>'^LEXIQUE'  12=>'XYZ'
C         13=>'TYPEOBJET'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  NOVEMBRE 1983
C.......................................................................
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
      IF( NUMTYV .LE. 0 .OR. NUMTYV .GT. NBTYPV ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUMTYV
         KERR(1) = 'TYPNUM :TYPE ERRONE DE VARIABLES='
     %           //KERR(MXLGER)(1:4)
         CALL LEREUR
         TYPNUM = 'INCONNU'
      ELSE
         TYPNUM = KTYPES( NUMTYV )
      ENDIF
      END
