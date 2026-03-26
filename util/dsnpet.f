      SUBROUTINE DSNPET( NTLXOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE LES TABLEAUX NOEUDS POINTS ELEMENTS TOPOLOGIE
C -----    DANS LE LEXIQUE NTLXOB
C
C ENTREE :
C --------
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C2345X7..............................................................012
      IMPLICIT          INTEGER (W)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
C
      CHARACTER*4       NOM4,CHARX
C
C     LE TABLEAU TOPOLOGIE
      CALL LXTSOU( NTLXOB , 'TOPOLOGIE' , NTTOPO , MNTOPO )
      IF( NTTOPO .LE. 0 ) RETURN
C
C     RECUPERATION DU TABLEAU NOEUDS
      CALL LXTSOU( NTLXOB , 'XYZNOEUD'  , NTNOEU , MNNOEU )
      IF( NTNOEU .GT. 0 ) THEN
         CALL LXTSDS( NTLXOB , 'XYZNOEUD' )
      ENDIF
C
C     RECUPERATION DU TABLEAU POINTS
      CALL LXTSOU( NTLXOB , 'XYZPOINT'  , NTPOGE , MNPOGE )
      IF( NTPOGE .GT. 0 ) THEN
         CALL LXTSDS( NTLXOB , 'XYZPOINT' )
      ENDIF
C
C     LES ELEMENTS SONT RETROUVES ET DETRUITS
      NBTYEL = MCN(MNTOPO + WBTYEL )
      DO  10 I=0,NBTYEL-1
C        LE NOM DU TABLEAU ELEMENTS"
         NOM4 = CHARX( MCN(MNTOPO+WMTYEL+I) )
         CALL LXTSOU( NTLXOB , 'ELEMENTS"'//NOM4 , NT , MN )
         IF( NT .GT. 0 ) THEN
            CALL LXTSDS( NTLXOB , 'ELEMENTS"'//NOM4 )
         ENDIF
 10   CONTINUE
C
C     DESTRUCTION DE TOPOLOGIE
      CALL LXTSDS( NTLXOB , 'TOPOLOGIE' )
      END
