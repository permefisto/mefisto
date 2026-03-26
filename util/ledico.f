      SUBROUTINE LEDICO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LE FICHIER ./td/dicotd DANS LE TABLEAU DICOTD
C ----- CORRESPONDANCE NUMERO ET NOM DES TABLEAUX DESCRIPTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      include"./incl/homdir.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*160     KNOM
C
C     RECHERCHE D'UN NUMERO D'UNITE
      CALL TRUNIT( NF )
C
C     LE NOM DU FICHIER
      KNOM = HOMDIR // '/td/dicotd'
      I    = NUDCNB( KNOM )
C     OUVERTURE DU FICHIER SUPPORT DU TD
      OPEN( FILE=KNOM(1:I), UNIT=NF, IOSTAT=J )
      IF( J .NE. 0 ) THEN
          NBLGRC(NRERR) = 1
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'PROBLEME pour OUVRIR le fichier ./td/dicotd'
          ELSE
             KERR(1) = 'PROBLEM to open the file ./td/dicotd'
          ENDIF
          CALL LEREUR
          RETURN
      ENDIF
C
C     LECTURE DU TABLEAU DICOTD
      DO 10 I=1,NBTD
         READ(NF,'(A48)',IOSTAT=J) DICOTD(I)
         IF( J .NE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='ERREUR en LECTURE de ./td/dicotd sur la ligne'
               KERR(2)=DICOTD(I)
            ELSE
               KERR(1)='ERROR during READ of ./td/dicotd on the line'
               KERR(2)=DICOTD(I)
            ENDIF
            CALL LEREUR
            RETURN
         ENDIF
C        LE NOM DICOTD(I) EST ECRIT EN minuscules
         CALL MINUSC( DICOTD(I) )
 10   CONTINUE
C
C     FERMETURE DU FICHIER SUPPORT DU TD
      CLOSE( UNIT=NF )
      RETURN
      END
