      SUBROUTINE NUMTMS( KNOMTS , NOTAMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETROUVER LE NUMERO DE TMS D'UN TMS CONNU PAR SON NOM
C -----
C ENTREE :
C --------
C KNOMTS : NOM DU TABLEAU TMS
C
C SORTIES :
C ---------
C NOTAMS : >0 NUMERO DU TMS S'IL EXISTE
C          =0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      CHARACTER*(*)     KNOMTS
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSSFTA / MSSF(28),NTADAM
C
      N1  = 1
      N2  = LEN( KNOMTS )
      NT1 = NTADAM
      IF( KNOMTS(N1:N1) .EQ. '~' ) N1 = 2
      IF( KNOMTS(N1:N1) .EQ. '>' ) N1 = N1 + 1
C
C     RECHERCHE DU SEPARATEUR > JUSQU'A NE PLUS EN TROUVER
 10   I = INDEX( KNOMTS(N1:N2) , '>' )
      IF( I .GT. 0 ) THEN
C        OUVERTURE DU LEXIQUE
         CALL LXLXOU( NT1 , KNOMTS(NT1:I-1) , NT2 , MN2 )
         IF( NT2 .GT. 0 ) THEN
C           PASSAGE AU SUIVANT
            N1  = I + 1
            NT1 = NT2
            GOTO 10
         ENDIF
      ENDIF
C
C     LE NOM FEUILLE DU TMS EST KNOMTS(N1:N2)
      CALL LXTSOU( NT1 , KNOMTS(N1:N2) , NOTAMS , MN2 )
      END
