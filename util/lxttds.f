      SUBROUTINE LXTTDS( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETRUIRE LE NOM KNOM TABLEAU DE TABLEAUX DU LEXIQUE NTLX
C ----- ET TOUS SES TABLEAUX
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C KNOM   : CHAINE DE CARACTERES NOM A DETRUIRE DU LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1985
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     KNOM
      COMMON / UNITES / LECTEU, IMPRIM,NUNITE(30)
C
C     CALCUL DU NUMERO DE KNOM DANS LE LEXIQUE ET DE SON CHAINAGE AMONT
      CALL LXNMNO( NTLX , KNOM , NONOM , MNLX )
C
C     OUVERTURE DU TABLEAU DE TABLEAUX
      IF( NONOM .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(2) = 'NOM INCONNU DANS LE LEXIQUE'
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE NUMERO DE TABLEAU MS DU TABLEAU DE TABLEAUX
      M1LX   = MCN( MNLX )
      NBENNM = MCN( MNLX + 2 )
      MN     = MNLX + M1LX * NONOM
      NTTATA = MCN( MN + NBENNM + 2 )
      CALL TAMSOU( NTTATA , MNTATA )
C
C     DESTRUCTION DES TABLEAUX DU TATA
      DO 10 N=4,MCN(MNTATA+3)+3
         NT = MCN( MNTATA + N )
         IF( NT .EQ. 0 ) GOTO 10
         CALL TAMSDS( ABS( NT ) )
 10   CONTINUE
C
C     DESTRUCTION DU TABLEAU DE TABLEAUX
      CALL LXNMDS( NTLX , KNOM )
      END
