      SUBROUTINE TTNODS( MNTATA,NOTATA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETRUIRE LE TABLEAU NUMERIQUE NOTATA DU TABLEAU DE TABLEAUX
C ----
C ENTREES :
C ---------
C MNTATA : ADRESSE MCN DU TABLEAU DES NUMEROS DE TAMS DES TABLEAUX
C NOTATA : NUMERO DU TABLEAU A DECLARER (1 A MXTATA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  OCTOBRE 1985
C.......................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     LE NUMERO NOTATA EST IL >0 ET < MXTATA ?
      MXTATA = MCN( MNTATA + 3 )
      IF( NOTATA .LE. 0 .OR. NOTATA .GT. MXTATA ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' NO TABLEAU INCORRECT'
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU EST IL DECLARE ?
      MN = MNTATA + 3 + NOTATA
      IF( MCN(MN) .NE. 0 ) THEN
C        OUI . LE TABLEAU EST DETRUIT
         CALL TAMSOU( ABS( MCN(MN) ) , MCNOTA )
         MCN( MN ) = 0
      ENDIF
      END
