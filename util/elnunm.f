      SUBROUTINE ELNUNM( NO, KNOMEL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   FOURNIR SELON LE NUMERO NO LE NOM DE L EF (KNOMEL(1),(2))
C -----
C
C ENTREE :
C --------
C NO     : NO DE L EF DANS LES SP UTILITAIRES
C
C SORTIE :
C --------
C KNOMEL : NOM SUR 2 MOTS DE 4 CARACTERES DU NOM DU TYPE D'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JUILLET 1989
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*4       KNOMEL(2)
      include"./incl/nomele.inc"
C
C     LE NUMERO EST IL DANS L'INTERVALLE CORRECT ?
C     ============================================
      IF( NO .LE. 0  .OR.  NO .GT. NBMXEL ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(5)(1:12),'(I12)') NO
         KERR(1) = 'ELNUNM: TYPE EF INCONNU NUMERO '//KERR(5)(1:12)
         CALL LEREUR
         KNOMEL( 1 ) = 'INCO'
         KNOMEL( 2 ) = 'NNU '
         RETURN
      ENDIF
C
C     LE NOM CORRECT EST FOURNI
C     =========================
      KNOMEL( 1 ) = NOMELE( 1, NO )
      KNOMEL( 2 ) = NOMELE( 2, NO )
      RETURN
      END
