      SUBROUTINE TTNOOU( MNTATA,NOTATA,MCNOTA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : OUVRIR LE TABLEAU NUMERIQUE NOTATA DU TABLEAU DE TABLEAUX
C ----
C ENTREES :
C ---------
C MNTATA : ADRESSE MCN DU TABLEAU DES NUMEROS DE TAMS DES TABLEAUX
C NOTATA : NUMERO DU TABLEAU A DECLARER (1 A MXTATA)
C
C SORTIE :
C --------
C MCNOTA : ADRESSE MC DU TABLEAU NOTATA DU TATA APRES OUVERTURE
C          0 SI LE TABLEAU NE PEUT ETRE OUVERT => TEST EN SORTIE
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
         WRITE(KERR(MXLGER)(1:10),'(I10)') NOTATA
         NBLGRC(NRERR) = 1
         KERR(1) = 'NUMERO TABLEAU INCORRECT'//KERR(MXLGER)(1:10)
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU EST IL DEJA DECLARE ?
      MN = MNTATA + 3 + NOTATA
      IF( MCN(MN) .EQ. 0 ) THEN
C        LE TABLEAU N'EST PAS DECLARE => ERREUR
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NOTATA
         KERR(1) = 'LE TABLEAU '//KERR(MXLGER)(1:10) //
     %             'DU TATA N''EST PAS DECLARE DONC NON OUVRABLE'
         CALL LEREUR
         MCNOTA = 0
      ELSE
C        LE TABLEAU EST OUVERT
         NT = ABS( MCN(MN) )
         CALL TAMSOU( NT , MCNOTA )
         MCN( MN ) = -NT
      ENDIF
      END
