      SUBROUTINE LXNMDS( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE LE NOM KNOM DANS LE LEXIQUE NTLX
C -----    ET LE TABLEAU TMS QUI LUI EST ASSOCIE S'IL EXISTE
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C KNOM   : CHAINE DE CARACTERES NOM A RETROUVER DANS LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS     OCTOBRE 1985
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     KNOM
C
      IF( NTLX .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'LXNMDS: LEXIQUE DE TMS NTLX<=0'
         KERR(2) = KNOM
         CALL LEREUR
         CALL XVPAUSE
      ENDIF
C
C     RECHERCHE DU NOM KNOM PARMI LES NOMS DU LEXIQUE NTLX
      CALL LXNMNU( NTLX , KNOM , NONOM0 , NONOM1 , MNLX )
C
C     LE NOM A T IL ETE RETROUVE ?
      IF( NONOM1 .LE. 0 ) THEN
C
C        NON. NOM NON RETROUVE
C        ===  ================
         NBLGRC(NRERR) = 2
         KERR(1) = 'LXNMDS: TABLEAU ' // KNOM
         KERR(2) = 'INCONNU DONC NON DETRUIT'
         CALL LEREUR
         CALL LXIM( NTLX )
         RETURN
C
      ENDIF
C
C     OUI. NOM RETROUVE
C     ===  ============
C     IL EST DETRUIT DES NOMS OCCUPES
C     POUR DEVENIR LE PREMIER NOM LIBRE
      M1LX   = MCN( MNLX )
      NBENNM = MCN( MNLX + 2 )
C     ADRESSE DU 1-ER MOT DU NOM KNOM DANS LE LEXIQUE
      MN = MNLX + M1LX * NONOM1
C
C     LE NOM EST EFFACE  15/5/2012
      MCN(MN) = 0
C
C     LE NUMERO DU TABLEAU MS ET S'IL EXISTE SA DESTRUCTION
      NOTAMS = MCN( MN + NBENNM + 2 )
      IF( NOTAMS .GT. 0 ) CALL TAMSDS( NOTAMS )
C
C     MISE A JOUR DU CHAINAGE
      IF( NONOM0 .GT. 0 ) THEN
C         IL EXISTE UN NOM PRECEDANT KNOM
C         ADRESSE DU 1-ER MOT DU NOM PRECEDANT KNOM
          MN0 = MNLX + M1LX * NONOM0
          MCN( MN0 + NBENNM ) = MCN( MN + NBENNM )
      ELSE
C         IL N'EXISTE PAS DE NOM PRECEDANT KNOM
          MCN( MNLX + 5 ) = MCN( MN + NBENNM )
      ENDIF
C
C     NONOM1 DEVIENT LE 1-ER NOM LIBRE
      MCN( MN + NBENNM ) = MCN( MNLX + 6 )
      MCN( MNLX + 6 ) = NONOM1
C
      RETURN
      END
