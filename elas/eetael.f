      SUBROUTINE EETAEL( NO,   NDSM,
     %                   NBDL, NCODEM, NTPOBA, NOMTAB, MOTAUX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DONNER LES CARACTERISTIQUES DES TABLEAUX ELEMENTAIRES DES
C -----    ELEMENTS FINIS D ELASTICITE
C
C ENTREES:
C --------
C NO     : NO CODE DE L ELEMENT FINI ( REMIS A JOUR SI INCORRECT )
C NDSM   : NOMBRE DE CAS DE CHARGE
C
C SORTIES:
C --------
C NBDL   : NOMBRE DE DEGRES DE LIBERTE DE L ELEMENT AVANT REDUCTION
C NCODEM : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE DE MASSE
C          -1:NON SYMETRIQUE, 0:DIAGONALE,+1:SYMETRIQUE PLEINE
C NTPOBA : NOMBRE DE TABLEAUX STOCKES SUR LE FICHIER POBA POUR CET EF
C NOMTAB : NOM DES NTPOBA TABLEAUX DE POBA ( 4 CARACTERES CHACUN )
C MOTAUX : NOMBRE DE REELS DOUBLE PRECISION AUXILIAIRES NECESSAIRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1991
C.......................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*4       NOMELE(2)
      INTEGER           NOMTAB(5)
C
C     A PRIORI MATRICE ELEMENTAIRE DE MASSE SYMETRIQUE PLEINE
      NCODEM = 1
      print *,'sp EETAEL: NDSM=',NDSM
C
C     LE NO CODE DE L ELEMENT NO
C     --------------------------
      GOTO( 100,  200,  300, 400,   1,   1,  1,   1,    1,    1,
     %        1,    1,  100,   1, 200, 300,  1, 400, 1900, 2000,
     %     2100, 2200, 2300,2400,   1,   1,  1,2800,  100,    1,
     %        1,    1, 3300,   1,   1) , NO
C
C     LE NOMBRE DE REEL2 AUXILIAIRES S ELEVE A
C     NPI * (4+NDIM*NBPOLY) +
C     MAX ( NBPOLY * ( (NBPOLY+1)/2 + NDIM*NBPOLY ),
C           NDIM * NBPOLY *(NBPOLY+1) + NBPOFA*(NBPOFA+1)/2)
C
C     1 REEL2 = 2 MOTS DE 32 BITS    A MULTIPLIER PAR 2
C     ============================================================
C
C     ERREUR
C     ------
 1    NBLGRC(NRERR) = 2
      CALL ELNUNM( NO , NOMELE )
      KERR(1) ='EETAEL: ELEMENT FINI INCONNU'
      KERR(2) = NOMELE(1) // ' ' // NOMELE(2)
      CALL LEREUR
      CALL ARRET(100)
C
C     TRIA AP1D ET TRIA 2P1D ET TRIA 2P1C
C     ===================================
C     MATRICE ELEMENTAIRE DE MASSE DIAGONALE POUR 2P1D
 100  IF( NO .EQ. 13 .OR. NO .EQ. 29 ) NCODEM = 0
      NBDL   = 6
      NTPOBA = 2
      NOMTAB(1) = ICHARX( '1P13'  )
      NOMTAB(2) = ICHARX( '2P12'  )
      MOTAUX = 72
      GOTO 9999
C
C     TRIA AP2C
C     =========
 200  NBDL   = 12
      NTPOBA = 2
      NOMTAB(1) = ICHARX( '1P25'  )
      NOMTAB(2) = ICHARX( '2P25'  )
      MOTAUX = 202
      GOTO 9999
C
C     QUAD AQ1C ET 2Q1C
C     =================
C     MATRICE ELEMENTAIRE DE MASSE DIAGONALE POUR 2Q1C
 300  NBDL   = 8
      NTPOBA = 2
      NOMTAB(1) = ICHARX( '1P13'  )
      NOMTAB(2) = ICHARX( '2Q13'  )
      MOTAUX = 91
      GOTO 9999
C
C     QUAD AQ2C
C     =========
 400  NBDL   = 16
      NTPOBA =  2
      NOMTAB(1) = ICHARX( '1P25'  )
      NOMTAB(2) = ICHARX( '2Q25'  )
      MOTAUX = 330
      GOTO 9999
C
C     TETR 3P1D
C     =========
C     MATRICE ELEMENTAIRE DE MASSE DIAGONALE POUR 3P1D
1900  NCODEM =   0
      NBDL   =  12
      MOTAUX = 130
      NTPOBA =   2
      NOMTAB(1) =  ICHARX('2P12')
      NOMTAB(2) =  ICHARX('3P11')
      GOTO 9999
C
C     TETR 3P2C
C     =========
2000  NBDL   = 30
      NTPOBA =  2
      NOMTAB(1) =  ICHARX('2P25')
      NOMTAB(2) =  ICHARX('3P25')
C     MOTAUX = NPI * ( 4 + NDIM*NBPOLY ) +
C              MAX ( NBPOLY * ( (NBPOLY+1)/2 + NDIM*NBPOLY ),
C                    NDIM * NBPOLY *(NBPOLY+1) + NBPOFA*(NBPOFA+1)/2 )
C     MOTAUX = 15 * ( 4 + 3*10 ) + MAX(  55 + 300  ,  30 * 11 + 21 )
      MOTAUX = 861
      GOTO 9999
C
C     PENT 3R1C
C     =========
2100  NBDL   = 18
      NTPOBA =  3
      NOMTAB(1) = ICHARX('2P12')
      NOMTAB(2) = ICHARX('3R12')
      NOMTAB(3) = ICHARX('2Q13')
C     MOTAUX = NPI * ( 4 + NDIM*NBPOLY ) +
C              MAX ( NBPOLY * ( (NBPOLY+1)/2 + NDIM*NBPOLY ),
C                    NDIM * NBPOLY *(NBPOLY+1) + NBPOFA*(NBPOFA+1)/2)
C     MOTAUX = 6 * ( 4 + 3*6 ) + MAX(  21 + 108,  3 * 6 * 7 + 10 )
      MOTAUX = 268
      GOTO 9999
C
C     PENT 3R2C
C     =========
2200  NBDL   = 45
      NTPOBA =  3
      NOMTAB(1) = ICHARX('2P25')
      NOMTAB(2) = ICHARX('3R25')
      NOMTAB(3) = ICHARX('2Q25')
C     MOTAUX = NPI * ( 4 + NDIM*NBPOLY ) +
C              MAX ( NBPOLY * ( (NBPOLY+1)/2 + NDIM*NBPOLY ),
C                    NDIM * NBPOLY *(NBPOLY+1) + NBPOFA*(NBPOFA+1)/2)
C     MOTAUX = 21 * ( 4 + 3*15 ) + MAX( 120 + 675, 3 * 15 * 16 + 45 )
      MOTAUX = 1785
      GOTO 9999
C
C     HEXA 3Q1C
C     =========
2300  NBDL   = 24
      NTPOBA =  2
      NOMTAB(1) = ICHARX('2Q13')
      NOMTAB(2) = ICHARX('3Q13')
C     MOTAUX = NPI * ( 4 + NDIM*NBPOLY ) +
C              MAX ( NBPOLY * ( (NBPOLY+1)/2 + NDIM*NBPOLY ),
C                    NDIM * NBPOLY *(NBPOLY+1) + NBPOFA*(NBPOFA+1)/2)
C     MOTAUX = 8 * ( 4 + 3*8 ) + MAX(  36 + 172,  216 + 10 )
      MOTAUX = 450
      GOTO 9999
C
C     HEXA 3Q2C
C     =========
2400  NBDL   = 60
      NTPOBA =  2
      NOMTAB(1) = ICHARX('2Q25')
      NOMTAB(2) = ICHARX('3Q25')
C     MOTAUX = NPI * ( 4 + NDIM*NBPOLY ) +
C              MAX ( NBPOLY * ( (NBPOLY+1)/2 + NDIM*NBPOLY ),
C                    NDIM * NBPOLY *(NBPOLY+1) + NBPOFA*(NBPOFA+1)/2)
C     MOTAUX = 27 * ( 4 + 3*20 ) + MAX( 210 + 1200 , 1260 + 45 )
      MOTAUX = 3024
      GOTO 9999
C
C     SEGM 1P1D
C     =========
 2800 NBDL   = 2
      NTPOBA = 1
      NOMTAB(1) = ICHARX( '1P13'  )
      MOTAUX = 16
      GOTO 9999
C
C     SEGM 1P2D
C     =========
 3300 NBDL   = 3
      NTPOBA = 1
      NOMTAB(1) = ICHARX( '1P25'  )
      MOTAUX = 32
      GOTO 9999
C
 9999 RETURN
      END
