      SUBROUTINE ETTAEL( NO,     NDSM,
     &                   NBPOLY, NCODEM, NTPOBA, NOMTAB, MOTAUX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR D APRES LE NO DU TYPE DE L'EF POUR LA THERMIQUE
C -----    LES CARACTERISTIQUES DE SA MATRICE DE CAPACITE, CONDUCTIVITE,
C          S-MEMBRES ET DES REMONTEES,
C          LE NOMBRE ET LE NOM DES TABLEAUX DU FICHIER POBA,
C          LE NOMBRE DE MOTS AUXILIAIRES NECESSAIRES
C
C ENTREES:
C --------
C NO     : NUMERO DU TYPE DE L'EF
C NDSM   : NOMBRE DE SECONDS MEMBRES OU CAS A TRAITER
C
C SORTIES:
C --------
C NBPOLY : NOMBRE DE DEGRES DE LIBERTE DE L ELEMENT FINI ou
C          NOMBRE DE POLYNOMES DE L'INTERPOLATION DE L'ELEMENT FINI
C NCODEM : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE DE CAPACITE
C          -1:NON SYMETRIQUE, 0:DIAGONALE,+1:SYMETRIQUE PLEINE
C NTPOBA : NOMBRE DE TABLEAUX DE L'EF SUR LE FICHIER POBA
C NOMTAB : NOMTAB(I)=NOM DU I-EME TABLEAU DE POBA DE CE TYPE D'EF
C MOTAUX : NOMBRE DE MOTS DES TABLEAUX AUXILIAIRES
C MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY)
C        + NDIM*(NDIM+1)/2 * ( NDIM + NDIM*(NDIM+1)/2 + NBPOLY ) )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     Octobre 1994
C MODIFS : ALAIN PERRONNET Texas A & M University           Juillet 2005
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Mars    2009
C23456---------------------------------------------------------------012
      PARAMETER  (MOREE2=2)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NOMTAB(5)
      CHARACTER*4       NOMELE(2)
C
C     CAPACITE A PRIORI SYMETRIQUE PLEINE
      NCODEM = 1
C
C     ******************************************************************
      GOTO (110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
     &      210, 220, 230, 240, 120, 130, 140, 140, 290, 300,
     &      310, 320, 330, 340, 350, 350,   1, 380, 390, 400,
     &      410, 420, 430,   1 ), NO
C     ******************************************************************
C
C     ERREUR
C     ------
 1    NBLGRC(NRERR) = 1
      CALL ELNUNM( NO , NOMELE )
      KERR(1) ='ETTAEL: ELEMENT FINI INCONNU'
      KERR(2) = NOMELE(1) // ' '// NOMELE(2)
      CALL LEREUR
      CALL ARRET(100)
C
C     ==================================================================
C     TRIA AP1D
C     ==================================================================
 110  NTPOBA    = 3
      NOMTAB(1) = ICHARX('1P13')
      NOMTAB(2) = ICHARX('2P12')
      NOMTAB(3) = ICHARX('3F13')
      NDIM   = 2
      NBPOLY = 3
      NPI    = 3
      GOTO 2000
C
C     ==================================================================
C     TRIA AP2C ET TRIA 2P2C
C     ==================================================================
 120  NTPOBA    = 3
      NOMTAB(1) = ICHARX('1P25')
      NOMTAB(2) = ICHARX('2P25')
      NOMTAB(3) = ICHARX('3F25')
      NDIM   = 2
      NBPOLY = 6
      NPI    = 7
      GOTO 2000
C
C     ==================================================================
C     QUAD AQ1C ET QUAD 2Q1C
C     ==================================================================
 130  NTPOBA    = 3
      NOMTAB(1) = ICHARX('1P13')
      NOMTAB(2) = ICHARX('2Q13')
      NOMTAB(3) = ICHARX('4F13')
      NDIM   = 2
      NBPOLY = 4
      NPI    = 4
      GOTO 2000
C
C     ==================================================================
C     QUAD AQ2C ET QUAD 2Q2C   npi=9 ndim=2 nbpoly=8
C     ==================================================================
 140  NTPOBA    =  3
      NOMTAB(1) = ICHARX('1P25')
      NOMTAB(2) = ICHARX('2Q25')
      NOMTAB(3) = ICHARX('4F25')
      NDIM   = 2
      NBPOLY = 8
      NPI    = 9
      GOTO 2000
C
C     ==================================================================
C     TRIA MT10
C     ==================================================================
 150  NTPOBA = 0
      MOTAUX = 8 * NDSM
      NBPOLY = 3
      RETURN
C
C     ==================================================================
C     TRIA MT21
C     ==================================================================
 160  NTPOBA    =  1
      NOMTAB(1) = ICHARX('MT21')
      MOTAUX    = 18 * NDSM
      NBPOLY    =  6
      RETURN
C
C     ==================================================================
C     QUAD MQ10
C     ==================================================================
 170  NTPOBA    =  1
      NOMTAB(1) = ICHARX('MQ10')
      MOTAUX    =  10 * NDSM
      NBPOLY    =  4
      RETURN
C
C     ==================================================================
C     QUAD MQ21
C     ==================================================================
 180  NTPOBA    =  1
      NOMTAB(1) = ICHARX('MQ21')
      MOTAUX    = 24 * NDSM
      NBPOLY    = 8
      RETURN
C
C     ==================================================================
C     TETR M3T1
C     ==================================================================
 190  NTPOBA    =  1
      NOMTAB(1) = ICHARX('M3T1')
      MOTAUX    = 10 * NDSM
      NBPOLY    = 4
      RETURN
C
C     ==================================================================
C     TETR M3T2
C     ==================================================================
 200  NTPOBA    =  1
      NOMTAB(1) = ICHARX('M3T2')
      MOTAUX    = 32 * NDSM
      NBPOLY    = 12
      RETURN
C
C     ==================================================================
C     HEXA M3H1
C     ==================================================================
 210  NTPOBA =  1
      MOTAUX = 14 * NDSM
      NBPOLY =  6
      RETURN
C
C     ==================================================================
C     HEXA M3H2
C     ==================================================================
 220  NTPOBA    =  1
      NOMTAB(1) = ICHARX('M3H2')
      MOTAUX    = 2244 + 64 * NDSM
      NBPOLY    = 24
      RETURN
C
C     ==================================================================
C     TRIA 2P1D
C     ==================================================================
C     MATRICE DE CAPACITE CALORIFIQUE (CAPACITE) DIAGONALE
 230  NCODEM = 0
      NTPOBA = 0
      MOTAUX = 24
      NBPOLY = 3
      RETURN
C
C     ==================================================================
C     TRIA 2P2D
C     ==================================================================
 240  NTPOBA = 0
      MOTAUX = 0
      NBPOLY = 6
      RETURN
C
C     ==================================================================
C     TETR 3P1D
C     ==================================================================
C     MATRICE DE CAPACITE CALORIFIQUE (CAPACITE) DIAGONALE
 290  NCODEM = 0
      MOTAUX = 48
      NBPOLY = 4
      NTPOBA = 2
      NOMTAB(1) = ICHARX('2P12')
      NOMTAB(2) = ICHARX('3P11')
      RETURN
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
 300  NTPOBA    =  3
      NOMTAB(1) =  ICHARX('2P25')
      NOMTAB(2) =  ICHARX('3P25')
      NOMTAB(3) =  ICHARX('5F25')
      NDIM   = 3
      NBPOLY = 10
      NPI    = 15
      GOTO 3000
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
 310  NTPOBA    =  4
      NOMTAB(1) = ICHARX('2P12')
      NOMTAB(2) = ICHARX('3R12')
      NOMTAB(3) = ICHARX('6F12')
      NOMTAB(4) = ICHARX('2Q13')
      NDIM   = 3
      NBPOLY = 6
      NPI    = 6
      GOTO 3000
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
 320  NTPOBA    =  4
      NOMTAB(1) = ICHARX('2P25')
      NOMTAB(2) = ICHARX('3R25')
      NOMTAB(3) = ICHARX('6F25')
      NOMTAB(4) = ICHARX('2Q25')
      NDIM   =  3
      NBPOLY = 18
      NPI    = 21
      GOTO 3000
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
 330  NTPOBA    =  3
      NOMTAB(1) = ICHARX('2Q13')
      NOMTAB(2) = ICHARX('3Q13')
      NOMTAB(3) = ICHARX('7F13')
      NDIM   =  3
      NBPOLY =  8
      NPI    =  8
      GOTO 3000
C
C     ==================================================================
C     HEXA 3Q2C     ndim=3 nbpoly=20   npi=27 INTEGRATION EXACTE Q5
C     ==================================================================
 340  NTPOBA    =  3
      NOMTAB(1) = ICHARX('2Q25')
      NOMTAB(2) = ICHARX('3Q25')
      NOMTAB(3) = ICHARX('7F25')
      NDIM   =   3
      NBPOLY =  20
      NPI    =  27
      GOTO 3000
C
cccC     ==================================================================
cccC     HEXA 3Q2C     ndim=3 nbpoly=20   npi=64 INTEGRATION EXACTE Q7
cccC     ==================================================================
ccc 340  NTPOBA    =  3
ccc      NOMTAB(1) = ICHARX('2Q27')
ccc      NOMTAB(2) = ICHARX('3Q27')
ccc      NOMTAB(3) = ICHARX('7F27')
ccc      NDIM   =   3
ccc      NBPOLY =  20
ccc      NPI    =  64
ccc      GOTO 3000
cccC
cccC     ==================================================================
cccC     HEXA 3Q2C      ndim=3 nbpoly=20   npi=125 INTEGRATION EXACTE Q9
cccC     ==================================================================
ccc 340  NTPOBA    =  3
ccc      NOMTAB(1) = ICHARX('2Q29')
ccc      NOMTAB(2) = ICHARX('3Q29')
ccc      NOMTAB(3) = ICHARX('7F29')
ccc      NDIM   =   3
ccc      NBPOLY =  20
ccc      NPI    = 125
ccc      GOTO 3000
C
C     ==================================================================
C     TRIA HD06 ET TRIA EQ06
C     ==================================================================
 350  NTPOBA = 1
      IF( NO-24 .EQ. 1 ) THEN
         NOMTAB(1) = ICHARX('HD06')
         MOTAUX = 72+24*NDSM
      ELSE
         NOMTAB(1) = ICHARX('EQ06')
         MOTAUX = 60+22*NDSM
      ENDIF
      NBPOLY = 6
      RETURN
C
C     ==================================================================
C     TRIA 2P1C   MEME QUE TRIA 2P1C AVEC 2 PT INTEGRATION PAR ARETE
C     ==================================================================
C     MATRICE DE CAPACITE DIAGONALE
 390  NCODEM    = 0
      NTPOBA    = 3
      NOMTAB(1) = ICHARX('1P13')
      NOMTAB(2) = ICHARX('2P12')
      NOMTAB(3) = ICHARX('3F13')
      NDIM   = 2
      NBPOLY = 3
      NPI    = 3
      GOTO 2000
C
C     ==================================================================
C     6CUB 6Q1C      npi=64 ndim=6 nbpoly=64
C     ==================================================================
 400  NTPOBA = 1
C     FORMULE D'INTEGRATION NUMERIQUE EN XYZ A 8 POINTS=SOMMETS  EXACTE Q1
C     FORMULE D'INTEGRATION NUMERIQUE EN UVW A 8 POINTS DE GAUSS EXACTE Q3
C     EN FAIT VU LE PRODUIT TENSORIEL => 8*8=64 POINTS D'INTEGRATION
C     2 FORMULES DIFFERENTES POUR EVITER UNE DISTANCE NULLE ENTRE ELECTRONS
      NOMTAB(1) = ICHARX('6Q13')
      NDIM   = 6
      NBPOLY = 64
      NPI    = 64
      GOTO 3000
C
C     ==================================================================
C     PYRA 3PY1
C     ==================================================================
 410  NTPOBA    =  4
      NOMTAB(1) = ICHARX('2P12')
      NOMTAB(2) = ICHARX('3PY1')
      NOMTAB(3) = ICHARX('FPY1')
      NOMTAB(4) = ICHARX('2Q13')
      NDIM   = 3
      NBPOLY = 5
      NPI    = 8
      GOTO 3000
C
C     ==================================================================
C     PYRA 3PY2
C     ==================================================================
 420  NTPOBA    =  4
      NOMTAB(1) = ICHARX('2P25')
      NOMTAB(2) = ICHARX('3PY2')
      NOMTAB(3) = ICHARX('FPY2')
      NOMTAB(4) = ICHARX('2Q25')
      NDIM   =  3
      NBPOLY = 13
      NPI    = 27
      GOTO 3000
C
C     ==================================================================
C     SEGM 1P1D
C     ==================================================================
 380  NTPOBA    = 2
      NOMTAB(1) = ICHARX('1P13')
      NOMTAB(2) = ICHARX('2F11')
      NDIM   = 1
      NBPOLY = 2
      NPI    = 2
      GOTO 3000
C
C     ==================================================================
C     SEGM 1P2D
C     ==================================================================
 430  NTPOBA    = 2
      NOMTAB(1) = ICHARX('1P25')
      NOMTAB(2) = ICHARX('2F21')
      NDIM   = 1
      NBPOLY = 3
      NPI    = 3
      GOTO 3000
C
C     ==================================================================
C     EN DIMENSION = 2 LE NOMBRE DE MOTS AUXILIAIRES S ELEVE A
C     MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY)
C            + NDIM*(NDIM+1)/2 * ( NDIM + NDIM*(NDIM+1)/2 + NBPOLY ) )
C     La SECONDE LIGNE EST NECESSAIRE AU CALCUL DE L'ESTIMATEUR D'ERREUR
C     POUR LES ELEMENTS FINIS 2D ET POUR LE PROBLEME STATIONNAIRE
C     ==================================================================
 2000 MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY)
     %       + NDIM*(NDIM+1)/2 * ( NDIM + NDIM*(NDIM+1)/2 + NBPOLY ) )
      RETURN
C
C     ==================================================================
C     EN DIMENSION 1 OU >=3 LE NOMBRE DE MOTS AUXILIAIRES S ELEVE A
C     MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY)
C     ==================================================================
 3000 MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY) )
C
C     A SUPPRIMER SI PAS DE CHANGEMENT DANS LES RESULTATS
      MOTAUX = 2* MOTAUX
      RETURN
      END
