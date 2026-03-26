      CHARACTER*48 FUNCTION FICHNM( NOFIMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DONNER UN NOM AU FICHIER NOFIMS DE LA MEMOIRE SECONDAIRE (MS)
C ----  ATTENTION CE NUMERO DOIT ETRE POSITIF ET INFERIEUR A 89
C***********************************************************************
C ATTENTION: SOUS-PROGRAMME DEPENDANT MACHINE ICI VERSION APOLLO
C***********************************************************************
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  DECEMBRE 1983
C.......................................................................
      include"./incl/gsmenu.inc"
      INTEGER NOFIMS
C     CE NOMBRE DOIT ETRE COMPRIS ENTRE 1 ET 89 POUR LA VERSION IBM
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSNMFI / NMFIC
      CHARACTER*16      NMFIC
      CHARACTER*2       KFICH
C
      IF( NOFIMS .LT. 0 .OR. NOFIMS .GT. 89 ) THEN
C        ERREUR
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:4),'(I4)') NOFIMS
         KERR(1) = 'FICHNM:NUMERO DE FICHIER MS='//KERR(MXLGER)(1:4)
         KERR(2) = 'NON COMPRIS ENTRE 0 (SF) ET 89'
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C
C     LE NOMBRE 10+NOFIMS EST TRADUIT EN CARACTERES
      WRITE(UNIT=KFICH,FMT='(I2)') 10 + NOFIMS
C
      FICHNM = 'MS' // KFICH
C
C     ECRITURE DES LETTRES EN MINUSCULES
      CALL MINUSC( FICHNM )
      END
