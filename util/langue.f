      SUBROUTINE LANGUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LE FICHIER $MEFISTO/td/d/debut ET EN DEDUIRE LA LANGUE
C ----- FRANCAIS OU ANGLAIS ET AFFECTER LA VARIABLE LANGAG
C       DU FICHIER INCLUDE $MEFISTO/incl/langue.inc
C
C L'IDEE EST D'ECRASER LES REPERTOIRES td/d ET td/m PAR LES
C REPERTOIRES CORRESPONDANTS DANS LA LANGUE DIFFERENTE DU FRANCAIS
C PAR EXEMPLE: POUR L'ANGLAIS LES REPERTOIRES td/da ET td/ma
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JUILLET 2006
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/homdir.inc"
      include"./incl/langue.inc"
      COMMON /UNITES/LECTEU,IMPRIM,NUNIT(30)
C
      CHARACTER*160  KNOM
      LOGICAL        FEXIST
C
C     LE NOM DU FICHIER
      KNOM = HOMDIR // '/td/m/anglais'
C     LE FICHIER anglais existe t il?
      INQUIRE( FILE=KNOM, EXIST=FEXIST )
      IF( FEXIST ) THEN
C        LANGUE ANGLAISE
         LANGAG = 1
         WRITE(IMPRIM,*) 'Mefisto speaks ENGLISH'
      ELSE
C        LANGUE (INITIALE) FRANCAISE
         LANGAG = 0
         WRITE(IMPRIM,*) 'Mefisto parle FRANCAIS'
      ENDIF
      RETURN
      END
