      SUBROUTINE AS1BEBG ( NBDL, NODL, BE, BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSEMBLER UN VECTEUR ELEMENTAIRE DANS UN VECTEUR GLOBAL
C -----
C
C ENTREES:
C --------
C NBDL   : NOMBRE DE DEGRES DE LIBERTE DU VECTEUR ELEMENTAIRE (de l'EF)
C NODL   : NO DES NBDL DEGRES DE LIBERTE DE L ELEMENT FINI
C BE     : VECTEUR ELEMENTAIRE (NBDL)
C
C SORTIES:
C --------
C BG     : VECTEUR SECOND MEMBRE GLOBAL (NTDL)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY      Mai 2010
C23456---------------------------------------------------------------012
      INTEGER           NODL(NBDL)
      DOUBLE PRECISION  BE(NBDL), BG(*)
C
      DO 20 I=1,NBDL
C
C        LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DANS L'EF
         IGLOB = NODL( I )
C
C        ASSEMBLAGE DE LA I-EME COMPOSANTE DE BE DANS BG
         BG( IGLOB ) = BG( IGLOB ) + BE( I )
C
 20   CONTINUE
C
      RETURN
      END
