      SUBROUTINE AS2BEBG ( NBDL, NODL, BE, NBNOEM, BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSEMBLER UN VECTEUR ELEMENTAIRE DANS UN VECTEUR GLOBAL
C -----    AVEC 2 DEGRES DE LIBERTE PAR NOEUD
C
C ENTREES:
C --------
C NBDL   : NOMBRE DE NOEUDS DE l'EF
C NODL   : NO DES NBDL NOEUDS DE L ELEMENT FINI
C BE     : VECTEUR ELEMENTAIRE (NBDL,2)
C NBNOEM : NOMBRE DE NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C BG     : VECTEUR SECOND MEMBRE GLOBAL  (NBNOEM,2)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      INTEGER           NODL(NBDL)
      DOUBLE PRECISION  BE(NBDL,2), BG(NBNOEM,2)
C
      DO I=1,NBDL
C
C        LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DANS L'EF
         IGLOB = NODL( I )
C
C        ASSEMBLAGE DE LA I-EME COMPOSANTE DE BE DANS BG
         BG( IGLOB, 1  ) = BG( IGLOB, 1 ) + BE( I, 1 )
         BG( IGLOB, 2  ) = BG( IGLOB, 2 ) + BE( I, 1 )
C
      ENDDO
C
      RETURN
      END
