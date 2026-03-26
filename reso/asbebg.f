      SUBROUTINE ASBEBG ( NTDL, NDSM, NBDL, NODL, BE, BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSEMBLER LE SECOND MEMBRE ELEMENTAIRE DANS LE VECTEUR GLOBAL
C -----
C
C ENTREES:
C --------
C NDSM   : NOMBRE DE CAS DE CHARGE
C NBDL   : NOMBRE DE DEGRES DE LIBERTE DE L ELEMENT
C NODL   : NO DES NBDL DEGRES DE LIBERTE DE L ELEMENT
C BE     : SECOND MEMBRE ELEMENTAIRE  (NDSM,NBDL)
C
C SORTIES:
C --------
C BG     : VECTEUR SECOND MEMBRE GLOBAL (NDSM,NTDL)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS FEVRIER 1998
C ......................................................................
      INTEGER           NODL(NBDL)
      DOUBLE PRECISION  BE(NDSM,NBDL) ,BG(NTDL,NDSM)
C
      DO 20 I=1,NBDL
C
C        LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
         IEG = NODL( I )
C
C        ASSEMBLAGE DE LA I-EME COMPOSANTE DE BE DANS BG
         DO 10 J=1,NDSM
            BG( IEG , J ) = BG( IEG , J ) + BE( J , I )
 10      CONTINUE
C
 20   CONTINUE
C
      RETURN
      END
