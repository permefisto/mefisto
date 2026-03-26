      SUBROUTINE ESTASF( IA, L, NELC, NFID, MOPAGE, NEFI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT : ENTREE SORTIE DU TABLEAU D ADRESSE IA DE M DE L MOTS DANS LES
C  ----- ENREGISTREMENTS NEFI+1,NEFI+2,... DU FICHIER NFID EN ACCES
C        DIRECT,CHAQUE PAGE A MOPAGE MOTS
C
C  PARAMETRES D ENTREE :
C  --------------------
C  IA     : ADRESSE DANS M DU TABLEAU A TRANSFERER
C  L      : NOMBRE DE SES MOTS
C  NELC   : NO OPTION : LECTURE=+1,ECRITURE=-1
C  NFID   : NO DU FICHIER EN ACCES DIRECT
C  MOPAGE : NOMBRE DE MOTS D UNE PAGE
C
C  PARAMETRES D ENTREE ET DE SORTIE :
C  ----------------------------------
C  NEFI   : NO-1 DU 1-ER ENREGISTREMENT A ECRIRE OU LIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN  PERRONNET ANALYSE NUMERIQUE UPMC PARIS  MAI 1989
C.......................................................................
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
C
      NFI=NFID
      IF(NFI.LT.0) NFI=-NFID
C
C     LE NOMBRE DE PAGES COMPLETES A TRANSFERER
      NBPAG  = L / MOPAGE
C     LE NOMBRE DE MOTS RESTANTS
      MDPAGE = L - NBPAG * MOPAGE
      IA1    = IA
C
      IF( NBPAG .GT. 0 ) THEN
C
C        TRANSFERT DES PAGES COMPLETES
C        =============================
         DO 4 I=1,NBPAG
            NEFI=NEFI+1
            IF( NELC .GE. 0 ) THEN
               CALL FIPALE( NFI , NEFI , MOPAGE , MCN(IA1) )
            ELSE
               CALL FIPAEC( NFI , NEFI , MOPAGE , MCN(IA1) )
            ENDIF
            IA1=IA1+MOPAGE
    4    CONTINUE
      ENDIF
C
C     TRANSFERT DES MDPAGE MOTS DE LA DERNIERE PAGE
C     =============================================
      IF( MDPAGE .GT. 0 ) THEN
          NEFI=NEFI+1
          IF( NELC .GE. 0 ) THEN
              CALL FIMOLE( NFI , NEFI , MOPAGE , MDPAGE , MCN(IA1) )
          ELSE
              CALL FIMOEC( NFI , NEFI , MOPAGE , MDPAGE , MCN(IA1) )
          ENDIF
      ENDIF

      RETURN
      END
