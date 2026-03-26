      SUBROUTINE LINVTX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LE FICHIER ./td/i des lignes de texte des invites
C -----
C L'IDEE EST D'ECRASER LES FICHIERS REPERTOIRES td/d ET td/m
C ET LE FICHIER td/i PAR LE FICHIER ET LES
C REPERTOIRES CORRESPONDANTS DANS LA LANGUE DIFFERENTE DU FRANCAIS
C PAR EXEMPLE: POUR L'ANGLAIS LES REPERTOIRES td/da ET td/ma ET td/ia
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1999
C23456---------------------------------------------------------------012
      include"./incl/iinvtx.inc"
      include"./incl/gsmenu.inc"
      include"./incl/homdir.inc"
      include"./incl/langue.inc"
      CHARACTER*160  KNOM
C
C     RECHERCHE D'UN NUMERO D'UNITE
      CALL TRUNIT( NF )
C
C     LE NOM DU FICHIER
      KNOM = HOMDIR // '/td/i'
C     OUVERTURE DU FICHIER SUPPORT DU TABLEAU DES INVITES
      OPEN( FILE=KNOM , UNIT=NF , IOSTAT=I )
      IF( I .NE. 0 ) THEN
          NBLGRC(NRERR) = 1
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'PROBLEME POUR OUVRIR LE FICHIER $MEFISTO/td/i'
          ELSE
             KERR(1) = 'PROBLEM TO OPEN THE FILE $MEFISTO/td/i'
          ENDIF
          CALL LEREUR
          RETURN
      ENDIF
C
C     LECTURE DU TABLEAU TABLEAU KINVTX
      DO 10 I=0,NINVTX
         READ(NF,'(A)',IOSTAT=J) KINVTX(I)
         IF( J .NE. 0 ) THEN
            WRITE(KERR(5)(1:4),'(I4)') I
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR LECTURE LIGNE'//KERR(5)(1:4)
     %                 //' DU FICHIER $MEFISTO/td/i'
               KERR(2) = 'MODIFIER CETTE LIGNE DE CE FICHIER TEXTE'
            ELSE
               KERR(1) = 'INPUT ERROR AT LINE'//KERR(5)(1:4)
     %                 // ' OF FILE $MEFISTO/td/i'
               KERR(2) = 'MODIFY THIS LINE OF THE TEXT FILE'
            ENDIF
            CALL LEREUR
            RETURN
         ENDIF
 10   CONTINUE
C
C     FERMETURE DU FICHIER SUPPORT DU TABLEAU DES INVITES
      CLOSE( UNIT=NF )
      RETURN
      END
