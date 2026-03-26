      SUBROUTINE FACESP2P1( KNOMOB, NBNOEU, NOSO, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRANSFORMER LE NO (1 a NBNOEU) DE NOEUD DES SOMMETS DES FACES
C -----     D'UN OBJET EN LE NUMERO DES SOMMETS (1 a NBSOMT )
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NBNOEU : NOMBRE DE NOEUDS  DU MAILLAGE DE L'OBJET
C NOSO   : NUMERO DE SOMMET DE CHAQUE NOEUD, 0 SI CE N'EST PAS UN SOMMET
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ST PIERRE DU PERRAY SEPTEMBRE  2013
C23456---------------------------------------------------------------012     
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      include"./incl/a___face.inc"

      CHARACTER*(*)     KNOMOB
      INTEGER           NOSO(NBNOEU)

      IERR = 0
C
C     CREATION OU REDECOUVERTE DU TMS OBJET>>>FACE
      CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
      IF( IERR .NE. 0 ) GOTO 9900
C
C     CREATION OU REDECOUVERTE DU HACHAGE DES ARETES
C     DES FACES FRONTALIERES DE L'OBJET
      CALL HACHAF( KNOMOB, 0, NTFAOB, MNFAOB,
     %             NTAFOB, MNAFOB, IERR )
      IF( IERR .NE. 0 ) GOTO 9900
C
C     LE NOMBRE D'ENTIERS PAR FACE
      MOFACE = MCN( MNFAOB + WOFACE )
C     LA MAJORATION DU NOMBRE DE FACES
      MXFACE = MCN( MNFAOB + WXFACE )
C
C     TRANSFORMATION DU NO DE NOEUD EN NO DE SOMMET
      MNF = MNFAOB + WFACES
      DO NF = 1, MXFACE
C
         IF( MCN( MNF ) .GT. 0 ) THEN
C
C           NOMBRE DE SOMMETS DE LA FACE
            IF( MCN( MNF + 3 ) .EQ. 0 ) THEN
               NBS = 3
            ELSE
               NBS = 4
            ENDIF
C
            DO N = 0, NBS-1
               NOEUD = MCN( MNF + N )
C              LE NUMERO DE SOMMET REMPLACE CELUI DE NOEUD
               MCN( MNF + N ) = NOSO( NOEUD )
            ENDDO
C
         ENDIF
C
         MNF = MNF + MOFACE
C
      ENDDO
C
 9900 RETURN
      END
