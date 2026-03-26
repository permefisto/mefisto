      SUBROUTINE GIBBB ( MNNOEU, MNRENU )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    REGENERATION DES COORDONNEES EN ACCORD AVEC LA NOUVELLE
C -----    RENUMEROTATION DES NOEUDS

C ENTREE :
C --------
C MNNOEU : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET
C MNRENU : ADRESSE MCN DU TABLEAU DE RENUMEROTATION DES NOEUDS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: DEFAIX THIERRY   ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1990
C2345X7..............................................................012
      include"./incl/a___xyznoeud.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))

C     INITIALISATION ET OUVERTURE DU TABLEAU INTERMEDIAIRE
C     ====================================================
      NBNOE = MCN( MNNOEU + WNBNOE )
      CALL TNMCDC( 'MOTS', 3*NBNOE, MNINTE )

C     TRANSFERT DES XYZ DANS LE TABLEAU INTERMEDIAIRE
C     ===============================================
      MNI = MNINTE
      MNN = MNNOEU + WYZNOE
      DO I=1,NBNOE
         DO J=0,2
            RMCN( MNI+J ) = RMCN( MNN+J )
         ENDDO
         MNI = MNI + 3
         MNN = MNN + 3
      ENDDO

C     TRANSFERT DU TABLEAU INTERMEDIAIRE EN LES XYZ RENUMEROTES
C     =========================================================
      MNRENU1 = MNRENU - 1
      MNI = MNINTE
      MNN = MNNOEU + WYZNOE
      DO I=1,NBNOE
         NEWI   = MCN( MNRENU1+I ) - 1
         MNNEWI = MNN + NEWI*3
         DO J=0,2
            RMCN( MNNEWI + J ) = RMCN( MNI + J )
         ENDDO
         MNI = MNI + 3
      ENDDO

C     DESTRUCTION DU TABLEAU NTERMEDIAIRE
C     ===================================
      CALL TNMCDS( 'MOTS', 3*NBNOE, MNINTE )

      RETURN
      END
