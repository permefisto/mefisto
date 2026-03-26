      SUBROUTINE NOBLNO( NBNOE, NBDLNO, NONOEF, NODLNO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PASSER DES NBNOE DL NUMEROTES DANS NBDLNO BLOCS
C -----    A LA NUMEROTATION DES DL RANGES PAR NOEUDS
C
C ENTREES:
C --------
C NBNOE  : NOMBRE DE NOEUDS DE L'EF
C NBDLNO : NOMBRE DE DL PAR NOEUD
C NONOEF : NUMERO DES NBNOE NOEUDS DE L'EF
C
C SORTIE :
C --------
C NODLNO : NUMERO DES DL RANGES PAR NOEUDS GLOBAUX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      INTEGER  NONOEF(NBNOE), NODLNO(NBNOE,NBDLNO)
C
      DO I = 1, NBNOE
C
C        NUMERO DU DL QUI PRECEDE LE PREMIER DL DU NOEUD NONOEF(I)
         NUMDL = NBDLNO * (NONOEF(I)-1)
C
         DO J = 1, NBDLNO
C
C           J EST LE NUMERO DU BLOC DANS LE VECTEUR ELEMENTAIRE
C
C           LE DL I DU BLOC J EST ENVOYE COMME DL J DU NOEUD NONOEF(I)
            NODLNO(I,J) = NUMDL + J
C
         ENDDO
C
      ENDDO
C
      RETURN
      END
