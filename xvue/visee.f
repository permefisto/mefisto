      SUBROUTINE VISEE( NMTCL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LES OPTIONS DE VISEE D'UN OBJET 1D 2D OU 3D
C -----
C
C SORTIE :
C --------
C NMTCL  : NUMERO DU DERNIER MOT CLE ACTIVE DANS LA VISEE
C          <0 DEMANDE UN ABANDON DU TRACE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C ......................................................................
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
C
C     DIMENSION DE L'ESPACE DE TRAVAIL
      IF( COOEXT(3,1) .NE. 0.0 .OR. COOEXT(3,2) .NE. 0.0 ) THEN
C
C        OBJET 3D
         NDIMLI = 3
         CALL VISE3D( NMTCL )
C
      ELSE IF( COOEXT(2,1) .NE. 0.0 .OR. COOEXT(2,2) .NE. 0.0 ) THEN
C
C        OBJET 2D
         NDIMLI = 2
         CALL VISE2D( NMTCL )
C
      ELSE
C
C        OBJET 1D
         NDIMLI = 1
         CALL VISE1D( NMTCL )
C
      ENDIF
C
      RETURN
      END
