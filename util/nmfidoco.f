      SUBROUTINE NMFIDOCO( NMAEVITER, EXIST, NMFILE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:  RECUPERER L'EVENTUEL NOM DU FICHIER DES DONNEES DE LA COMMANDE
C ----  APPARAISSANT DANS LA LISTE DES ARGUMENTS DE LA COMMANDE
C ENTREE :
C --------
C NMAEVITER: NOM A EVITER DANS NMFILE (NOM DE COMMANDE PRINCIPAL NON ddd)
C
C SORTIES:
C --------
C EXIST  : 1 SI UN NOM DE FICHIER DE DONNEES EST TROUVE
C          0 SI PAS DE NOM
C NMFILE : NOM DU FICHIER DES DONNEES DE LA COMMANDE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & St PIERRE DU PERRAY Decembre 2012
C2345X7..............................................................012
      INTEGER       EXIST
      CHARACTER*(*) NMFILE, NMAEVITER
C
C     NOMBRE DES ARGUMENTS DE LA COMMANDE
      NBARGS = IARGC()
      print *,'Nombre des arguments de la commande Mefisto=',NBARGS
C
      DO 10 N=0,NBARGS
C
C        EXPLORATION DES ARGUMENTS
         CALL GETARG( N, NMFILE )
         print *,'Mefisto: Argument',N,'=',NMFILE
C
C        ELIMINATION DU NOM A EVITER
         IF( INDEX( NMFILE, NMAEVITER ) .GT. 0 ) GOTO 10
         IF( NMFILE(1:3) .EQ. 'ddd' ) GOTO 10
         IF( NMFILE(1:1) .EQ. '<'   ) GOTO 10
         IF( NMFILE(1:1) .EQ. '>'   ) GOTO 10
         IF( NMFILE(1:1) .EQ. '&'   ) GOTO 10
C
C        SEMBLE ETRE LE NOM DU FICHIER
         EXIST = 1
         print *,'Mefisto: Nom du fichier=',NMFILE
         GOTO 9900

 10   CONTINUE
C     AUCUN ARGUMENT SEMBLE ETRE UN NOM DE FICHIER DE DONNEES
      EXIST = 0
C
 9900 RETURN
      END
