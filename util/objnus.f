      SUBROUTINE OBJNUS( NMTOBJ , MXNBNO , MNNBNO , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE DES NOMS D' OBJETS ET STOCKER LEUR NUMERO
C ----- DANS LE TABLEAU D'ADRESSE MNNBNO DANS LE SUPER TABLEAU MCN
C       LA LECTURE SE TERMINE PAR LECTURE DU CARACTERE '@'
C
C ENTREE :
C --------
C NMTOBJ : NOM DU TYPE DES OBJETS A LIRE 'POINT' 'LIGNE' ...
C
C ENTREES ET SORTIES :
C --------------------
C MXNBNO : NOMBRE MAXIMAL DE MOTS DU TABLEAU NBNO
C MNNBNO : ADRESSE MCN DU TABLEAU NBNO DANS MCN
C          MCN( MNNBNO ) = NB NOMBRE DE NUMEROS STOCKES
C          MCN( MNNBNO+1 : MNNBNO + NB ) = LES NB NUMEROS
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR , >0 SINON
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS           AVRIL 1987
C......................................................................
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NMTOBJ
      CHARACTER*25      NOMOBJ
C
C     ADRESSE MCN DU LEXIQUE DE CES OBJETS
      N = NTYOBJ( NMTOBJ )
      NTOBJT = NTMN( N )
C
C     NOMBRE DES OBJETS LUS
      NB     = 0
C
C     LECTURE DU NOM DE L'OBJET
 10   GOTO( 21, 22, 23, 24, 25, 26, 27 ),N
 21   CALL INVITE( 51 )
      GOTO 29
 22   CALL INVITE( 40 )
      GOTO 29
 23   CALL INVITE( 42 )
      GOTO 29
 24   CALL INVITE( 60 )
      GOTO 29
 25   CALL INVITE( 45 )
      GOTO 29
 26   CALL INVITE( 38 )
      GOTO 29
 27   CALL INVITE( 44 )
C
 29   NCVALS = 0
      CALL LIRCAR( NCVALS , NOMOBJ )
      IF( NCVALS .EQ. -1 ) THEN
C        NOMOBJ EST '@' ALORS RETOUR
         MCN( MNNBNO ) = NB
         IERR = 0
         RETURN
      ENDIF
C
C     RECHERCHE DU NOM DANS LE LEXIQUE
      CALL LXNMNO( NTOBJT , NOMOBJ , NUMOBJ , MNOBJT )
C
C     L'OBJET A T IL ETE RETROUVE ?
      IF( NUMOBJ .LE. 0 ) THEN
C
C        NON.L'OBJET N'EXISTE PAS
C        -------------------------
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  NMTOBJ//' INCONNU :'//NOMOBJ
         ELSE
            KERR(1) =  NMTOBJ//' UNKNOWN :'//NOMOBJ
         ENDIF
         CALL LEREUR
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19000) NMTOBJ,NOMOBJ
19000       FORMAT(1X,A7,1X,A,' NON RETROUVE PARMI ...')
         ELSE
            WRITE(IMPRIM,29000) NMTOBJ,NOMOBJ
29000       FORMAT(1X,A7,1X,A,' NOT RECOVERED AMONG ...')
         ENDIF
         CALL LXIM0( MNOBJT )
         IF( INTERA .GE. 3 ) GOTO 10
         IERR = 1
         GOTO 10
C
      ELSE
C
C        OUI.LE NUMERO DE L'OBJET EST AJOUTE
C        -----------------------------------
         IF( NB .GE. MXNBNO-1 ) THEN
            CALL TNMCAU( 'ENTIER',MXNBNO,MXNBNO+64,MXNBNO,MNNBNO)
            MXNBNO = MXNBNO + 64
         ENDIF
C        UN NUMERO DE PLUS
         NB = NB + 1
C        EST STOCKE
         MCN( MNNBNO + NB ) = NUMOBJ
         GOTO 10
      ENDIF
      END
