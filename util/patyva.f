      SUBROUTINE PATYVA( NOFICH , NOIDEN , IMPOSS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE PARAMETER DE L'IDENTIFICATEUR NOIDEN
C -----     DU TABLEAU IDENT
C
C ENTREES :
C ---------
C NOFICH  : NUMERO D'UNITE DU FICHIER DES PARAMETER
C NOIDEN  : NUMERO DE L'IDENTIFICATEUR
C IMPOSS  : 3 SI LE PARAMETER NE DOIT PAS ETRE CREE PAR CAUSE
C             DE BORNES D'INDICES DEFINIES PAR IDENTIFICATEUR
C           0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      include"./incl/msvaau.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*80      KNOM
      CHARACTER*12      KENTIE
C
C     RECHERCHE DU NOMBRE DE VARIABLES SI C'EST UN TABLEAU
      CALL NVARID( NOIDEN , NBV , NOTYPE )
C     LE TYPE DE L'IDENTIFICATEUR EST NOTYPE
      IF(  NOTYPE .LE. 0 .OR.
     %   ( NOTYPE .GT. NBTYPV .AND. NOTYPE .NE. 21 ) ) THEN
          NBLGRC(NRERR) = 1
          WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
          KERR(1) = 'ERREUR PATYVA:'//KERR(MXLGER)(1:4)
     %            //' TYPE INCORRECT'
          CALL LEREUR
C         ON BLOQUE LA GENERATION PAR IMPOSS=10
          IMPOSS = 10
          RETURN
      ENDIF
C
      IF( NOMUET .EQ. 0 .OR. IMPOSS .GT. 2 ) GOTO 100
C     L'IMPOSSIBILITE D'ADRESSER LA SUITE EST RETARDEE
      IF( IMPOSS .EQ. 1 .OR. IMPOSS .EQ. 2 ) IMPOSS = 3
C
C     GENERATION DU NOM DU PARAMETER
      WRITE(KENTIE,'(I12)' ) LDTS( LHTMS )
      I = INDEX( KIDENT(NOIDEN) , ' ' )
      IF( I .EQ. 7 ) THEN
C        LE 1-ER CARACTERE DE L'IDENTIFICATEUR EST ECRASE PAR W
         I1 = 2
         I2 = 6
      ELSE
C        W EST AJOUTE DEVANT LE NOM DE L'IDENTIFICATEUR
         I1 = 1
         I2 = I - 1
      ENDIF
C
C     LE PARAMETER EST PORTE SUR LE FICHIER
      KNOM = '     % , W' // KIDENT(NOIDEN)(I1:I2) // ' = ' // KENTIE
C
C     LA ' EST SUPPRIMEE LORS DU PREMIER PARAMETER
      IF( NOIDEN .EQ. 3 ) KNOM(8:8) = ' '
C
      WRITE( NOFICH , '(A)' ) KNOM
C
C     MISE A JOUR DU TABLEAU LDTS TANT QUE IMPOSS=0
C     LE NOMBRE DE MOTS DE DECALAGE
 100  NBV = NBV * MOTVAR( NOTYPE )
C
C     MISE A JOUR DU DECALAGE DANS LE TMS
      LDTS(LHTMS) = LDTS(LHTMS) + NBV
C
      RETURN
      END
