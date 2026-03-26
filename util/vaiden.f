      SUBROUTINE VAIDEN( NOIDEN , NVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOURNER LA VALEUR ENTIERE DANS LE TABLEAU MS
C -----     DE L'IDENTIFICATEUR NOIDEN DU TABLEAU IDENT
C
C ENTREES :
C ---------
C NOIDEN  : NUMERO DE L'IDENTIFICATEUR
C
C SORTIES :
C ---------
C NVAL    : VALEUR ENTIERE DE L'IDENTIFICATEUR
C           LE PLUS GRAND ENTIER >0 SI L'IDENTIFICATEUR
C           N'EST PAS ENTIER OU ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      include"./incl/pp.inc"
      COMMON MCN(MOTMCN)
C
C     RECHERCHE DES CARACTERISTIQUES DE L'IDENTIFICATEUR
      IF( IDENT(1,NOIDEN) .EQ. 4 .OR. IDENT(1,NOIDEN) .EQ. 11) THEN
         DO 10  I = LHPILE , 1 , -1
            IF( LAPILE(1,I) .EQ. 1 ) THEN
C              LE NUMERO DE L'ARTICLE EN COURS DE TRAITEMENT
               LHA = LAPILE(0,LHPILE)
C              L'ADRESSE MCN DU DEBUT DU TABLEAU MS
               MN = MCTAMS( LHA )
C              L'ADRESSE DANS LE TABLEAU MCN DU 1-ER MOT DE L'ENTIER
               MN = MN + IDENT(3,NOIDEN)
C              LA VALEUR
               NVAL = MCN( MN )
               RETURN
            ENDIF
 10      CONTINUE
      ENDIF
C
C     IDENTIFICATEUR NON RETROUVE
      NVAL = IINFO( 'GRAND' )
      END
