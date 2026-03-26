      SUBROUTINE LXNMNO( NTLX, KNOM, NONOM, MNLX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO NONOM D'UN NOM KNOM DANS LE LEXIQUE NTLX
C -----

C ENTREES:
C --------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C KNOM   : CHAINE DE CARACTERES NOM A RETROUVER DANS LE LEXIQUE

C SORTIE :
C --------
C NONOM  : NUMERO DU NOM KNOM DANS LE LEXIQUE S'IL A ETE RETROUVE
C          =0 S'IL N'A PAS ETE RETROUVE
C MNLX   : ADRESSE MCN DU LEXIQUE NTLX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS     OCTOBRE 1984
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           INOM(18)
      CHARACTER*(*)     KNOM

      NBNOMS = 0
      MNLX   = 0

C     OUVERTURE DU TABLEAU MS = LEXIQUE NTLX
      CALL TAMSOU( NTLX, MNLX )
C     MNLX : ADRESSE MCN DU LEXIQUE NTLX
      IF( MNLX .LE. 0 ) THEN
         NONOM = 0
         GOTO 9999
      ENDIF

C     NOMBRE D ENTIERS POUR UN NOM ET SES ATTRIBUTS DANS LE LEXIQUE
      M1LX = MCN( MNLX )

C     LE NOMBRE MAXIMAL DE NOMS DU LEXIQUE
      MXNOMS = MCN( MNLX + 1 )

C     NBENNM : NOMBRE D'ENTIERS POUR STOCKER LES CARACTERES D'UN NOM
      NBENNM = MCN( MNLX + 2 )

      IF( NBENNM .LE. 0 .OR. NBENNM .GT. 16 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            print*,'lxnmno: NTLX=',NTLX,' Recherche de KNOM=',KNOM
            print*,'lxnmno: NBENNM=',NBENNM,
     %             ' A ETE ECRASE. A REDEMARRER'
         ELSE
            print*,'lxnmno: NTLX=',NTLX,' SEARCH of KNOM=',KNOM
            print*,'lxnmno: NBENNM=',NBENNM,
     %             ' CRUSHED. RESTART'
         ENDIF
ccc         CALL LXIM( NTLX )  2/10/2020 TROP D'AFFICHAGE...
C        LEXIQUE INCORRECT. NOM KNOM NON RETROUVABLE
         GOTO 50 
      ENDIF

C     CONVERSION DE KNOM EN ENTIERS DE INOM
      CALL NOMENT( NBENNM, KNOM, INOM )

C     DEBUT DU CHAINAGE DES NOMS OCCUPES
      NONOM = MCN( MNLX + 5 )

C     LA BOUCLE SUR LES NOMS OCCUPES DU LEXIQUE
C     =========================================
 10   IF( NONOM .GT. 0 ) THEN

C        IL EXISTE ENCORE UN NOM
         NBNOMS = NBNOMS + 1
         IF( NBNOMS .GT. MXNOMS ) THEN
C           POUR EVITER UNE BOUCLE INFINIE
            IF( LANGAG .EQ. 0 ) THEN
               print*,'lxnmno: MXNOMS=',MXNOMS,
     %                ' TROP PETIT pour NBNOMS=',NBNOMS,' A CORRIGER'
            ELSE
               print*,'lxnmno: MXNOMS=',MXNOMS,
     %                ' TOO SMALL for NBNOMS=',NBNOMS,' TO BE CORRECTED'
            ENDIF
ccc            CALL LXIM( NTLX )
            GOTO 50
         ENDIF

C        ADRESSE DU NOM DANS LE TABLEAU DES NOMS DU LEXIQUE
         MN = MNLX + M1LX * NONOM

C        LE NOM OCCUPE EST IL KNOM?
         DO 20 N=1,NBENNM
            IF( INOM(N) .NE. MCN( MN + N - 1) ) GOTO 30
 20      CONTINUE

C        C'EST LE NOM
         GOTO 9999

C        CE N'EST PAS LE NOM KNOM
 30      NONOM = MCN( MN + NBENNM )
         GOTO 10

      ENDIF

C     IL N'EXISTE PLUS DE NOMS. NONOM=0 EN RETOUR
 50   NONOM = 0

 9999 RETURN
      END
