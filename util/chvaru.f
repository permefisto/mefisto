      SUBROUTINE CHVARU( NLD , NCD , NLF , NCF ,
     %                   NOVARU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE A PARTIR DU CARACTERE (NLD,NCD) D'UNE VARIABLE
C ----- UTILISATEUR
C       CE PEUT NE PAS ETRE UNE TELLE VARIABLE
C
C ENTREES :
C ---------
C NLD,NCD : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C
C SORTIES :
C ---------
C NLF,NCF : SI NOVARU>0 POSITION DANS KLG DU DERNIER CARACTERE
C                       DE LA VARIABLE UTILISATEUR
C           SINON NLF=0
C NOVARU :  0 PAS DE VARIABLE RETROUVEE
C          >0 LE NUMERO DE LA VARIABLE UTILISATEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      CHARACTER*1        CAR
      CHARACTER*(NBCAVA) NOMVAR
C
C     RECHERCHE DU NOM DE LA VARIABLE EVENTUELLE COMMENCANT EN NCD
C     ------------------------------------------------------------
C     LE PREMIER CARACTERE
      N   = NCD
      CAR = KLG(NLD)(N:N)
C     LE PREMIER CARACTERE DOIT ETRE UNE LETTRE OU '_'
      IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %            CAR .EQ. '_' ) ) GOTO 9000
C     LES CARACTERES SUIVANTS DE L'EVENTUELLE VARIABLE
      DO 5 N = NCD+1 , NCD+NBCAVA
         CAR = KLG(NLD)(N:N)
         IF(.NOT. ( (CAR .GE. 'A' .AND. CAR .LE. 'Z') .OR.
     %              (CAR .GE. '0' .AND. CAR .LE. '9') .OR.
     %               CAR .EQ. '_' ) ) GOTO 8
C        LE CARACTERE EST ILLICITE
C        REMARQUE : LES MINUSCULES N'EXISTENT PAS CAR TOUTE LIGNE
C                   EST MISE EN MAJUSCULES CF SP LIRLIG
 5    CONTINUE
C     LA POSITION DU DERNIER CARACTERE DE L'EVENTUELLE VARIABLE
 8    NCF = N - 1
C     LE NOM DE L'EVENTUELLE VARIABLE
      NOMVAR = KLG(NLD)(NCD:NCF)
C
c     BOUCLE SUR LES NOMS DES VARIABLES ACTUELLEMENT RECENSEES
C     --------------------------------------------------------
      DO 10 N=NBVARU,1,-1
C        LECTURE RETROGRADE POUR RETROUVER LES VARIABLES LOCALES
C        AVANT LES VARIABLES GLOBALES
         IF( NOMVAR .EQ. KVARU(N) ) THEN
C           LA VARIABLE N EST RETROUVEE
            NLF    = NLD
            NOVARU = N
            RETURN
         ENDIF
 10   CONTINUE
C
C     VARIABLE NON RETROUVEE
 9000 NOVARU = 0
      NLF    = 0
      END
