      SUBROUTINE LICHAI( NL , NC , NLF , NCF , NCVALS , CHAINE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE UNE CHAINE DE CARACTERES COMMENCANT PAR ' ET FINISSANT
C ----- AU PROCHAIN ' NON PRECEDE DE '
C
C ENTREES ET SORTIES :
C --------------------
C NL,NC : EN ENTREE POSITION DANS KLG DU PREMIER '
C         EN SORTIE POSITION DANS KLG DU DERNIER ' NON PRECEDE DE '
C
C SORTIES :
C ---------
C NLF,NCF : POSITION DU DERNIER ' DE LA CHAINE
C NCVALS  : 0 SI PAS DE CHAINE  ET NLF=0
C           2 SINON
C CHAINE  : LA CHAINE RETROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*(*)     CHAINE
C
      NLF    = 0
      LCHAI  = LEN( CHAINE )
C     LE NOMBRE DE DOUBLE APOSTROPHES
      NCVALS = 0
      NC1    = 0
      I      = NC
C
 10   I = I + 1
      IF( I .LE. NBCALI ) THEN
C        LE CARACTERE I EST SUR LA LIGNE
         IF( KLG(NL)(I:I) .EQ. '''' ) THEN
C           ' RETROUVE . EST ELLE SUIVIE DE ' ?
            IF( I .LT. NBCALI ) THEN
               IF( KLG(NL)(I+1:I+1) .EQ. '''' ) THEN
C                 NOMBRE DE CARACTERES
                  NBC = I-NC
                  IF( NC1+NBC .GT. LCHAI ) THEN
C                     CHAINE EST TROP PETITE
                      GOTO 9000
                  ENDIF
                  CHAINE(NC1+1:NC1+NBC) = KLG(NL)(NC+1:I)
                  NC1    = NC1 + NBC
                  NC     = I + 1
C                 OUI : PASSAGE AU SUIVANT
                  I      = I + 1
                  GOTO 10
               ENDIF
            ENDIF
         ELSE
C           CE N'EST PAS '
            GOTO 10
         ENDIF
C
C        LE CARACTERE I EST ' ET NON SUIVI DE '
C        NOMBRE DE CARACTERES
         NBC = I-NC-1
         IF( NC1+NBC .GT. LCHAI ) THEN
C           CHAINE EST TROP PETITE
            GOTO 9000
         ENDIF
         CHAINE(NC1+1:NC1+NBC) = KLG(NL)(NC+1:I-1)
         NC  = I
C        LE DERNIER CARACTERE DE LA CHAINE
         NLF = NL
         NCF = I
      ELSE
C        FIN DE LIGNE
         NBC = NBCALI - NC
         IF( NC1+NBC .GT. LCHAI ) THEN
C           CHAINE EST TROP PETITE
            GOTO 9000
         ENDIF
         CHAINE(NC1+1:NC1+NBC) = KLG(NL)(NC+1:NBCALI)
         NC1 = NC1 + NBC
C        PASSAGE A LA LIGNE SUIVANTE
         NL  = NL + 1
         I   = 0
         GOTO 10
      ENDIF
C
C     CHAINE CORRECTE
      NCVALS = 2
      RETURN
C
C     CHAINE INCORRECTE
 9000 NBLGRC(NRERR) = 1
      KERR(1) =  'LU: CHAINE INCORRECTE'
      CALL LEREUR
      NCVALS = 0
      END
