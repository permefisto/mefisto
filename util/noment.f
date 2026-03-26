      SUBROUTINE NOMENT( NBENNM , KNOM , INOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER LA CHAINE DE CARACTERES KNOM EN UN TABLEAU D ENTIERS
C ----- INOM DE NBENNM ENTIERS
C
C ENTREES :
C ---------
C NBENNM : NOMBRE D ENTIERS DU TABLEAU INOM
C KNOM   : CHAINE D'AU PLUS 72 CARACTERES
C
C SORTIE  :
C ---------
C INOM   : TABLEAU DE NBENNM ENTIERS CONTENANT KNOM EN SORTIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1984
C.......................................................................
      include"./incl/nbcamo.inc"
      CHARACTER*(*)      KNOM
      CHARACTER*(NBCAMO) BUFFER,BLANC
      INTEGER            INOM(NBENNM)
      INTRINSIC          LEN
      DATA               BLANC/'    '/
C
C     LES CARACTERES DE KNOM SONT TRADUITS EN ENTIER
      L = LEN( KNOM )
C
C     NOMBRE DE MOTS OCCUPES PAR KNOM
      MOTS = ( L - 1 ) / NBCAMO + 1
      IF( MOTS .LE. 0 ) MOTS = 1
C
C     SI LE NOMBRE DE MOTS DE KNOM EXCEDE NBENNM IL EST RAMENE A NBENNM
      IF( MOTS .GT. NBENNM ) MOTS = NBENNM
C
C     CONVERSION DES MOTS DE KNOM
      L1 = 1
      L2 = MIN( NBCAMO , L )
      DO 10 I=1,MOTS
C        BUFFERISATION POUR GENERER DES BLANCS EN FIN DE MOT
         BUFFER    = KNOM(L1:L2)
         INOM( I ) = ICHARX( BUFFER )
         L1        = L1 + NBCAMO
         L2        = MIN( L2 + NBCAMO , L )
 10   CONTINUE
C
C     COMPLETION PAR DES BLANCS
C
C     ATTENTION SI PLUS DE 4 CARACTERES PAR MOTS
C     ========= METTRE A JOUR LA DECLARATION ET LE DATA
C                             DE BLANC CI DESSOUS
      DO 20 I=MOTS+1 , NBENNM
         INOM( I ) = ICHARX( BLANC )
 20   CONTINUE
C
      RETURN
      END
