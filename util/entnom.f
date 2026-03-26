      SUBROUTINE ENTNOM( NBENNM, INOM, KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER LES NBENNM ENTIERS DE INOM EN UNE CHAINE DE
C ----- CARACTERES KNOM
C
C ENTREES :
C ---------
C NBENNM  : NOMBRE D ENTIERS DU TABLEAU INOM
C INOM    : TABLEAU DE NBENNM ENTIERS
C
C SORTIE  :
C ---------
C KNOM    : CHAINE D'AU PLUS 72 CARACTERES TRADUCTION CARACTERES DE INOM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS     OCTOBRE 1984
C.......................................................................
      include"./incl/langue.inc"
      include"./incl/nbcamo.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*(NBCAMO) CHARX
      CHARACTER*(*)      KNOM
      INTEGER            INOM(NBENNM)
C
C     LE NOMBRE DE CARACTERES DE LA CHAINE KNOM
      L = LEN( KNOM )
C
C     LE NOMBRE D'ENTIERS POUR CODER DE NOM
      MOTS = ( L - 1 ) / NBCAMO + 1
      IF( MOTS .LT. NBENNM ) THEN
C        PROBLEME:KNOM NE PEUT RECEVOIR LE CODAGE DES NBENNM ENTIERS
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ENTNOM: CHAINE KNOM TROP COURTE'
         ELSE
            KERR(1) = 'ENTNOM: STRING KNOM TOO SHORT'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE DECODAGE DES NBENNM ENTIERS EN CARACTERES
      L1 = 1
      L2 = NBCAMO
      DO 10 I=1,NBENNM
         KNOM(L1:L2) = CHARX( INOM(I) )
         L1 = L1 + NBCAMO
         L2 = L2 + NBCAMO
 10   CONTINUE
C
C     COMPLETION DES CARACTERES SUPPLEMENTAIRES DE KNOM PAR DES BLANCS
      DO 20 I=L2-NBCAMO+1,L
         KNOM(I:I) = ' '
 20   CONTINUE
C
      RETURN
      END
