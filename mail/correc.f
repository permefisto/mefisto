      SUBROUTINE CORREC( MOT,LONG )
C ......................................................................
C BUT : CORRIGER(DANS UNE CERTAINE MESURE) LA SYNTAXE DES MOTS ENTRES
C       AU CLAVIER ET EXTRAIRE LEURS LONGUEURS
C ......................................................................
C ENTREE : MOT  : CHAINE DE CARACTERES CONTENANT LE MOT
C          LONG : LONGUEUR MAXIMALE DE CETTE CHAINE
C SORTIE : MOT  : LE MOT CORRIGE
C          LONG : LONGUEUR DU MOT
C ......................................................................
C AUTEURS : A.GOLGOLAB ET X.DENG  ENS-CACHAN  DEC 1987
C ......................................................................
C
      CHARACTER*(*) MOT
      INTEGER LONG,I,POSD,IC,OLDLON
      CHARACTER*1 C
C
      POSD = 0
C
C - Convertir en majuscule et la position du dernier caractere -
      OLDLON = LONG
      DO I= 1 , LONG
         C = MOT(I:I)
         IF ( C .NE. ' ' ) POSD = I
         IC = ICHAR(C)
         IF ( IC.GE.97 .AND. IC.LE.122 ) MOT(I:I) = CHAR(IC-32)
      ENDDO
      LONG = POSD
C - remplacement des blancs par un '%' -
      DO I= 1 , LONG
         IF ( MOT(I:I) .EQ. ' ' )  MOT(I:I) = '%'
      ENDDO
C - Elimination des blancs du debut -
      DOWHILE ( MOT(1:1) .EQ. '%' )
         LONG = LONG -1
         MOT(1:LONG) = MOT(2:LONG+1)
      ENDDO
C - Mise a blanc de la fin du mot -
      IF ( LONG+1 .LE. OLDLON ) MOT(LONG+1:OLDLON) = ' '
C - Elimination des blancs intermediares inutiles -
      POSD = 2
      DOWHILE ( POSD + 1 .LT. LONG )
         IF ( MOT(POSD:POSD+1) .EQ. '%%' ) THEN
            MOT(POSD+1:LONG) = MOT(POSD+2:LONG) // ' '
            LONG = LONG - 1
         ELSE
            POSD = POSD + 1
         ENDIF
      ENDDO
C - ELIMINATION DES % -
      DO I = 1 , LONG
         IF ( MOT(I:I) .EQ. '%' ) MOT (I:I) = ' '
      ENDDO
      END
