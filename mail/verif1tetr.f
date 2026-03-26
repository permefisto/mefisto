      SUBROUTINE VERIF1TETR( NT, NOTETR, NBSOMM, PTXYZD, V )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    VERIFIER LE NUMERO DES 4 SOMMETS DU TETRAEDRE NT
C -----    CALCULER SON VOLUME
C
C ENTREES:
C --------
C NT     : NUMERO DU TETRAEDRE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NBSOMM : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES SOMMETS
C
C SORTIES:
C --------
C V      : VOLUME SIGNE DU TETRAEDRE NT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    AOUT 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NOTETR(8,*)
      DOUBLE PRECISION  PTXYZD(4,NBSOMM), V, VOLTET
C
C     LES NUMEROS DES 4 SOMMETS DU TETRAEDRE NT DE NOTETR
      DO 20 I=1,4
         NS = NOTETR(I,NT)
         IF( NS .LE. 0 .OR. NS .GT. NBSOMM ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10020) NT,(NOTETR(J,NT),J=1,4)
            ELSE
               WRITE(IMPRIM,20020) NT,(NOTETR(J,NT),J=1,4)
            ENDIF
         ENDIF
 20   CONTINUE
C
C     LE VOLUME DU TETRAEDRE NT EST-IL POSITIF ?
      V = VOLTET( PTXYZD(1,NOTETR(1,NT) ) ,
     %            PTXYZD(1,NOTETR(2,NT) ) ,
     %            PTXYZD(1,NOTETR(3,NT) ) ,
     %            PTXYZD(1,NOTETR(4,NT) ) )
C
      IF( V .LE. 0.D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10020) NT,(NOTETR(J,NT),J=1,4),V
         ELSE
            WRITE(IMPRIM,20020) NT,(NOTETR(J,NT),J=1,4),V
         ENDIF
10020 FORMAT(/' TETRAEDRE ',I8,' NO SOMMETS =',4I8,' DE VOLUME ',
     %          G15.6,' <0 INCORRECT'/
     %        ' EXECUTER APRES AMELIORATION D''UNE TETRAEDRISATION'/)
20020 FORMAT(/' TETRAHEDRON ',I8,' VERTEX NUMBERS =',4I8,' of VOLUME ',
     %          G15.6,' <0 INCORRECT'//
     %        ' EXECUTE AFTER the TETRAHEDRIZATION IMPROVEMENT'/)
C
CCC      SINON PERMUTATION DE 2 SOMMETS
CCC      J         = NOTETR(3)
CCC      NOTETR(3) = NOTETR(4)
CCC      NOTETR(4) = J
      ENDIF
C
      RETURN
      END
