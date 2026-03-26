      SUBROUTINE ELCGNM( NCOGEL, NBMOTS, NOMELE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR SELON LE CODE GEOMETRIQUE  LE NOM DE L ELEMENT FINI
C -----    SUR 1 A 3 MOTS DE 4 CARACTERES
C
C ENTREES :
C ---------
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT
C          1:NOEUD 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE 5:TETRAEDRE
C          6:PENTAEDRE 7:HEXAEDRE 8:6-CUBE
C NBMOTS : OPTION DE TRAITEMENT DU SP
C           1               NCOGEL    => NOMELE(1) EN SORTIE
C           2               NCOGEL    => NOMELE(1),(2)
C           3               NCOGEL    => NOMELE(1),(2),(3)
C
C RESULTATS :
C -----------
C NOMELE : 1 OU 2 OU 3 MOTS CONTENANT LES CARACTERES(4 PAR MOT)
C          DU NOM DE L ELEMENT GEOMETRIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS     JANVIER 1991
C ......................................................................
      include"./incl/gsmenu.inc"
      CHARACTER*4       LISTEL(3,8)
      CHARACTER*4       NOMELE(*)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DATA LISTEL / 'NOEU' , 'D   ' , '    ' ,
     &              'SEGM' , 'ENT ' , '    ' ,
     &              'TRIA' , 'NGLE' , '    ' ,
     &              'QUAD' , 'RANG' , 'LE  ' ,
     &              'TETR' , 'AEDR' , 'E   ' ,
     &              'PENT' , 'AEDR' , 'E   ' ,
     &              'HEXA' , 'EDRE' , '    ' ,
     &              '6-CU' , 'BE  ' , '    ' /
C
      IF( NCOGEL .LE. 0  .OR.  NCOGEL .GT. 8 ) THEN
          NBLGRC(NRERR) = 1
          KERR(1) ='ELCGNM:CODE GEOMETRIE INCONNU'
          CALL LEREUR
          NOMELE(1) = 'ERRE'
          RETURN
      ENDIF
C
C     PROTECTION DE NBMOTS
      N = MIN( 3 , NBMOTS )
      N = MAX( 1 , N      )
      DO 220 I=1,N
         NOMELE(I) = LISTEL(I,NCOGEL)
  220 CONTINUE
      RETURN
      END
