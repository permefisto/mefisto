      SUBROUTINE TRAJFA( NS1, NS2, NS3,
     %                   NF1, NF2, NF3,
     %                   N1TRVI, NOTRIA, NOTRSO, NF )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LA FACE DE SOMMET NS1 2 3 ,
C -----    ET DE FACES ADJACENTES NF1 2 3 DANS NOTRIA
C
C ENTREES:
C --------
C NS1,NS2,NS3 : LES 3 NUMEROS DES SOMMETS DE LA FACE
C NF1,NF2,NF3 : LES 3 NUMEROS DES FACES ADJACENTES PAR LES ARETES
C
C MODIFIES :
C ----------
C NOTRIA : LES 3 SOMMETS, 3 FACES VOISINES ET CHAINAGE
C          DES FACES TRIANGULAIRES DU CONTOUR ET INTERFACES
C NOTRSO : NOTRSO(NS)=NUMERO (DANS NOTRIA) D'UNE FACE DE SOMMET NS
C
C SORTIES:
C --------
C NF     : NUMERO DANS NOTRIA DE LA FACE AJOUTEE
C          0 SI PLUS DE FACE LIBRE DANS NOTRIA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C....................................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NOTRIA(6,1:*),NOTRSO(*),
     %                  NSF(3),NFA(3)
C
C     MISE EN TABLEAUX
      NSF(1) = NS1
      NSF(2) = NS2
      NSF(3) = NS3
      NFA(1) = NF1
      NFA(2) = NF2
      NFA(3) = NF3
C
C     RECHERCHE D'UNE FACE LIBRE DANS NOTRIA
      IF( N1TRVI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SATURATION DES TRIANGLES'
         CALL LEREUR
         NF = 0
         RETURN
      ENDIF
C
C     MISE A JOUR DU 1-ER TRIANGLE VIDE
      NF     = N1TRVI
      N1TRVI = NOTRIA(4,N1TRVI)
C
C     LA FACE NF EST VIDE. AJOUT DES DONNEES
      NOTRIA(1,NF) = NSF(1)
      NOTRIA(2,NF) = NSF(2)
      NOTRIA(3,NF) = NSF(3)
C     LE NUMERO DES FACES ADJACENTES
      NOTRIA(4,NF) = NFA(1)
      NOTRIA(5,NF) = NFA(2)
      NOTRIA(6,NF) = NFA(3)
C
C     LE NUMERO DE FACE NF DANS LES FACES ADJACENTES
      DO 80 I=1,3
         NF0 = NFA(I)
         IF( NF0 .LE. 0 ) GOTO 80
         IF( I .NE. 3 ) THEN
            I1 = I + 1
         ELSE
            I1 = 1
         ENDIF
C        RECHERCHE DU NUMERO D'ARETE DE NFA(I) ADJACENTE A NF
         DO 70 J=1,3
            IF( J .NE. 3 ) THEN
               J1 = J + 1
            ELSE
               J1 = 1
            ENDIF
            IF( (NSF(I).EQ.NOTRIA(J1,NF0).AND.NSF(I1).EQ.NOTRIA(J,NF0)))
     %      THEN
C              L'ARETE I DE NF ET J DE NF0 SONT COMMUNES
               NOTRIA(3+I,NF ) = NF0
               NOTRIA(3+J,NF0) = NF
               GOTO 80
            ENDIF
 70      CONTINUE
 80   CONTINUE
C
C     LE CHAINAGE A ETE MIS A JOUR
C     MISE A JOUR DU TABLEAU NO D'UNE FACE AYANT CES SOMMETS
      NOTRSO( NSF(1) ) = NF
      NOTRSO( NSF(2) ) = NF
      NOTRSO( NSF(3) ) = NF
      END
