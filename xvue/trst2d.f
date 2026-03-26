      SUBROUTINE TRST2D( NCOUL, NOST, XYZSOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    en 2D TRACE des CARACTERES '+NOST' DU SOMMET NOST de XYZSOM
C -----

C ENTREES:
C --------
C NCOUL  : NUMERO DE LA COULEUR DE TRACE
C NOST   : NUMERO XYZSOM du SOMMET NOST
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C2345X7..............................................................012
      REAL           XYZSOM(3,*)
      CHARACTER*12   CARSYM

C     CONSTRUCTION DU MOT '+NOST'
      CARSYM = '+ '
      WRITE( CARSYM(2:11),'(I9)') NOST

C     SUPPRESSION DES BLANCS
      CALL SANSBL( CARSYM, NBC )

cccC     COORDONNEES ECRAN DE L'ITEM
ccc      NX = NUPXEX( XYZSOM(1,NOST) )
ccc      NY = NUPXEY( XYZSOM(2,NOST) )

C     TRACE en NCOUL du MOT '+NOST'
      CALL SYMBOLE2D( NCOUL,XYZSOM(1,NOST),XYZSOM(2,NOST),CARSYM(1:NBC))

      RETURN
      END
