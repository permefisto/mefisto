      SUBROUTINE TRST3D( NCOUL, NOST, XYZSOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    en 3D TRACE des CARACTERES '+NOST' DU SOMMET NOST de XYZSOM
C -----

C ENTREES:
C --------
C NCOUL  : NUMERO DE LA COULEUR DE TRACE
C NOST   : NUMERO XYZSOM du SOMMET NOST
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C2345X7..............................................................012
      REAL           XYZSOM(3,*), XYZA(3)
      CHARACTER*12   CARSYM

C     SI COULEUR NEGATIVE PAS DE TRACE
      IF( NCOUL .LT. 0 ) RETURN

C     CONSTRUCTION DU MOT '+NOST'
      CARSYM = '+ '
      WRITE( CARSYM(2:11),'(I9)') NOST

C     SUPPRESSION DES BLANCS
      CALL SANSBL( CARSYM, NBC )

C     COORDONNEES AXONOMETRIQUES DU SOMMET NOST
      CALL XYZAXO( XYZSOM(1,NOST), XYZA )

C     TRACE en NCOUL du MOT '+Nost'
      CALL SYMBOLE2D( NCOUL, XYZA(1), XYZA(2), CARSYM(1:NBC) )

      RETURN
      END
