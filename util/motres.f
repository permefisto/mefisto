      SUBROUTINE MOTRES( MOT, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER 0 SI MOT N'EST PAS UN MOT RESERVE ET 1 SINON
C -----
C CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREE :
C --------
C MOT    : LE MOT A COMPARER A LA LISTE DES MOTS RESERVES
C
C SORTIE :
C --------
C NONOUI : 1 SI LE MOT   EST     UN MOT RESERVE
C          0 SI LE MOT N'EST PAS UN MOT RESERVE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1995
C23456---------------------------------------------------------------012
      PARAMETER          (NCMOLU=8, NBMOLU=2)
      CHARACTER*(*)       MOT
      CHARACTER*(NCMOLU)  MOTDLU(1:NBMOLU)
C
      DATA      MOTDLU / 'AFFICHER', 'DISPLAY ' /
C
C     LE NOMBRE DE CARACTERES DU MOT
      NBCMOT = NUDCNB( MOT )
      DO 10 I=1,NBMOLU
C        LE NOMBRE DE CARACTERES DU MOT RESERVE I
         NBCI = NUDCNB( MOTDLU(I) )
         IF( MOT(1:NBCMOT) .EQ. MOTDLU(I)(1:NBCI) ) THEN
            NONOUI = 1
            RETURN
         ENDIF
 10   CONTINUE
C
C     MOT NON RESERVE
      NONOUI = 0
      RETURN
      END
