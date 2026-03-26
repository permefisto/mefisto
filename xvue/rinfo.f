      REAL FUNCTION RINFO ( MOTCLE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER UNE INFORMATION DE TYPE REAL
C -----
C
C ENTREE :
C --------
C MOTCLE : NATURE DE L'INFORMATION A RETOURNER, EN CLAIR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C2345X7..............................................................072
      CHARACTER*(*) MOTCLE
C
C     PLUS PETIT REEL POSITIF
      IF (MOTCLE .EQ. 'PETIT') THEN
         RINFO = 1.E-37
C
C     PLUS PETIT REEL NEGATIF
      ELSE IF (MOTCLE .EQ. '-PETIT') THEN
         RINFO = -1.E-37
C
C     PLUS GRAND REEL POSITIF
      ELSE IF (MOTCLE .EQ. 'GRAND') THEN
         RINFO =  1.E37
C
C     PLUS GRAND REEL NEGATIF
      ELSE IF (MOTCLE .EQ. '-GRAND') THEN
         RINFO = -1.E37
C
C     PRECISION D'UN REEL
      ELSE IF (MOTCLE .EQ. 'PRECISION') THEN
         RINFO = 1.E-6
C
C     ERREUR
      ELSE
         WRITE (IINFO('I'),*) 'RINFO : MOT-CLE INCONNU : ', MOTCLE
         RINFO = 0.
      END IF
      END
