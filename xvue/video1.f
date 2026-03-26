      SUBROUTINE VIDEO1( NOMFIC, NOIMAG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  METTRE LA PIXMAP DE LA FENETRE ACTUELLE SUR UN FICHIER DE NOM
C -----  NomficNoImag.jpg avec NoImag forme de 4 CHIFFRES
C
C ENTREES:
C --------
C NOMFIC : NOM GENERIQUE DU FICHIER (Ensuite suivi de 4 chiffres et .jpg)
C NOIMAG : NUMERO<10000 DE L'IMAGE (4 chiffres seulement)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Janvier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/nmproj.inc"
      CHARACTER*(*)  NOMFIC
      CHARACTER*48   NOMFIMAG
      CHARACTER*4    NMIMAG
C
C     SI PAS DE DEMANDE DE FICHIERS VIDEO => RETOUR
      IF( LVIDEO .EQ. 0 ) RETURN
C
C     NOM DU FICHIER DE L'IMAGE NOIMAG SELON LE NO DE L'IMAGE
      NBC = NUDCNB( NOMFIC )
      IF( NOIMAG .LE. 9 ) THEN
         WRITE(NMIMAG(4:4),'(I1)') NOIMAG
         NMIMAG(1:3) = '000'
      ELSE IF( NOIMAG .LE. 99 ) THEN
         WRITE(NMIMAG(3:4),'(I2)') NOIMAG
         NMIMAG(1:2) = '00'
      ELSE IF( NOIMAG .LE. 999 ) THEN
         WRITE(NMIMAG(2:4),'(I3)') NOIMAG
         NMIMAG(1:1) = '0'
      ELSE IF( NOIMAG .LE. 9999 ) THEN
         WRITE(NMIMAG(1:4),'(I4)') NOIMAG
      ENDIF
C
C     NOM DU FICHIER IMAGE.xwd
      NOMFIMAG = NOMFIC(1:NBC) // NMIMAG(1:4) // '.xwd'
      NBC = NBC + 4
C
C     NOMBRE DE CARACTERES DU NOM DU PROJET
      LL = NUDCNB( NMPROJ )
C
C     A PARTIR DE LA PIXMAP de la FENETRE X11 ACTUELLE CREATION DU FICHIER.xwd
      print *
      print *, 'xwd -xy -name Mefisto -out '   //NOMFIMAG(1:NBC+4),
     %         ' EST EXECUTE'
      CALL SYSTEM('xwd -xy -name Mefisto -out '//NOMFIMAG(1:NBC+4) )
C
C     CONVERSION DU FICHIER.xwd en le FICHIER.jpg
      IF( LANGAG .EQ. 0 ) THEN
         print *,    'convert ' // NOMFIMAG(1:NBC+4) // ' ' //
     %  NMPROJ(1:LL)  // '_' // NOMFIMAG(1:NBC)   // '.jpg  EST EXECUTE'
      ELSE
         print *,    'convert ' // NOMFIMAG(1:NBC+4) // ' ' //
     %  NMPROJ(1:LL)  // '_' // NOMFIMAG(1:NBC)   // '.jpg is EXECUTING'
      ENDIF
C
      CALL SYSTEM('convert ' // NOMFIMAG(1:NBC+4) // ' ' //
     %             NMPROJ(1:LL)  // '_' // NOMFIMAG(1:NBC)   // '.jpg' )
C
      RETURN
      END
