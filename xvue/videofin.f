      SUBROUTINE VIDEOFIN( NOMFIC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE le FICHIER Nomfic.gif A PARTIR DES FICHIERS
C -----    CONSTRUITS de NOMS NomficNoImag.jpg
C
C ENTREES:
C --------
C NOMFIC : NOM DU FICHIER puis suivi de 4 CHIFFRES
C NOIMAG : NUMERO DE L'IMAGE  <10000
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR  Janvier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/nmproj.inc"
      CHARACTER*(*) NOMFIC
C
C     SI PAS DE DEMANDE DE FICHIERS VIDEO => RETOUR
      IF( LVIDEO .EQ. 0 ) RETURN
C
C     NOM DU FICHIER DE L'IMAGE NOIMAG SELON LE NO DE L'IMAGE
      NBC = NUDCNB( NOMFIC )
C     NOMBRE DE CARACTERES DU NOM DU PROJET
      LL = NUDCNB( NMPROJ )
C
C     GENERATION DU FICHIER NOMFIC.gif A PARTIR DES FICHIERS NOMFIC*.jpg
      IF( LANGAG .EQ. 0 ) THEN
         print *,'convert ' // NMPROJ(1:LL) // '_' // NOMFIC(1:NBC)
     %// '*.jpg -delay 100 ' // NMPROJ(1:LL) // '_' // NOMFIC(1:NBC)//
     %'.gif    EST EXECUTE'
      ELSE
         print *,'convert ' // NMPROJ(1:LL) // '_' // NOMFIC(1:NBC)
     %// '*.jpg -delay 100 ' // NMPROJ(1:LL) // '_' // NOMFIC(1:NBC)//
     %'.gif  is EXECUTING'
      ENDIF
C
      CALL SYSTEM( 'convert ' // NMPROJ(1:LL) // '_' // NOMFIC(1:NBC)
     %//   '*.jpg -delay 100 ' // NMPROJ(1:LL) // '_' // NOMFIC(1:NBC)//
     %'.gif' )
C
C     DESTRUCTION DES FICHIERS NOMFIC.xwd
      CALL SYSTEM( 'rm -Rf ' // NOMFIC(1:NBC) // '*.xwd' )
C
C     VIDEO TRAITEE => NON DEMANDE DE VIDEO ENSUITE
      LVIDEO = 0
C
C     AFFICHAGE DE LA CONSTRUCTION DU FICHIER VIDEO.gif
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='Le FICHIER VIDEO ' // NOMFIC(1:NBC)//'.gif  EST CREE'
      ELSE
         KERR(1)='The VIDEO FILE ' // NOMFIC(1:NBC)//'.gif  is CREATED'
      ENDIF
      CALL LERESU
C
      RETURN
C
cccC     AUTRE POSSIBILITE DE FAIRE UN FILM apres des XVSAUVERPS(NOMFIC)
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,*) 'Il est possible de faire a partir
ccc          des fichiers.eps un film avec la commande'
ccc      ELSE
ccc         WRITE(IMPRIM,*)
ccc     %'A video from files.eps may be realized with the command'
ccc      ENDIF
ccc
ccc      WRITE(IMPRIM,*)
ccc     %'convert -rotate 90 zfxy*.eps -delay 10 -extent 980x550 zfxy.gif'
C
      END
