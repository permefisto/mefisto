      FUNCTION LAPXMENU( NBL, TEXTE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   NOMBRE DE PIXELS DE LA LARGEUR DU TRACE DU TEXTE POUR UN MENU
C -----           VERSION xvue et sans Motif
C
C ENTREE :
C --------
C NBL    : NOMBRE DE LIGNES DU TEXTE
C TEXTE  : LES LIGNES DE TEXTE A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1996
C23456---------------------------------------------------------------012
      CHARACTER*(*)   TEXTE(NBL)
C
C     LE NOMBRE MAXIMAL DE CARACTERES NON BLANCS
      MXC = MXCANB( NBL, TEXTE )
C
C     LE NOMBRE MAXIMAL DE PIXELS OCCUPES PAR LES LIGNES DU TEXTE
      R = LAMXPXTXT( NBL, TEXTE )
C
C     4 CARACTERES DE PLUS POUR LE NUMERO D'OPTION
      IF( R .GT. 0 ) THEN
         LAPXMENU = NINT( R * (4.0+MXC) / MXC )
      ELSE
         LAPXMENU = 0
      ENDIF
      END
