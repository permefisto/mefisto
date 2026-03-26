      SUBROUTINE LINVIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES LIGNES D'INVITE DANS UN RECTANGLE GRAPHIQUE
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1995
C2345X7--------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      CHARACTER*(NBCAIN+8) KINVIT(MXLGIN)
      INTRINSIC     INT
C
      IF( LHLECT .EQ. 1 ) THEN
C
C        LECTURE MANUELLE
C        LE NOMBRE MAXIMAL DE CARACTERES NON BLANCS
         MDRECT(NRINVI) = MXCANB( NBLGRC(NRINVI) , KINVI )
C
C        SELON L'INTERACTIVITE TRACE OU AFFICHAGE
         IF( INTERA .GE. 3 ) THEN
C
C           SI PAS DE CARACTERES NON BLANCS DANS TEXTE: RETOUR
            IF( MDRECT(NRINVI) .LE. 0 ) THEN
               NBLGRC(NRINVI) = 0
               DXRECT(NRINVI) = 0
                XRECT(NRINVI) = 0
               DYRECT(NRINVI) = 0
                YRECT(NRINVI) = 0
               RETURN
            ENDIF
C
C           AJOUT DES CARACTERES 'INVITE:' EN DEBUT DE LIGNE
            IF( LANGAG .EQ. 0 ) THEN
C              FRANCAIS
               KINVIT(1) = 'INVITE:' // KINVI(1)
               DO 5 I=2,NBLGRC(NRINVI)
                  KINVIT(I) = 'INVITE:' // KINVI(I)
 5             CONTINUE
            ELSE
C              ENGLISH
               KINVIT(1) = 'PROMPT:' // KINVI(1)
               DO 6 I=2,NBLGRC(NRINVI)
                  KINVIT(I) = 'PROMPT:' // KINVI(I)
 6             CONTINUE
            ENDIF
C
            LAPX = LAMXPXTXT( NBLGRC(NRINVI), KINVIT )
            DXRECT(NRINVI) = ECARLR(NRINVI) * 2 + LAPX
CCC            DXRECT(NRINVI) = MAX( DXRECT(NRERR),
            DXRECT(NRINVI) = MAX( DXRECT(NRMENU),
     %                            DXRECT(NRLGLU),
     %                            DXRECT(NRINVI) )
             XRECT(NRINVI) = LAPXFE - 2*ECARRC - DXRECT(NRINVI)
C
            DYRECT(NRINVI) = NBLGRC(NRINVI) * DYLGRC(NRINVI)
            IF( DYRECT(NRINVI) .GT. 0 ) THEN
              IF( DYRECT(NRLGLU) .GT. 0 ) THEN
              YRECT(NRINVI)=INT(YRECT(NRLGLU)+DYRECT(NRLGLU)+1.5*ECARRC)
              ELSE IF( DYRECT(NRMENU) .GT. 0 ) THEN
              YRECT(NRINVI)=INT(YRECT(NRMENU)+DYRECT(NRMENU)+1.5*ECARRC)
              ELSE IF( DYRECT(NRERR) .GT. 0 ) THEN
              YRECT(NRINVI)=INT(YRECT(NRERR) +DYRECT(NRERR) +1.5*ECARRC)
              ELSE
                 YRECT(NRINVI)=ECARRC
              ENDIF
            ELSE
               YRECT(NRINVI) = 0
            ENDIF
C
C           LE TRACE DU TEXTE
            CALL RECTTX( NRINVI, KINVIT, 0, NA )
C
C           SAUVEGARDE DANS LE FICHIER FRAPPE
            DO 10 I=1,NBLGRC( NRINVI )
               KINVIT(MXLGIN)='{ ' // KINVI(I)(1:MDRECT(NRINVI)) // ' }'
               CALL SANSDBL( KINVIT(MXLGIN), L )
               WRITE(NFFRAP,*) KINVIT(MXLGIN)(1:L)
 10         CONTINUE
C
         ELSE IF( INTERA .GT. 0 ) THEN
C
C           AFFICHAGE
            DO 20 I=1,NBLGRC( NRINVI )
               WRITE(IMPRIM,*) KINVI(I)(1:MDRECT(NRINVI))
 20         CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
