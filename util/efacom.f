      SUBROUTINE EFACOM( NLMIN  , NCMIN  , NLMAX  , NCMAX , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EFFACER DANS KLG LES COMMENTAIRES { ... } A PARTIR DE
C -----   (NLMIN,NCMIN) JUSQU'A (NLMAX,NCMAX)
C
C ENTREES :
C ---------
C NLMIN,NCMIN : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C NLMAX,NCMAX : POSITION DANS KLG DU DERNIER CARACTERE A TRAITER
C
C SORTIE :
C --------
C NRETOU  : =0 SI PAS D'ERREUR
C           >0 SI UNE ERREUR A ETE RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      NL0 = NLMIN
      NC0 = NCMIN
C
  20  IF( NL0 .LT. NLMAX ) THEN
         NC1 = NBCALI
      ELSE
         NC1 = NCMAX - 1
      ENDIF
      I = INDEX( KLG(NL0)(NC0:NC1) , '{' )
      IF( I .GT. 0 ) THEN
C        DEBUT DE COMMENTAIRE . RECHERCHE DU } FINAL
         J = INDEX( KLG(NL0)(NC0+I:NC1) , '}' )
         IF( J .GT. 0 ) THEN
C           LES CARACTERES NC0-1+I A NC0+I+J-1 SONT BLANCHIS
            KLG(NL0)(NC0-1+I:NC0+I-1+J) = ' '
C           RECHERCHE AU DELA DE CE CARACTERE
            NC0 = NC0+I+J-1
            GOTO  20
         ELSE
C           COMMENTAIRE INCOMPLET
C           LES CARACTERES NC0-1+I A NC1 SONT BLANCHIS
            KLG(NL0)(NC0-1+I:NC1) = ' '
C           TRAITEMENT SUR LA LIGNE SUIVANTE
  40        NL0 = NL0 + 1
            IF( NL0 .GT. NLMAX ) THEN
C              PLUS DE COMMENTAIRE ET COMMENTAIRE INACHEVE
               NBLGRC(NRERR) = 1
               KERR(1) =  'LU: } OUBLIEE'
               CALL LEREUR
               GOTO 9000
            ENDIF
            NC0 = 1
            IF( NL0 .LT. NLMAX ) THEN
               NC1 = NBCALI
            ELSE
               NC1 = NCMAX - 1
            ENDIF
            J = INDEX( KLG(NL0)(NC0:NC1) , '}' )
            IF( J .GT. 0 ) THEN
               NC1 = NC0 + J - 1
            ENDIF
C           LES CARACTERES NC0 A NC1 SONT BLANCHIS
            KLG(NL0)(NC0:NC1) = ' '
            IF( J .EQ. 0 ) GOTO  40
C           RECHERCHE AU DELA DE CE CARACTERE
            NC0 = NC1
            GOTO  20
         ENDIF
      ENDIF
C
C     EXISTE T IL D'AUTRES LIGNES A TRAITER ?
      IF( NL0 .LT. NLMAX ) THEN
C        OUI
         NL0 = NL0 + 1
         NC0 = 1
         GOTO 20
      ENDIF
C
C     TOUTES LES LIGNES SONT TRAITEES
      NRETOU = 0
      RETURN
C
C     ERREUR } OUBLIEE
 9000 NRETOU = 1
      END
