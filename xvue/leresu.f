      SUBROUTINE LERESU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPRIMER OU TRACER UN RESULTAT STOCKE DANS KERR
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           NAT(1)
C
C     SUR L'IMPRIMANTE
C     ================
      NBL = NBLGRC(NRERR)
      IF( NBL .LE. 0 ) RETURN
      DO 10 I=1,NBL
         WRITE(IMPRIM,'(1X,A)' ) KERR(I)
 10   CONTINUE

      IF( INTERA .LE. 1 ) RETURN
ccc      IF( INTERA .LT. 3 ) RETURN
C
cccC     POUR TRACER ET AFFICHER LE RESULTAT DANS TOUS LES CAS
ccc      INTER0 = INTERA
ccc      IF( INTERA .LT. 3 ) THEN
ccc         INTERA = 3
ccc      ENDIF
C
C     SUR L'ECRAN
C     ===========
C     L'ANCIEN RESULTAT EST EFFACE
      CALL RECTEF( NRERR0 )
C
C     LE NUMERO MAXIMAL DU DERNIER CARACTERE NON BLANC DE KERR
      MDRECT(NRERR) = MXCANB( NBLGRC(NRERR) , KERR )
      IF( MDRECT(NRERR) .GT. 0 ) THEN
C        AJOUT DE 'RESULTAT:' SI CE N'EST PAS FAIT
         L = NUDCNB( KERR(1) )
         IF( LANGAG .EQ. 0 ) THEN
            IF( INDEX( KERR(1), 'Resultat:' ) .NE. 1 ) THEN
               KERR(1) = 'Resultat:' // KERR(1)(1:L)
            ENDIF
         ELSE
            IF( INDEX( KERR(1), 'Result:' ) .NE. 1 ) THEN
               KERR(1) = 'Result:' // KERR(1)(1:L)
            ENDIF
         ENDIF
         MDRECT(NRERR) = MXCANB( NBLGRC(NRERR) , KERR )
      ENDIF
C
C     LE COIN SUPERIEUR GAUCHE ET LA LARGEUR ET HAUTEUR
      LAPX = LAMXPXTXT( NBLGRC(NRERR), KERR )
      DXRECT(NRERR) = ECARLR(NRERR) * 2 + LAPX
      XRECT (NRERR) = LAPXFE - 2*ECARRC - DXRECT(NRERR)
C
      DYRECT(NRERR) = DYLGRC(NRERR) * NBLGRC(NRERR)
      YRECT (NRERR) = ECARRC
C
      IF( DXRECT(NRERR) .GT. 0 .AND. DYRECT(NRERR) .GT. 0 ) THEN
C
C        LES COULEURS DU TRACE : FOND BORD CARACTERES
         NFRECT(NRERR) = NCCYAN
         NBRECT(NRERR) = NCVERT
         NKRECT(NRERR) = NCGRIS
C
C        LE TRACE DU TEXTE
         CALL RECTTX( NRERR, KERR, 0, NAT )
C
C        COPIE DES PARAMETRES POUR LE FUTUR EFFACEMENT
         NBLGRC( NRERR0 ) = NBLGRC( NRERR )
         NFRECT( NRERR0 ) = NFRECT( NRERR )
         NBRECT( NRERR0 ) = NBRECT( NRERR )
         NKRECT( NRERR0 ) = NKRECT( NRERR )
         XRECT ( NRERR0 ) = XRECT ( NRERR )
         YRECT ( NRERR0 ) = YRECT ( NRERR )
         DXRECT( NRERR0 ) = DXRECT( NRERR )
         DYRECT( NRERR0 ) = DYRECT( NRERR )
      ENDIF
C
ccc      INTERA = INTER0
      RETURN
      END
