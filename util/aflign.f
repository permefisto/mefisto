      SUBROUTINE AFLIGN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER LE CONTENU DE LA LIGNE KLIGNE ET METTRE A JOUR LCLIGN
C ----- LE POINTEUR SUR LE DERNIER CARACTERE OCCUPE DE LA KLIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC   PARIS      MARS 1990
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
C.......................................................................
      COMMON /UNITES/LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NOAFTS,NUNIT(26)
C
C     LA LIGNE EST BUFFERISEE A L'AIDE DU POINTEUR LCLIGN DERNIER
C     CARACTERE ENTRE DANS LA LIGNE
      IF( IMPRES*NOMUET .LE. 0 ) RETURN
      IF( KLIGNE .EQ. ' ' ) GOTO 1000
C
c     REMPLACEMENT DES '' PAR '
 10   LCLIGN = INDEX( KLIGNE , '''' // '''' )
      IF( LCLIGN .GT. 0 ) THEN
         LCLIGN = LCLIGN + 1
         KLIGNE( LCLIGN:LCLIGN ) = ' '
         GOTO 10
      ENDIF
C
      IF( NOAFTS .GT. 0 ) THEN
C         AFFICHAGE EXPLICITE D'UN TMS
C         LE FORMAT EVITE UN SECOND BLANC EN TETE DE LIGNE
          WRITE(IMPRIM,10010) KLIGNE
10010 FORMAT(A)
C
      ELSE IF( INTERA .GE. 3 ) THEN
C
C        MENU ET SOURIS
C        ==============
         IF( LHLECT .NE. 1 ) GOTO 1000
C
C        EN INTERACTIF => CONSTRUCTION DU MENU
         DO 15 N=1,NCLIGN
C           REMPLACEMENT DE ' PAR UN BLANC
            IF( KLIGNE(N:N) .EQ. '''' ) KLIGNE(N:N)=' '
 15      CONTINUE
         DO 20 N=NCLIGN,1,-1
C           RECHERCHE DU DERNIER CARACTERE NON BLANC
            IF( KLIGNE(N:N) .NE. ' ' ) GOTO 30
  20     CONTINUE
C        LA LIGNE EST BLANCHE
         GOTO 1000
C
C        COMPRESSION DE KLIGNE  RECHERCHE DES ' '
  30     LCLIGN = INDEX( KLIGNE(1:N) , '  ' )
         IF( LCLIGN .GT. 0 ) THEN
            KLIGNE(LCLIGN:NCLIGN) = KLIGNE(LCLIGN+1:N)
            N = N - 1
            GOTO 30
         ENDIF
C        ICI PAS 2 BLANCS DE SUITE
C
         IF( NBLGRC(NRMENU) .GE. MXLGME ) THEN
C           CONTROLE DE DEBORDEMENT DES LIGNES DE KMENU
C           SI DEBORDEMENT RETOUR A ZERO
            NBLGRC(NRMENU) = 0
         ENDIF
C
C        N EST LE NOMBRE DE CARACTERES A AJOUTER
         IF( NBLGRC(NRMENU) .LE. 0 ) THEN
C           UNE LIGNE EST AJOUTEE
            NBLGRC(NRMENU) = 1
            MDLGRC(NRMENU) = N
            KMENU (NBLGRC(NRMENU)) = KLIGNE(1:N)
            NAMENU(NBLGRC(NRMENU)) = INMENU
         ELSE
C           KMENU (1:MDLGRC(NRMENU)) DEJA ENREGISTRES
            IF( MDLGRC(NRMENU) + N .LT. NBCAME ) THEN
C             IL Y A DE LA PLACE AU DELA DE MDLGRC(NRMENU)
              KMENU(NBLGRC(NRMENU))(MDLGRC(NRMENU)+1:NBCAME)=KLIGNE(1:N)
              MDLGRC(NRMENU) = MDLGRC(NRMENU) + 1 + N
            ELSE
C              IL N'Y A PAS ASSEZ DE PLACE AU DELA DE MDLGRC(NRMENU)
C              UNE LIGNE EST AJOUTEE
               NBLGRC(NRMENU) = NBLGRC(NRMENU) + 1
               MDLGRC(NRMENU) = N
               KMENU (NBLGRC(NRMENU)) = KLIGNE(1:N)
               NAMENU(NBLGRC(NRMENU)) = INMENU
            ENDIF
         ENDIF
C
C        SUPPRESSION DE L'EVENTUEL BLANC EN PREMIER CARACTERE
         IF( KMENU(NBLGRC(NRMENU))(1:1) .EQ. ' ' ) THEN
            KLIGNE = KMENU(NBLGRC(NRMENU))(2:MDLGRC(NRMENU))
            KMENU(NBLGRC(NRMENU)) = KLIGNE
            MDLGRC(NRMENU) = MDLGRC(NRMENU) - 1
         ENDIF
      ENDIF
C
C     LE POINTEUR LCLIGN
 1000 LCLIGN = 0
C     LA LIGNE EST RENDUE BLANCHE
      KLIGNE = ' '
      END
