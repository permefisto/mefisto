      SUBROUTINE QUIFIC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  QUITTER LE LECTEUR SUR FICHIER POUR REVENIR AU LECTEUR CLAVIER
C -----  SAUF EN CAS DE BATCH OU LE CALCUL EST ARRETE
C++++++++++++++++++++++++ +++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1992
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       MARS 1996
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      IF( LHLECT.GT. 1 .AND. LHLECT.LE.MXLECT ) THEN
C
C        LECTURE ACTUELLE SUR UN FICHIER AVEC UNE ERREUR
C        => FERMETURE DE TOUS LES FICHIERS OUVERTS ET RETOUR AU CLAVIER
C        MODIFICATION OCTOBRE 1994
         DO I=LHLECT,2,-1
            CLOSE( LPLECT(I) )
         ENDDO

         IF( INTERA .LE. 1 ) THEN
C            EN LECTURE DE FICHIER ARRET DU CALCUL
             CALL ARRET( 100 )
         ENDIF

C        AFFICHAGE DU RETOUR AU CLAVIER
         LHLECT = 1
         LECTEU = LPLECT( LHLECT )
C        LE NIVEAU D'INTER-ACTION EST REGENERE
         INTERA = INTERB( LHLECT )
C
         IF( INTERA .GE. 3 ) THEN
C           LE RECTANGLE ERREUR EST EFFACE
            CALL RECTEF( NRERR  )
C           LE RECTANGLE INVITE EST EFFACE
            CALL RECTEF( NRINVI )
C           LE RECTANGLE MENU EST EFFACE
            CALL RECTEF( NRMENU )
         ENDIF
C
C        REMISE A ZERO DES BUFFERS
         NBLGRC(NRMENU) = 0
         NBLGRC(NRLGLU) = 0
         NBLGRC(NRINVI) = 0
         NBLGRC(NRLGSA) = 0
C
C        L'AFFICHAGE DE L'ERREUR
         NBLGRC(NRERR ) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'DONNEE INCORRECTE LUE DANS LE FICHIER '
            KERR(2) = 'CORRIGER LE FICHIER ET LE SAUVEGARDER '
            KERR(3) = 'PUIS, RELANCER PAR LA FRAPPE DE '
            KERR(4) = 'liref nom_du_fichier;'
         ELSE
            KERR(1) = 'INCORRECT DATA READ ON THIS FILE'
            KERR(2) = 'RECTIFY the FILE and SAVE IT'
            KERR(3) = 'THUS, Type again'
            KERR(4) = 'readf File_Name ;'
         ENDIF
C
C        LES DONNEES REPARTENT DE ZERO
         NLPTV1 = 0
      ENDIF
CCCC
CCCC     =====================================================================
CCCC         PAS DE SORTIE DE LA LECTURE DU FICHIER
CCCC     =====================================================================
CCCC     L'AFFICHAGE DE L'ERREUR
CCC      NBLGRC(NRERR ) = 2
CCC      KERR(1) = 'DONNEE INCORRECTE LUE DANS LE FICHIER'
CCC      KERR(2) = 'CORRIGER LE FICHIER'

      RETURN
      END
