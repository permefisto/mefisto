      SUBROUTINE ARRET( NOCODE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SAUVEGARDER LES TMS ET AUTRES DONNEES
C -----
C ENTREE :
C --------
C NOCODE : =  0 ARRET DEMANDE POUR FINIR NORMALEMENT L'EXECUTION
C          >  0 ARRET DEMANDE PAR CAUSE D'ERREUR
C          = -1 PAS DE SAUVEGARDE DES DONNEES . ARRET IMMEDIAT
C          < -1 SAUVEGARDE DES DONNEES ET RETOUR AU CALCUL
C               UNE REINITIALISATION EST OBLIGATOIRE ( CF OUVRMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MARS 1990
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      COMMON /UNITES/ LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
C     ESSAI POUR ENRAYER LA RECURSIVITE
      INTEGER         NBPASS
      SAVE            NBPASS
      DATA            NBPASS / 0 /
C
      IF( NBPASS .GT. 0 ) STOP
C
C     NOMBRE DE PASSAGE DANS ARRET
      NBPASS = NBPASS + 1
C
C     AFFICHAGE DE LA SAUVEGARDE OU DE L'ARRET D'URGENCE
C     ==================================================
      IF( NOCODE .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'CORRECTE SAUVEGARDE des DONNEES sur FICHIERS ...'
         ELSE
            KERR(1) = 'OK, All DATA are SAVING on DISK FILES ...'
         ENDIF
         CALL LERESU
      ELSE
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'Desole, ARRET suite a une ERREUR DETECTEE'
            KERR(2) = 'ESSAI de SAUVEGARDE des DONNEES sur DISQUE ...'
         ELSE
            KERR(1) = 'SORRY, STOP after an ERROR DETECTION'
            KERR(2) = 'TRYing to SAVE DATA on DISK FILES ...'
         ENDIF
         CALL LEREUR
      ENDIF
      WRITE(IMPRIM,*)
C
C     ARRET IMMEDIAT ?
C     ================
      IF( NOCODE .EQ. -1 ) STOP
C
C     SAUVEGARDE DES VARIABLES DU LANGAGE UTILISATEUR
C     ===============================================
      CALL LUFE
C
C     LA MS EST FERMEE
C     ================
      CALL MSFE
C
C     FERMETURE DE TOUTES LES RESSOURCES X11 UTILISEES PAR MEFISTO
C     ============================================================
      IF( INTERA .GE. 1 ) CALL XVFERMER
C
C     FERMETURE DU FICHIER FRAPPE
      CLOSE( UNIT=NFFRAP, STATUS='KEEP' )
C
      IF( NOCODE .EQ. 0 ) THEN
C
C        LE TRAITEMENT EST NORMALEMENT ARRETE
C        ====================================
         IF( LANGAG .EQ. 0 ) THEN
            STOP ': FIN CORRECTE d''EXECUTION MEFISTO'
         ELSE
            STOP ': CORRECT END of MEFISTO EXECUTION'
         ENDIF
C
      ELSE IF( NOCODE .GT. 0 ) THEN
C
C        STOP A AJUSTER SELON LE SITE
         IF( LANGAG .EQ. 0 ) THEN
            STOP 'Desole, CRASH du logiciel MEFISTO'
         ELSE
            STOP 'Sorry, CRASH of MEFISTO software'
         ENDIF
CCC         CALL TILT

      ENDIF
      WRITE(IMPRIM,*)
C
C     APRES SAUVEGARDE DE LA MS RETOUR AU CALCUL
C     ==========================================
      NBPASS = 0
      RETURN
      END
