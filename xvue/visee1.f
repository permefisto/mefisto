      SUBROUTINE VISEE1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  ACTIVER LES PARAMETRES DE / TRVARI / POUR LA VISEE DE L'OBJET
C -----  TOUS LES PARAMETRES SONT SUPPOSES CORRECTS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C ......................................................................
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      IF( NOTYVI .EQ. 0 ) THEN
C
C        VISEE PAR DEFAUT
C        ================
         CALL VISEE0
         GOTO 1000
C
      ENDIF
C
      IF( NOTYVI .EQ. 1 ) THEN
C
C        PARAMETRES DE LA VISEE INITIALISES
C        ==================================
         IF( NDIMLI .EQ. 1 ) THEN
C
C           VISEE 1D
C           --------
C           L'ECRAN VU EN COORDONNEES OBJET 2D EN UNITES UTILISATEUR
            CALL FENETRE( AXOPTV(1)-AXOLAR, AXOPTV(1)+AXOLAR,
     %                    AXOPTV(2)-AXOHAU, AXOPTV(2)+AXOHAU )
C
         ELSE
C
C           VISEE 2D
C           --------
C           L'ECRAN VU EN COORDONNEES OBJET 2D EN UNITES UTILISATEUR
            CALL ISOFENETRE( AXOPTV(1)-AXOLAR, AXOPTV(1)+AXOLAR,
     %                       AXOPTV(2)-AXOHAU, AXOPTV(2)+AXOHAU )
         ENDIF
C
      ELSE IF( NOTYVI .EQ. 11 ) THEN
C
C        VISEE 3D: AXONOMETRIE PAR AXOPTV AXOEIL ET AXOARR AXOAVA
C        --------------------------------------------------------
C        LA MATRICE DE L'AXONOMETRIE
         CALL MATAXO
C
      ELSE
C
C        AUCUNE AUTRE AXONOMETRIE N'EST PROGRAMMEE
C        -----------------------------------------
         WRITE(IMPRIM,*) 'VISEE1:',NOTYVI,' CODE VISEE NON PROGRAMME'
C
      ENDIF
C
C     LA MEMOIRE PIXELS EST EFFACEE POUR UN NOUVEAU TRACE
C     ===================================================
ccc 1000 CALL MEMPXFENETRE    09/12/2015

 1000 RETURN
      END
