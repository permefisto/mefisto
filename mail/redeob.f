      SUBROUTINE REDEOB( KNOMOB, NUOB, MNOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NOM DU DERNIER OBJET UTILISE
C -----    ET LE TRACER SI MODE INTERACTIF GRAPHIQUE AVEC X11
C
C SORTIES:
C --------
C KNOMOB : NOM DU DERNIER OBJET DECLARE ET 'INCONNU' ou 'UNKNOWN'' SINON
C NUOB   : NUMERO DU NOM DE L'OBJET DANS KNOM LE LEXIQUE S'IL EN EXISTE UN
C          =0 SI PAS DE NOM DANS CE LEXIQUE
C MNOB   : ADRESSE MCN DU LEXIQUE DE L'OBJET
C          =0 SI LEXIQUE INCONNU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC  &  St PIERRE du PERRAY  JUILLET 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/xyzext.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*24      KNOMOB
C
C     RETROUVER LE NOM DU DERNIER OBJET UTILISE DANS NTOBJE
      INIEXT = 0
      MNOB   = 0
      CALL NOMDER( NTOBJE,  KNOMOB, NUOB, MNOB )
ccc      print *,'redeob: ntobje=',ntobje,' knomob=',knomob,
ccc     %        ' nuob=',nuob,' mnob=',mnob
C
      IF( MNOB .LE. 0 ) THEN
C
C        PAS DE LEXIQUE DES OBJETS
         NBC = NUDCNB(KNOMOB)
         NBLGRC(NRHIST) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KHIST(1) = 'OBJET:' // KNOMOB(1:NBC)
         ELSE
            KHIST(1) = 'OBJECT:' // KNOMOB(1:NBC)
         ENDIF
         CALL LHISTO
C
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: AUCUN OBJET CONSTRUIT'
            KERR(2) = 'A CREER AVEC LA COMMANDE MAILLER'
         ELSE
            KERR(1) = 'ERROR: ZERO DECLARED OBJECT'
            KERR(2) = 'CREATE IT WITH MESHER COMMAND'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
      IF( NUOB .GT. 0 ) THEN
C
C        TRACE DE L'OBJET POUR VOIR CE QUI EXISTE
         CALL LXNLOU( NTOBJE, NUOB, NTLXOB, MNLXOB )
C
C        RECHERCHE DES COORDONNEES DES POINTS
         CALL LXTSOU( NTLXOB, 'XYZPOINT', NTXYZP, MNXYZP )
         IF( NTXYZP .LE. 0 ) THEN
            CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZP, MNXYZP )
            IF( NTXYZP .LE. 0 ) THEN
               NBLGRC(NRERR) = 3
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'ERREUR: OBJET ' // KNOMOB
                  KERR(2) = 'SANS TMS XYZSOMMET et XYZPOINT'
                  KERR(3) = 'A CREER AVEC MAILLER'
               ELSE
                  KERR(1) = 'ERROR: OBJECT ' // KNOMOB
                  KERR(2) = 'WITHOUT TMS XYZSOMMET and XYZPOINT'
                  KERR(3) = 'CREATE IT WITH MAILLER'
               ENDIF
               CALL LEREUR
               GOTO 9900
            ENDIF
         ENDIF
         CALL MAJEXT( MNXYZP )
C
C        TRACE DU DERNIER OBJET
         IF( INTERA .GE. 1 ) THEN
            CALL EFFACEMEMPX
            CALL ITEMS0
            CALL VISEE0
C
C           TYPE OBJET ET SON NOM
            NBC = NUDCNB(KNOMOB)
            NBLGRC(NRHIST) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KHIST(1) = 'OBJET:' // KNOMOB(1:NBC)
            ELSE
               KHIST(1) = 'OBJECT:' // KNOMOB(1:NBC)
            ENDIF
            CALL LHISTO
C
C           TRACE ET AFFICHAGE DE LA QUALITE DU MAILLAGE DE L'OBJET 
C           mettre ccc au dessous pour recuperer plus vite l'execution des 6cubes
C           TRACE DE LA QUALITE DES EF et TRANSLATION ROTATION ZOOM
            LORBITE = 0
            CALL T1OBJE( KNOMOB )
         ENDIF
C
      ENDIF
C
 9900 RETURN
      END
