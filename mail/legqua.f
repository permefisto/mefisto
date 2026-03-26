      SUBROUTINE LEGQUA
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA LEGENDE DE LA QUALITE SUR L'ECRAN
C------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1994
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
C     LES VARIABLES LOCALES
      CHARACTER*12 KNOM
C
C     TRACE DE LA LEGENDE DES QUALITES
C     --------------------------------
      IF( INTERA .GE. 1 ) THEN
C
C        EFFACEMENT DE LA LEGENDE DE LA QUALITE SUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
           IF ( LASOPS.EQ.1 ) THEN
             LASOPS = -11
           ELSE
             IF ( LASOPS.EQ.2 ) THEN
               LASOPS = -12
             ELSE
               LASOPS = 0
               NBLGRC(NRERR) = 2
               KERR(1) = 'LEGQUA: MAUVAISE VALEUR DE LASOPS'
               KERR(2) = '        ARRET DU POSTSCRIPT'
               CALL LEREUR
             ENDIF
           ENDIF
           CALL XVPOSTSCRIPT(LASOPS)
           LASOPS = - LASOPS
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
C        TRACE DE LA LEGENDE DE LA QUALITE DES EF
C        ----------------------------------------
         NX   = LAPXFE - 170
         NY   = 500
         KNOM = 'QUALITES EF'
         CALL XVCOULEUR( NCBLAN )
         CALL XVTEXTE( KNOM(1:11), 11, NX, NY )
         NY = NY + 20
         DO 930 II=1,10
            NCOUL = N1COUL + 10 - II
            CALL XVCOULEUR( NCOUL )
            CALL XVRECTANGLE( NX, NY, 30, 10 )
            IF( II .NE. 10 ) THEN
               WRITE( KNOM(3:12), '(G10.2)' ) (10.0 - II ) * 0.1
            ELSE
               KNOM(3:12) = '   .00    '
            ENDIF
            KNOM(1:2) = ' >'
            CALL XVTEXTE( KNOM(1:12), 12, NX+30, NY+10 )
            NY = NY + 20
 930     CONTINUE
C
C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
           LASOPS = LASOPS -10
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
CCCC        BLOCAGE POUR LIRE TRANQUILLEMENT LE TRACE
CCC         CALL CLICSO
      ENDIF

      RETURN
      END
