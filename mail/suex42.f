      SUBROUTINE SUEX42( NTLXSU, LADEFI,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PROJECTION DES TANGENTES D'UN MEME SOMMET SUR LE PLAN
C -----    A DISTANCE MINIMALE DES EXTREMITES DES TANGENTES
C          => SURFACE MAILLEE G1-CONTINUE SI TANGENTES EN TOUT SOMMET
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF $MEFISTO/td/d/a_surface__definition
C
C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C          CF $MEFISTO/td/d/a___nsef
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF $MEFISTO/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        JUIN 1999
C ...................................................................012
      include"./incl/a___xyzsommet.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      INTEGER        LADEFI(0:*)
C
C     LA SURFACE A TRAITER
C     ====================
C     LE NUMERO DE CETTE SURFACE DANS LE LEXIQUE DES SURFACES
      NUSUIN = LADEFI(WUSUIN)
      CALL SUTOTG( NUSUIN, NTLXSU,
     %             NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     PROJECTION DES TGS SUR LE PLAN A DISTANCE MINIMALE
C     POUR TOUS LES SOMMETS DU MAILLAGE DE LA SURFACE
C     ==================================================
      CALL SUC0G1( NTLXSU, 0,      NUSUIN,
     %             NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      IF( IERR .EQ. 0 .AND. MNXYZS .GT. 0 ) MCN( MNXYZS + WBCOOR ) = 3
C
      RETURN
      END
