      SUBROUTINE VOEX54( NUVOIN,
     %                   NTNSEV, MNNSEV, NTXYZV, MNXYZV, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     1 VOLUME MULTI-MATERIAUX, UNION DE PLUSIEURS VOLUMES,
C -----     DEVIENT UN SEUL VOLUME MONO-MATERIAU
C           IL SUFFIT DE DETRUIRE les tableaux a_volume__materiaux et union
C           ATTENTION: Cette OPERATION N'EST PAS REVERSIBLE!
C
C ENTREES:
C --------
C NUVOIN : NUMERO DU VOLUME INITIAL et FINAL DANS LE LEXIQUE DES VOLUMES
C          Le tableau de definition est ici inutile.
C          Le no du volume initial (=final) suffit
C
C SORTIES:
C --------
C NTNSEV : NUMERO      DU TMS 'NSEF' DU VOLUME MONO-MATERIAU
C MNNSEV : ADRESSE MCN DU TMS 'NSEF' DU VOLUME MONO-MATERIAU
C          CF '$MEFISTO/td/d/a___nsef'
C NTXYZV : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME MONO-MATERIAU
C MNXYZV : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME MONO-MATERIAU
C          CF '~/td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC  & St PIERRE du PERRAY   avril 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___materiaux.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*24      NMVOIN
C
C     NUMERO DU VOLUME MULTI-MATERIAUX A TRAITER
      IF( NUVOIN .LE. 0 ) GOTO 9991
C
C     NOM DU VOLUME
      CALL NMOBNU( 'VOLUME', NUVOIN, NMVOIN )

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'voex54: ',NMVOIN,
     %   ' VOLUME MULTI-MATERIAUX DEVIENT un VOLUME MONO-MATERIAU'
      ELSE
         PRINT*,'voex54: ',NMVOIN,
     %   ' MULTI-MATERIAL VOLUME BECOMES a MONO-MATERIAL VOLUME'
      ENDIF
C
C     RECUPERATION DES TABLEAUX DU VOLUME MULTI-MATERIAUX
      CALL LXNLOU( NTVOLU, NUVOIN, NTVOIN, MNVOIN )
      IF( NTVOIN .LE. 0 ) GOTO 9991
C
C     LE TABLEAU a_volume__materiaux EST SUPPRIME
      CALL LXTSOU( NTVOIN, 'MATERIAUX', NTMATE, MNMATE )
      IF( NTMATE .LE. 0 ) THEN
C        UN SEUL MATERIAU => PAS DE SUPPRESSION
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = NMVOIN // ' VOLUME INITIAL DEJA MONO-MATERIAU'
         ELSE
            KERR(1) = NMVOIN // ' INITIAL VOLUME ALREADY MONO-MATERIAL'
         ENDIF
         CALL LERESU
      ELSE
C        SUPPRESSION DU TABLEAU a_volume__materiaux
         CALL LXTSDS( NTVOIN, 'MATERIAUX' )
      ENDIF
C
C     LE TABLEAU a___union EST SUPPRIME
      CALL LXTSOU( NTVOIN, 'UNION', NTUNION, MNUNION )
C     SUPPRESSION DU TABLEAU a_volume__union
      IF( NTUNION .GT. 0 ) CALL LXTSDS( NTVOIN, 'UNION' )
C
C     OUVERTURE DU TABLEAU XYZSOMMET DU VOLUME
      CALL LXTSOU( NTVOIN, 'XYZSOMMET', NTXYZV, MNXYZV )
      IF( NTXYZV .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME INITIAL ' // NMVOIN
            KERR(2) = 'est un VOLUME SANS le TMS XYZSOMMET'
         ELSE
            KERR(1) = 'INITIAL VOLUME NAME: ' // NMVOIN
            KERR(2) = 'is a VOLUME WITHOUT the XYZSOMMET TMS'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ELSE
C        LA DATE DE CREATION DU TABLEAU XYZSOMMET DU VOLUME MONO-MATERIAU
         CALL ECDATE( MCN(MNXYZV) )
      ENDIF
C
C     LE TABLEAU 'NSEF' DU VOLUME
      CALL LXTSOU( NTVOIN, 'NSEF', NTNSEV, MNNSEV )
      IF( NTNSEV .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME INITIAL ' // NMVOIN
            KERR(2) = 'est un VOLUME SANS le TMS NSEF'
         ELSE
            KERR(1) = 'INITIAL VOLUME NAME: ' // NMVOIN
            KERR(2) = 'is a VOLUME WITHOUT the NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ELSE
C        LA DATE DE CREATION DU TABLEAU NSEF DU VOLUME MONO-MATERIAU
         CALL ECDATE( MCN(MNNSEV) )
      ENDIF
C
C     SORTIE SANS ERREUR RENCONTREE
      IERR = 0
      RETURN
C
C     ERREUR
 9991 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = NMVOIN // ' VOLUME INITIAL MULTI-MATERIAUX INCONNU'
      ELSE
         KERR(1) = NMVOIN // ' UNKNOWN INITIAL MULTI-MATERIAL VOLUME'
      ENDIF
      CALL LEREUR
      NTNSEV = 0
      MNNSEV = 0
      NTXYZV = 0
      MNXYZV = 0
      IERR   = 1
C
      RETURN
      END
