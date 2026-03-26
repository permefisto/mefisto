      SUBROUTINE AFFLUXFR( KNOMOB )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   IMPRESSION DES NDSM FLUX DE CHALEUR A TRAVERS LES FRONTIERES
C -----   AVEC CONDITION AUX LIMITES DE L'OBJET C'EST A DIRE
C         EN 1D : FLUX AUX 2 POINTS EXTREMITES DE L'OBJET SI AVEC CL
C         EN 2D : FLUX INTEGRE SUR CHAQUE LIGNE   DE L'OBJET
C         EN 3D : FLUX INTEGRE SUR CHAQUE SURFACE DE L'OBJET (=>CL)
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET DE FLUX A AFFICHER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC PARIS St Pierre du Perray JUIN 2009
C23456---------------------------------------------------------------012
      include"./incl/ntmnlt.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___fluxfr.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      CHARACTER*(*)     KNOMOB
      CHARACTER*24      KNOM
      CHARACTER*8       KTYPE,KTYPEA
C
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 9000
      ENDIF
C
      CALL LXLXOU( NTLXOB, 'FLUXFR', NTFLFR, MNFLFR )
      IF( NTFLFR .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET ' // KNOMOB // ' SANS FLUXFR'
         ELSE
            KERR(1) = 'OBJECT ' // KNOMOB // ' WITHOUT FLUXFR'
         ENDIF
         CALL LEREUR
         GOTO 9000
      ENDIF
C
C     NOMBRE DE CAS CALCULES
      NBCAFF = MCN( MNFLFR + WBCAFF )
C     NUMERO DU TYPE DES FRONTIERES: 1 POINT 2 LIGNE 3 SURFACE
      NTPLSF = MCN( MNFLFR + WTPLSF )
C     NOMBRE DE FRONTIERES PLSF
      NBPLSF = MCN( MNFLFR + WBPLSF )
C
      IF( NTPLSF .EQ. 3 ) THEN
         KTYPE  = 'SURFACE '
         KTYPEA = 'SURFACE '
      ELSE IF( NTPLSF .EQ. 2 ) THEN
         KTYPE  = 'LIGNE '
         KTYPEA = 'LINE '
      ELSE IF( NTPLSF .EQ. 1 ) THEN
         KTYPE  = 'POINT '
         KTYPEA = 'POINT '
      ENDIF
C
      MOREE2 = MOTVAR(6)
      MOFLFR = MOREE2 * NBPLSF * NBCAFF
      MNNO   = MNFLFR + WLUXFR + MOFLFR - 1
      MNFLUX = ( MNFLFR + WLUXFR - 1 ) / MOREE2
      DO 20 N = 1, NBPLSF
C
C        LE NOM DU PLS
         CALL NMOBNU( KTYPE, MCN(MNNO+N), KNOM )
C
C        AFFICHAGE DU FLUX A TRAVERS CETTE SURFACE ou LIGNE ou POINT
         MN = MNFLUX + N - NBPLSF
         IF( LANGAG .EQ. 0 ) THEN
C           AFFICHAGE EN FRANCAIS
            IF( NBCAFF .LE. 5 ) THEN
               WRITE(IMPRIM,10020) KTYPE, KNOM,
     %               (K, DMCN(MN+K*NBPLSF),K=1,NBCAFF)
            ELSE
               WRITE(IMPRIM,10021) KTYPE, KNOM,
     %               (K, DMCN(MN+K*NBPLSF),K=1,NBCAFF)
            ENDIF
         ELSE
C           AFFICHAGE EN ANGLAIS
            IF( NBCAFF .LE. 5 ) THEN
              WRITE(IMPRIM,20020) KTYPEA, KNOM,
     %               (K, DMCN(MN+K*NBPLSF),K=1,NBCAFF)
            ELSE
              WRITE(IMPRIM,20021) KTYPEA, KNOM,
     %               (K, DMCN(MN+K*NBPLSF),K=1,NBCAFF)
            ENDIF
         ENDIF
C
 20   CONTINUE
C
10020    FORMAT(1X,A,': ',A,' FLUX | ',5('cas',I5,':',G14.6,' |') )
10021    FORMAT(1X,A,': ',A / (' FLUX | ',5('cas',I5,':',G14.6,' |')) )
C
20020    FORMAT(1X,A,': ',A,' FLUX | ',5('case',I5,':',G14.6,' |') )
20021    FORMAT(1X,A,': ',A / (' FLUX | ',5('case',I5,':',G14.6,' |')) )
C
 9000    RETURN
         END
