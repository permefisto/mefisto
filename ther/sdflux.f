      SUBROUTINE SDFLUX( NBDOBJ, NTYOBJ, NDIM, NDSM,
     %                   NUMINI, NUMAXI, FLUXTO , NBLIF , LIGNE , FLUX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LES CONTRIBUTIONS DES FLUX LOCAUX AUX FLUX GLOBAUX
C -----
C
C ENTREES:
C --------
C NBDOBJ : NOMBRE DE PLSV DE L'OBJET
C NTYOBJ : TYPE OBJET DES PLSV DE L'OBJET  (NO TYPE 1 A 4 ET NO DANS LX )
C NDIM   : DIMENSION DE L'ESPACE (2 OU 3)
C NDSM   : NOMBRE DE CARTES DE TEMPERATURE OU CAS TRAITES
C NUMINI : NUMERO MINIMAL DES SURFACES (LIGNES EN 2D)
C NUMAXI : NUMERO MAXIMAL DES SURFACES (LIGNES EN 2D)
C FLUXTO : FLUX(I,J) J-EME FLUX A TRAVERS LA SURFACE DE NUMERO I
C LIGF   : NOMBRE DE LIGNES DE L'OBJET GLOBAL
C LIGNE  : NUMEROS DES LIGNES DE L'OBJET GLOBAL
C
C SORTIES :
C ---------
C FLUX   : LES FLUX ASSOCIES AUX  LIGNES DE L'OBJET GLOBAL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY    ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      DOUBLE PRECISION  FLUXTO( NUMINI:NUMAXI, 1:NDSM )
      DOUBLE PRECISION  FLUX( 1:NBLIF , 1:NDSM )
      INTEGER           NTYOBJ(2,NBDOBJ),LIGNE(NBLIF)
C
      IF( NDIM .NE. 2 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) ='ERREUR: APPEL DE SDFLUX EN DIMENSION 3'
         KERR(2) ='CAS NON PROGRAMME'
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C
      DO 20 N=NUMINI, NUMAXI
C        LE NUMERO DE LA LIGNE A RECHERCHER
         DO 10 I=1,NBDOBJ
            IF( NTYOBJ(2,I) .EQ. N ) THEN
C               IL EXISTE UN PLSV DE NUMERO N DANS SON LEXIQUE
C               EST-IL UNE SURFACE EN 3D OU UNE LIGNE EN 2D?
                NULIRE = N
                IF( NTYOBJ(1,I) .EQ. NDIM ) GOTO 15
            ENDIF
 10      CONTINUE
C        LE LEXIQUE DE LA LIGNE A RETROUVER
 15      CALL LXNLOU( NTLIGN , NULIRE , NTLIRE , MNLIRE )
C        LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
         CALL LXTSOU( NTLIRE , 'XYZSOMMET' , NTSORE , MNSORE )
         NBSORE = MCN(MNSORE+WNBSOM)
         NXYZRE = MNSORE+WYZSOM-1
C
C        RECHERCHE PARMI LES OBJETS AUX LIMITES GLOBAUX
C        ----------------------------------------------
         DO 100 NBL = 1 , NBLIF
            NULITR = LIGNE(NBL)
C           LE LEXIQUE DE CETTE LIGNE
            CALL LXNLOU( NTLIGN , NULITR , NTLITR , MNLITR )
C           LE TABLEAU 'XYZSOMMET' DE CET OBJET
            CALL LXTSOU( NTLITR , 'XYZSOMMET' , NTSOTR , MNSOTR )
            NBSOTR = MCN(MNSOTR+WNBSOM)
            NXYZTR = MNSOTR+WYZSOM-1
C           COMPARAISON DES POINTS
            NB = 0
            DO 101 NBSRE = 1 , NBSORE
               NUM1 = NXYZRE + ( NBSRE - 1 ) * 3
               DO 102 NBSTR = 1 , NBSOTR
                  NUM2 = NXYZTR + ( NBSTR - 1 ) * 3
                  DO 103 J=1,3
                     R1 = RMCN(NUM1+J)
                     R2 = RMCN(NUM2+J)
                     IF     ( ABS(R1) .LE. EPZERO ) THEN
                        IF  ( ABS(R2) .GT. EPZERO ) GOTO 102
                     ELSE IF( ABS(R2) .LE. EPZERO ) THEN
                         GOTO 102
                     ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
                         GOTO 102
                     ENDIF
 103              CONTINUE
C                 LE SOMMET EST RETROUVE
                  NB = NB + 1
                  GO TO 101
 102          CONTINUE
 101        CONTINUE
            IF ( NB .EQ. NBSORE ) THEN
C              TOUS LES SOMMETS SONT RETROUVES
               DO 16 ND = 1 , NDSM
                  FLUX(NBL,ND) = FLUX(NBL,ND) + FLUXTO(N,ND)
 16            CONTINUE
               GOTO 20
            ELSE
C              ON POURSUIT LA RECHERCHE
               GO TO 100
            END IF
 100     CONTINUE
C
 20   CONTINUE
C
      RETURN
      END
