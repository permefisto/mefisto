      SUBROUTINE SDDEF5( NBLIMI, LIMIOB, NUTYOB, NTOBRE, NTOBTR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RECHERCHE D'UN OBJET AUX LIMITES DEFINI DANS UN OBJET
C -----   DANS LA LISTE DES OBJETS DEFINI POUR L'OBJET SD
C          ( RECHERCHE COMME SOUS-ENSEMBLE )
C
C ENTREE :
C --------
C NBLIMI : NOMBRE D'OBJETS AUX LIMITES DEFINIS PUR L'OBJET SD
C LIMIOB : LISTE DE CES OBJETS
C NUTYOB : NUMERO DU TYPE DE L'OBJET RECHERCHE
C NTOBRE : SON NUMERO DE LEXIQUE
C
C SORTIE :
C --------
C NTOBTR : NUMERO DU LEXIQUE DE L'OBJET AUX LIMITES SD QUI LE CONTIENT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  FEVRIER 1990
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
C
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      INTEGER           LIMIOB(2,NBLIMI)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE TABLEAU 'XYZSOMMET' DE L'OBJET A RETROUVER
      CALL LXTSOU( NTOBRE , 'XYZSOMMET' , NTSORE , MNSORE )
      IF( NTSORE .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR SDDEF5: TMS XYZSOMMET NON RETROUVE'
         ELSE
            KERR(1) = 'ERROR SDDEF5: TMS XYZSOMMET UNKNOWN'
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBRE )
         GOTO 9000
      ENDIF
      NBSORE = MCN(MNSORE+WNBSOM)
      NXYZRE = MNSORE+WYZSOM-1
C
C     RECHERCHE PARMI LES OBJETS AUX LIMITES
C     --------------------------------------
      DO 100 NBL = 1 , NBLIMI
         NUTYTR = LIMIOB(1,NBL)
         NTOBTR = LIMIOB(2,NBL)
         IF (NUTYOB .LE. NUTYTR) THEN
C           LE LEXIQUE DE L'OBJET
C           LE TABLEAU 'XYZSOMMET' DE CET OBJET
            CALL LXTSOU( NTOBTR , 'XYZSOMMET' , NTSOTR , MNSOTR )
            IF( NTSOTR .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'ERREUR SDDEF5: TMS XYZSOMMET NON RETROUVE'
               ELSE
                  KERR(1) = 'ERROR SDDEF5: TMS XYZSOMMET UNKNOWN'
               ENDIF
               CALL LEREUR
               CALL LXIM( NTOBTR )
               GOTO 9000
            ENDIF
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
 102           CONTINUE
 101        CONTINUE
            IF ( NB .EQ. NBSORE ) THEN
C             TOUS LES SOMMETS SONT RETROUVES : OBRE EST SOUS-ENSEMBLE DE OBTR
              RETURN
            ELSE
C             TOUS LES SOMMETS NE SONT PAS RETROUVES : ON POURSUIT LA RECHERCHE
              GO TO 100
            END IF
         END IF
 100  CONTINUE
C
C     OBRE N'EST PAS UN SOUS-ENSEMBLE DE OBTR : PASSAGE A L'OBJET SUIVANT
 9000 NTOBTR = 0
      END
