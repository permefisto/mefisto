      SUBROUTINE LLIGNE( KNOMLG , DLIGNE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA LONGUEUR DES ARETES DE LA LIGNE KNOMLG
C -----
C
C ENTREES:
C --------
C KNOMLG : NOM DE LA LIGNE
C
C SORTIES:
C --------
C DLIGNE : LONGUEUR REELLE DOUBLE PRECISION DES ARETES DE LA LIGNE
C          -1D0 SI LIGNE INCONNUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     FEVRIER 1994
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMLG
      DOUBLE PRECISION  DLIGNE
      INTEGER           NOSOEL(4)
C
C     LEXIQUE DE LA LIGNE
      IF( NTLIGN .LE. 0 ) GOTO 9999
      CALL LXLXOU( NTLIGN, KNOMLG, NTLIG, MN )
      IF( NTLIG .LE. 0 ) GOTO 9999
C
C     TMS NSEF DE LA LIGNE
      CALL LXTSOU( NTLIG, 'NSEF',  NTARLG, MNARLG )
      IF( NTARLG .LE. 0 ) GOTO 9999
C
C     TMS XYZSOMMET DE LA LIGNE
      CALL LXTSOU( NTLIG, 'XYZSOMMET', NTSOM, MNSOM )
      IF( NTSOM .LE. 0 ) GOTO 9999
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNARLG) ,
     %             NUTYMA , NBSOEL , NBSOEF , NBTGEF,
     %             LDAPEF , LDNGEF , LDTGEF, NBEFOB ,
     %             NX     , NY     , NZ     ,
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     LA BOUCLE SUR LES ARETES DU MAILLAGE
C     ------------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      DLIGNE = 0D0
      DO 100 N=1,NBEFOB
C
C        LE NUMERO DES 2 SOMMETS DE L'ARETE N
         CALL NSEFNS( N      , NUTYMA , NBSOEF , NBTGEF,
     %                LDAPEF , LDNGEF , LDTGEF ,
     %                MNARLG , NX , NY , NZ ,
     %                NCOGEL , NUGEEF , NUEFTG, NOSOEL , IERR )
C
C        LES ADRESSES DES COORDONNEES DES 2 SOMMETS DE L'ARETE N
         MN1 = MNSOM + WYZSOM + 3 * NOSOEL(1) - 3
         MN2 = MNSOM + WYZSOM + 3 * NOSOEL(2) - 3
C
C        DISTANCE ENTRE CES 2 POINTS
         DLIGNE = DLIGNE + SQRT(
     %                     (RMCN(MN1  )-RMCN(MN2  ))**2
     %                   + (RMCN(MN1+1)-RMCN(MN2+1))**2
     %                   + (RMCN(MN1+2)-RMCN(MN2+2))**2 )
C
 100  CONTINUE
      RETURN
C
C     ERREUR
 9999 DLIGNE = -1D0
      END
