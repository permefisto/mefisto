      SUBROUTINE IDLISU( MNXYLI, MNXYSU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER ET IMPOSER LES 3 COORDONNEES DES SOMMETS
C -----    DE LA LIGNE AUX SOMMETS DE LA SURFACE SACHANT QU'ILS
C          DOIVENT TOUS ETRE IDENTIFIES CAR CES DERNIERS SONT
C          SUPPOSES ENTACHES D'UNE ERREUR D'ARRONDI
C
C ENTREES:
C --------
C MNXYLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C MNXYSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C
C SORTIES:
C --------
C MNXYSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          LES 3 COORDONNEES DES SOMMETS DE LA LIGNE SONT IMPOSES
C          AUX SOMMETS LES PLUS PROCHES DE LA SURFACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS        DECEMBRE 1997
C234567..............................................................012
      include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     PAS DE TEST D'EXISTENCE DES TMS
C     ILS SONT SUPPOSES EXISTER ET ETRE CORRECT
C     (CAS DES MAILLAGES ALGEBRIQUES)
C     LE NOMBRE DE SOMMETS DE LA LIGNE ET DE LA SURFACE
      NBSOLI = MCN( MNXYLI + WNBSOM )
      NBSOSU = MCN( MNXYSU + WNBSOM )
C
C     ADRESSE DE DEBUT DES COORDONNEES
      MNSL = MNXYLI + WYZSOM
      MNSS = MNXYSU + WYZSOM
C
C     ALGORITHME BESTIAL EN NBSOLI*NBSOSU OPERATIONS !
C     ================================================
      DO 50 I=1,NBSOLI
C        TRAITEMENT DU SOMMET I DE LA LIGNE
         DMIN  = 1E38
         JDMIN = 0
         MN    = MNSS
         DO 10 J=1,NBSOSU
            D = ( RMCN(MNSL  )-RMCN(MN  ) ) ** 2
     %        + ( RMCN(MNSL+1)-RMCN(MN+1) ) ** 2
     %        + ( RMCN(MNSL+2)-RMCN(MN+2) ) ** 2
            IF( D .LT. DMIN ) THEN
               DMIN  = D
               JDMIN = J
            ENDIF
            MN = MN + 3
 10      CONTINUE
C
C        LE SOMMET JDMIN DE LA SURFACE EST LE SOMMET I DE LA LIGNE
         MN = MNSS + 3 * JDMIN - 3
         RMCN(MN  ) = RMCN(MNSL  )
         RMCN(MN+1) = RMCN(MNSL+1)
         RMCN(MN+2) = RMCN(MNSL+2)
C
C        PASSAGE AU SOMMET SUIVANT DE LA LIGNE
         MNSL = MNSL + 3
 50   CONTINUE
      RETURN
      END
