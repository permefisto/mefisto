      SUBROUTINE TG2ARL( NI, XYZ1, XYZ2, MNXYTG, MNNUTG, TG2AR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      CALCUL DES 3 COMPOSANTES DES 2 TANGENTES CROISEES
C -----      DE L'ARETE NI DE LA LIGNE
C
C ENTREES:
C --------
C NI     : LE  NUMERO DE L'ARETE A TRAITER
C XYZ1   : LES 3 COORDONNEES DU PREMIER SOMMET DE L'ARETE NI DE LA LIGNE
C XYZ2   : LES 3 COORDONNEES DU SECOND  SOMMET DE L'ARETE NI DE LA LIGNE
C MNXYTG : L'ADRESSE MCN DES 3 COMPOSANTES DE TOUTES LES TANGENTES
C          DE LA LIGNE
C MNNUTG : L'ADRESSE MCN DES 2 NUMEROS DES TANGENTES DE CHAQUE ARETE
C          DE LA LIGNE
C
C SORTIES:
C --------
C TG2AR  : LES 3 COMPOSANTES DES 2 TANGENTES CROISEES DE L'ARETE NI
C          DE LA LIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      REAL              XYZ1(3), XYZ2(3), TG2AR(3,2)
C
      IF( MNXYTG .LE. 0 .OR. MNNUTG .LE. 0 ) THEN
C        PAS DE TANGENTES: LES TANGENTES SONT LES ARETES
         TG2AR(1,1) = XYZ2(1) - XYZ1(1)
         TG2AR(1,2) = -TG2AR(1,1)
         TG2AR(2,1) = XYZ2(2) - XYZ1(2)
         TG2AR(2,2) = -TG2AR(2,1)
         TG2AR(3,1) = XYZ2(3) - XYZ1(3)
         TG2AR(3,2) = -TG2AR(3,1)
         RETURN
      ENDIF
C
C     LES TABLEAUX DES TANGENTES EXISTENT
      DO 10 K=1,2
C
C        +- LE NUMERO DE LA TANGENTE K DE L'ARETE NI
         NTG = MCN(MNNUTG-3+K+2*NI)
C
         IF( NTG .NE. 0 ) THEN
C
C           IL EXISTE UNE TANGENTE AU SOMMET K DE L'ARETE NI
            IF( NTG .GT. 0 ) THEN
C              SENS A CONSERVER
               MN = MNXYTG + 3 * NTG - 3
               TG2AR(1,K) =  RMCN(MN  )
               TG2AR(2,K) =  RMCN(MN+1)
               TG2AR(3,K) =  RMCN(MN+2)
            ELSE
C              SENS A INVERSER
               MN = MNXYTG - 3 * NTG - 3
               TG2AR(1,K) = -RMCN(MN  )
               TG2AR(2,K) = -RMCN(MN+1)
               TG2AR(3,K) = -RMCN(MN+2)
            ENDIF
C
         ELSE
C
C           PAS DE TANGENTE: COTE DROIT
            IF( K .EQ. 1 ) THEN
               LSIGNE = 1
            ELSE
               LSIGNE = -1
            ENDIF
            TG2AR(1,K) = LSIGNE * ( XYZ2(1) - XYZ1(1) )
            TG2AR(2,K) = LSIGNE * ( XYZ2(2) - XYZ1(2) )
            TG2AR(3,K) = LSIGNE * ( XYZ2(3) - XYZ1(3) )
C
         ENDIF
 10   CONTINUE
      END
