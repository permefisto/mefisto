      SUBROUTINE CHARIS( NBISO, LTARCH, NSARCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CHAINER PAR CONTINUITE LES ARETES DE CHAQUE ISOVALEUR
C ----- SI CONTINUITE EFFECTIVE NSARCH(2,NAR)=NSARCH(1,NSARCH(4,NAR))
C       SINON L'ISOVALEUR N'EST PAS CONTINUE EN CES 2 POINTS
C
C ENTREES:
C --------
C NBISO  : NOMBRE D'ISOVALEURS
C
C MODIFIES:
C ---------
C LTARCH : LA TETE DE LISTE DE CHAQUE ISOVALEUR (POINTEUR SUR LA 1ERE ARETE)
C NSARCH : NUMERO DES 2 SOMMETS DE L'ARETE, ARETE PRECEDENTE ET SUIVANTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1995
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER  LTARCH(NBISO), NSARCH(4,*)
C
C     BOUCLE SUR LES ISOVALEURS
      DO 100 ISO = 1 , NBISO
C
C        LA TETE DE LA LISTE DES ARETES DE L'ISO
         NA0 = LTARCH( ISO )
         IF( NA0 .LE. 0 ) GOTO 100
C
C        CONSTRUCTION DE LA COMPOSANTE CONNEXE DE L'ARETE NA0
         NBC = 0
 1       NA  = NA0
C
C        LE SECOND SOMMET DE L'ARETE
 5       NS2 = NSARCH(2,NA)
         NAS = NA
C
C        RECHERCHE PARMI LES ARETES CHAINEES DE L'ARETE SUIVANTE NS2...
 10      NAS = NSARCH(4,NAS)
         IF( NAS .GT. 0 ) THEN
C
            IF( NSARCH(2,NAS) .EQ. NS2 ) THEN
C              SOMMET RETROUVE. LE SENS DE L'ARETE EST INVERSE
               N             = NSARCH(1,NAS)
               NSARCH(1,NAS) = NSARCH(2,NAS)
               NSARCH(2,NAS) = N
            ENDIF
            IF( NSARCH(1,NAS) .NE. NS2 ) GOTO 10
C
C           SOMMET RETROUVE EN POSITION 1 DE L'ARETE
C           CHAINAGE AVEC L'ARETE NA
            IF( NSARCH(3,NAS) .NE. NA ) THEN
C
C              NAS N'EST PAS L'ARETE SUIVANTE DE NA
C              NAS EST RETIREE DU CHAINAGE DES ARETES
               MAP = NSARCH(3,NAS)
               MAS = NSARCH(4,NAS)
               IF( MAP .GT. 0 ) NSARCH(4,MAP) = MAS
               IF( MAS .GT. 0 ) NSARCH(3,MAS) = MAP
C
C              L'ARETE SUIVANTE DE NA DEVIENT L'ARETE SUIVANTE DE NAS
               MAS = NSARCH(4,NA)
               NSARCH(4,NAS) = MAS
               NSARCH(3,MAS) = NAS
C              NAS EST CHAINEE COMME ARETE SUIVANTE DE NA
               NSARCH(3,NAS) = NA
               NSARCH(4,NA ) = NAS
            ENDIF
C
C           ICI NA EST L'ARETE PRECEDENTE DE NAS. PASSAGE A L'ARETE SUIVANTE
            NA = NAS
            GOTO 5
         ENDIF
C
C        PLUS D'ARETE DE SOMMET NS2
C        NAF EST LA DERNIERE ARETE EN CONTINU
         NAF = NA
C        EXISTE-T-IL ENCORE UNE ARETE NON CHAINEE?
         NAS = NSARCH(4,NAF)
         IF( NAS .GT. 0 ) THEN
C
C           IL RESTE AU MOINS UNE ARETE
C           RECHERCHE DU SOMMET 1 DE LA TETE D'ISO
            NA = NA0
C
C           LE PREMIER SOMMET DE L'ARETE
 25         NS1 = NSARCH(1,NA)
            GOTO 32
C
C           RECHERCHE PARMI LES ARETES CHAINEES DE L'ARETE ...NS1
 30         NAS = NSARCH(4,NAS)
 32         IF( NAS .GT. 0 ) THEN
C
               IF( NSARCH(1,NAS) .EQ. NS1 ) THEN
C                 SOMMET RETROUVE. LE SENS DE L'ARETE EST INVERSE
                  N             = NSARCH(1,NAS)
                  NSARCH(1,NAS) = NSARCH(2,NAS)
                  NSARCH(2,NAS) = N
               ENDIF
               IF( NSARCH(2,NAS) .NE. NS1 ) GOTO 30
C
C              SOMMET RETROUVE EN POSITION 2 DE L'ARETE NAS
C              CHAINAGE PRECEDENT AVEC L'ARETE NA
               MAS = NSARCH(4,NAS)
               IF( MAS .NE. NA ) THEN
C
C                 NA N'EST PAS L'ARETE SUIVANTE DE NAS
C                 NAS EST RETIREE DU CHAINAGE DES ARETES
                  MAP = NSARCH(3,NAS)
                  IF( MAP .GT. 0 ) NSARCH(4,MAP) = MAS
                  IF( MAS .GT. 0 ) NSARCH(3,MAS) = MAP
C
C                 LA PRECEDENTE DE NA EST CHAINEE EN PRECEDENTE DE NAS
                  MAP = NSARCH(3,NA)
                  NSARCH(3,NAS) = MAP
                  IF( MAP .GT. 0 ) NSARCH(4,MAP) = NAS
C                 NAS EST CHAINEE COMME ARETE PRECEDENTE DE NA
                  NSARCH(4,NAS) = NA
                  NSARCH(3,NA ) = NAS
               ENDIF
C
C              ICI NA EST L'ARETE PRECEDENTE DE NAS. PASSAGE A L'ARETE SUIVANTE
               NA  = NAS
               NAS = MAS
               GOTO 25
            ENDIF
C
C           PLUS D'ARETE DE SOMMET NS1
C           NAD EST LA PREMIERE ARETE EN CONTINU DE CETTE COMPOSANTE
            NAD = NA
C
            IF( NBC .EQ. 0 ) THEN
C              MISE A JOUR DE LA TETE DE L'ISO
               LTARCH(ISO) = NAD
            ENDIF
            NBC = NBC + 1
         ENDIF
C
C        LA COMPOSANTE CONNEXE DE NA0 VA DE L'ARETE NAD A L'ARETE NAF
         NA0 = NSARCH(4,NAF)
         IF( NA0 .GT. 0 ) THEN
C           IL RESTE DES ARETES => TRAITONS UNE AUTRE COMPOSANTE CONNEXE
            GOTO 1
         ENDIF
         WRITE(IMPRIM,*) 'ISOVALEUR ',ISO,' AVEC ',NBC,' COMPOSANTES'
 100  CONTINUE
      END
