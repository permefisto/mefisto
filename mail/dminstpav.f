      SUBROUTINE DMINSTPAV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                      XYZ,    PTXYZD,
     %                      MXSMIN, NBSMIN, NOSMIN, NODSNO, DISTMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETROUVER LES PLUS PROCHES POINTS DE XYZ RANGES DANS LES PAVES
C -----
C ENTREES :
C ---------
C HEXAPAVE: MIN ET MAX DES COORDONNEES DU PAVAGE
C NBIPAV  : NOMBRE D'ARETES DANS LA DIRECTION I
C ECHPAV  : ECHELLE DANS LA DIRECTION I
C N1SPAVE : NO DU 1-ER SOMMET DANS PTXYZD DU PAVE
C NOPTSUIV: NO DU POINT SUIVANT DANS LE CHAINAGE DES POINTS DES PAVES
C XYZ     : X Y Z DU POINT A TRAITER
C PTXYZD  : X Y Z DISTANCE SOUHAITEE DES POINTS
C MXSMIN  : NOMBRE MAXIMAL DE SOMMETS VOISINS DE XYZ

C SORTIES :
C ---------
C NBSMIN  : NOMBRE DE SOMMETS PROCHES
C NOSMIN  : >0 NO DANS PTXYZD DES NBSMIN PLUS PROCHES POINT DE XYZ
C           =0 SI PAS DE POINT PROCHE DANS LE PAVAGE
C NODSNO  : NODSNO(1)=NO DANS NOSMIN DU PLUS PROCHE POINT DE PTXYZD
C           APRES TRI CROISSANT DES DISTANCES A XYZ
C DISTMIN : DISTANCE AU CARRE DU POINT DU PAVAGE LE PLUS PROCHE DE XYZ
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 2005
C2345X7..............................................................012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NOSMIN(MXSMIN), NODSNO(MXSMIN)
      DOUBLE PRECISION  DISTMIN(MXSMIN)
      INTEGER           NBIPAV(3), N1SPAVE(0:*), NOPTSUIV(*)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3)
      DOUBLE PRECISION  XYZ(3), PTXYZD(4,*)
      INTEGER           N1PAV(3), N2PAV(3)

      NBSMIN = 0
C     NOMBRE D'INTERVALLES POUR LA BOUCLE DE RECHERCHE DES SOMMETS PROCHES
      NBP    = 1

C     LE NUMERO EN X Y Z DES PAVES SUSCEPTIBLES DE CONTENIR CE POINT XYZ
 1    DO K=1,3
C        RECHERCHE DES INDICES DU PAVAGE EN XK POUVANT CONTENIR LE SOMMET
         M = INT( ( XYZ(K) - HEXAPAVE(K,1) ) * ECHPAV(K) )
         IF( M .LE. 0 ) THEN
            N1PAV(K) = 0
            N2PAV(K) = MIN( NBP, NBIPAV(K)-1 )
         ELSE IF( M .GE. (NBIPAV(K)-1) ) THEN
            N1PAV(K) = MAX( 0, NBIPAV(K)-1-NBP )
            N2PAV(K) = NBIPAV(K)-1
         ELSE
            N1PAV(K) = M - NBP
            N2PAV(K) = M + NBP
         ENDIF
C        PROTECTION DES PAVES NEGATIFS
         N1PAV(K) = MAX( 0, N1PAV(K) )
      ENDDO

C     PARCOURS DES PAVES POSSIBLES
      DO NZ = N1PAV(3), N2PAV(3)
         MZ = NZ * NBIPAV(2)
         DO NY = N1PAV(2), N2PAV(2)
            MY = ( NY + MZ ) * NBIPAV(1)
            DO NX = N1PAV(1), N2PAV(1)
C              LE NO DU PAVE >=0
               NUPAVE = NX + MY

C              LE PREMIER SOMMET DU PAVE NUPAVE
               NOSOMM = N1SPAVE( NUPAVE )

 10            IF( NOSOMM .GT. 0 ) THEN
C                 LE SOMMET NOSOMM EXISTE : DISTANCE A XYZ
                  IF( NBSMIN .LT. MXSMIN ) THEN
                      NBSMIN = NBSMIN + 1
                      NOSMIN( NBSMIN) = NOSOMM
                      NODSNO( NBSMIN) = NBSMIN
                      DISTMIN(NBSMIN) = (XYZ(1)- PTXYZD(1,NOSOMM))**2
     %                                + (XYZ(2)- PTXYZD(2,NOSOMM))**2
     %                                + (XYZ(3)- PTXYZD(3,NOSOMM))**2
                   ELSE
                      GOTO 50
                   ENDIF

C                 PASSAGE AU SOMMET SUIVANT
                  NOSOMM = NOPTSUIV( NOSOMM )
                  GOTO 10
C
               ENDIF
            ENDDO
         ENDDO
      ENDDO

 50   IF( NBSMIN .EQ. 0 ) THEN

C        PAS DE POINT RETROUVE. ON AUGMENTE LE CUBE DE RECHERCHE
C        AUGMENTATION DU CUBE DE PAVES DE RECHERCHE?
         IF( NBP .GE. 3 ) THEN
C           NON: 7 CUBES DE LARGE NE DOIVENT PAS ETRE DEPASSES
C           PAS DE PAVE CONTENANT LE POINT XYZ

ccc            PRINT*,'dminstpav: PAS de PAVE CONTENANT XYZ=',XYZ,
ccc     %             ' pour', 2*NBP+1,' PAVES de RECHERCHE en xyz'
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            PRINT 10000, ((HEXAPAVE(I,J),I=1,3),J=1,2)
ccc         ELSE
ccc            PRINT 20000, ((HEXAPAVE(I,J),I=1,3),J=1,2)
ccc         ENDIF

            GOTO 9999
         ENDIF

C        OUI: NBP2 PAVES DE RECHERCHE EN X Y Z
         NBP   = NBP + 1
         NBP2  = 2*NBP+1
         MIN2P = MIN( NBIPAV(1), NBIPAV(2), NBIPAV(3) )
         IF( NBP2 .LE. MIN2P ) THEN
ccc            PRINT*,'dminstpav: Nombre de PAVES de RECHERCHE=',NBP2,
ccc     %             ' en xyz pour TROUVER',XYZ
            GOTO 1
         ENDIF

      ENDIF

ccc10000 FORMAT(
ccc     % ' PAVAGE de l''HEXAEDRE: COORDONNEES MIN-MAX'/
ccc     % ' XMIN=',G15.7,'  YMIN=',G15.7,'  ZMIN=',G15.7 /
ccc     % ' XMAX=',G15.7,'  YMAX=',G15.7,'  ZMAX=',G15.7 )

ccc20000 FORMAT(
ccc     % ' HEXAHEDRON PAVEMENT: MIN-MAX COORDINATES'/
ccc     % ' XMIN=',G15.7,'  YMIN=',G15.7,'  ZMIN=',G15.7 /
ccc     % ' XMAX=',G15.7,'  YMAX=',G15.7,'  ZMAX=',G15.7 )

      IF( NBSMIN .GT. 1 ) THEN
C        TRI CROISSANT DES DISTANCES
         CALL TRITRD( NBSMIN, DISTMIN, NODSNO )
      ENDIF

ccc      IF( NBSMIN .GE. 1 .AND. NBP .GE. 3 ) THEN
ccc         PRINT*,'dminstpav: XYZ=',XYZ,' RETROUVE avec',
ccc     %           NBP2,' PAVES de RECHERCHE en xyz'
ccc      ENDIF

 9999 RETURN
      END
