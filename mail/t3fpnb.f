      SUBROUTINE T3FPNB( NMOBJT , NUOBJT ,
     %                   NBMOFA , NBFACE , LFACES , XYZSOM ,
     %                   NBFAFR , NOFAFR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN NOIR ET BLANC LES FACES DU CONTOUR D'UN VOLUME
C -----    UNE FACE EST TRACEE SI ELLE APPARTIENT A UN SEUL ELEMENT
C          ET SI ELLE EST ECLAIREE PAR LA DIRECTION DE VISEE
C
C ENTREES:
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DU VOLUME
C NBMOFA : NOMBRE DE MOTS PAR FACE
C NBFACE : NOMBRE DE FACES
C LFACES : TABLEAU ENTIER DU NO DES SOMMETS DES FACES DU VOLUME
C XYZSOM : COORDONNEES DES SOMMETS
C NBFAFR : NOMBRE DE FACES FRONTALIERES
C NOFAFR : NUMERO DES FACES SELON LEUR DISTANCE CROISSANTE A L'OEIL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1991
C ...................................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      INTEGER           LFACES (NBMOFA,NBFACE),NOFAFR(1:NBFAFR)
      CHARACTER*(*)     NMOBJT
      CHARACTER*8       NMSOMM
      REAL              XYZSOM(3,*)
      REAL              XYZ(3,4),XYZP(3)
C
C     TRACE NOIR ET BLANC
C     REDUCTION DES FACES
      REDUCF = PREDUF * 0.01
      REDUC1 = 1.0 - REDUCF
C
C     LE TRACE DES FACES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     =======================================================
      DO 200 NF = 1, NBFAFR
C
C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NF1 = NOFAFR( NF )
C
C        LE NOMBRE DE SOMMETS DE LA FACE
         IF( LFACES(4,NF1) .GT. 0 ) THEN
            NBSF = 4
         ELSE
            NBSF = 3
         ENDIF
C
         DO 110 J=1,NBSF
C           LES COORDONNEES DU SOMMET J
            NSOM = LFACES(J,NF1)
            XYZ(1,J) = XYZSOM(1,NSOM)
            XYZ(2,J) = XYZSOM(2,NSOM)
            XYZ(3,J) = XYZSOM(3,NSOM)
 110     CONTINUE
C
         IF( IAVNEF .NE. 0 .OR. REDUCF .GT. 0.0 ) THEN
C
C           REDUCTION DES FACES
C           -------------------
C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
            CALL COBAPO( NBSF, XYZ, XYZP )
C           L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
            DO 160 J=1,NBSF
               XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
               XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
               XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
 160        CONTINUE
         ENDIF
C
C        LE TRACE DE LA FACE EN NOIR
         CALL FACE3D( NCNOIR, NCBLAN, NBSF, XYZ )
C
C        TRACE EVENTUEL DU NO DE L'EF
         IF( IAVNEF .NE. 0 ) THEN
            IF( REDUCF .EQ. 0.0 ) THEN
C              CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
               CALL COBAPO( NBSF, XYZ, XYZP )
            ENDIF
            NE = ABS( LFACES(6,NF1) )
            WRITE( NMSOMM , '(I8)' ) NE
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, XYZP, NMSOMM(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            DO 180 J=1,NBSF
               NSOM = LFACES(J,NF1)
               WRITE( NMSOMM , '(I8)' ) NSOM
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, XYZSOM(1,NSOM), NMSOMM(1:L) )
 180        CONTINUE
         ENDIF
 200  CONTINUE
C
C     LE TRACE DE LA POIGNEE DU VOLUME
C     CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
      CALL COBAPO( NBSF, XYZ, XYZP )
      CALL ITEMV3( XYZP , NMOBJT , NUOBJT )
      END
