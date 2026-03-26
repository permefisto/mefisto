      SUBROUTINE T3TOUR( NMOBJT , NUOBJT ,
     %                   NBMOFA , NBFACE , LFACES , XYZSOM ,
     %                   NBFAFR , NOFAFR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER FIL DE FER LES FACES DU CONTOUR D UN VOLUME
C ----- UNE FACE EST TRACEE SI ELLE APPARTIENT A UN SEUL ELEMENT
C
C PARAMTRES DONNEES :
C -------------------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DU VOLUME
C NBMOFA : NOMBRE DE MOTS PAR FACE
C NBFACE : NOMBRE DE FACES
C LFACES : TABLEAU DU NO DES SOMMETES DES FACES DU MAILLAGE
C XYZSOM : COORDONNEES DES SOMMETS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C ...................................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      INTEGER           LFACES (NBMOFA,NBFACE),NOFAFR(1:NBFAFR)
      REAL              XYZSOM(3,*),XYZP(3),XYZ(3,4)
      EQUIVALENCE      (XYZP(1),XX),(XYZP(2),YY),(XYZP(3),ZZ)
      CHARACTER*(*)     NMOBJT
      CHARACTER*8       NMSOMM
C
C     DESCRIPTION DU VOLUME
C     =====================
C     REDUCTION DES FACES
      REDUCF = PREDUF * 0.01
      REDUC1 = 1.0 - REDUCF
C
C     LE TRACE DES FACES FRONTALIERES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     ====================================================================
      NBS = 0
      NF1 = 0
      DO 200 NF = 1, NBFAFR
C
C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NF1 = NOFAFR( NF )
C
C        LE NOMBRE DE SOMMETS DE LA FACE
         IF( LFACES(4,NF1) .GT. 0 ) THEN
            NBS = 4
         ELSE
            NBS = 3
         ENDIF
C
         IF( IAVNEF .NE. 0 .OR. REDUCF .GT. 0.0 ) THEN
C
C           REDUCTION DES FACES
C           -------------------
C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE
            XX = 0
            YY = 0
            ZZ = 0
            DO 110 J=1,NBS
               NS = LFACES(J,NF1)
               XYZ(1,J) = XYZSOM(1,NS)
               XX   = XX + XYZ(1,J)
               XYZ(2,J) = XYZSOM(2,NS)
               YY   = YY + XYZ(2,J)
               XYZ(3,J) = XYZSOM(3,NS)
               ZZ   = ZZ + XYZ(3,J)
 110        CONTINUE
            XYZP(1) = XX / NBS
            XYZP(2) = YY / NBS
            XYZP(3) = ZZ / NBS
            DO 120 J=1,NBS
               XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
               XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
               XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
 120        CONTINUE
C
         ELSE
C
            DO 130 J=1,NBS
               NS = LFACES(J,NF1)
               XYZ(1,J) = XYZSOM(1,NS)
               XYZ(2,J) = XYZSOM(2,NS)
               XYZ(3,J) = XYZSOM(3,NS)
 130        CONTINUE
         ENDIF
C
C        LE TRACE DES NBS ARETES AVEC LA COULEUR NCOUAF
         NS = NBS
         DO 140 J=1,NBS
            CALL TRAIT3D( NCOUAF, XYZ(1,NS), XYZ(1,J) )
            NS = J
 140     CONTINUE
C
C        TRACE EVENTUEL DU NO DE L'EF
         IF( IAVNEF .NE. 0 ) THEN
            NE = ABS( LFACES(6,NF1) )
            WRITE( NMSOMM , '(I8)' ) NE
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, XYZP, NMSOMM(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            DO 158 J=1,NBS
               NS = LFACES(J,NF1)
               WRITE( NMSOMM , '(I8)' ) NS
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, XYZ(1,J), NMSOMM(1:L) )
 158        CONTINUE
         ENDIF
C
 200  CONTINUE
C
C     LE TRACE DE LA POIGNEE DU VOLUME DANS LA DERNIERE FACE TRACEE
      IF( IAVNEF .EQ. 0 .AND. REDUCF .EQ. 0.0 ) THEN
         XX = 0.
         YY = 0.
         ZZ = 0.
         DO 310 J=1,NBS
            NS = LFACES(J,NF1)
            XX = XX + XYZSOM(1,NS)
            YY = YY + XYZSOM(2,NS)
            ZZ = ZZ + XYZSOM(3,NS)
 310     CONTINUE
         XX = XX / NBS
         YY = YY / NBS
         ZZ = ZZ / NBS
      ENDIF
      CALL ITEMV3( XYZP , NMOBJT , NUOBJT )
      END
