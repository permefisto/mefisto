      SUBROUTINE TRATRINT( NBPT, NONOPT, XYZPT, NBTRA0, NBTRA1, NUSTRA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LES SOUS-TRIANGLES NBTRA0 a NBTRA1
C -----
C ENTREES:
C --------
C NBPT   : NOMBRE DE POINTS
C NONOPT : NOUVEAU NUMERO DE CHAQUE POINT (-NS NOUVEAU SI SUPPRIME)
C XYZPT  : 3 XYZ DES POINTS-SOMMETS
C NBTRA0 : NUMERO DU PREMIER TRIANGLE A TRACER
C NBTRA1 : NUMERO DU DERNIER TRIANGLE A TRACER
C NUSTRA : (3,NBTRA) NO DANS XYZPT DES 3 SOMMETS DES TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray Decembre 2011
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
C
      DOUBLE PRECISION  XYZPT(3,*), QUALIT
      INTEGER           NUSTRA(3,*), NONOPT(*)
C
      REAL              XYZB(3)
      REAL              XYZ1(3), XYZ2(3), XYZ3(3), XYZTR(3,3)
      EQUIVALENCE      (XYZ1(1),XYZTR(1,1)),
     %                 (XYZ2(1),XYZTR(1,2)),
     %                 (XYZ3(1),XYZTR(1,3))
      INTEGER           NOSOTR(3)
      EQUIVALENCE      (NOSOTR(1),NS1),(NOSOTR(2),NS2),(NOSOTR(3),NS3)
C
C     LE CADRE MIN-MAX
      DO K=1,3
         COOEXT(K,1) = REAL( XYZPT(K,1) )
         COOEXT(K,2) = REAL( XYZPT(K,1) )
      ENDDO
      DO M=1,NBPT
         DO K=1,3
            S = REAL( XYZPT(K,M) )
            IF( S .LT. COOEXT(K,1) ) THEN
               COOEXT(K,1) = S
            ENDIF
            IF( S .GT. COOEXT(K,2) ) THEN
               COOEXT(K,2) = S
            ENDIF
         ENDDO
      ENDDO
C
C     LES TRACES SONT DEMANDES
      LORBITE = 1
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9900
      ENDIF
C
C     LA MEMOIRE PIXEL EST EFFACEE
 10   CALL EFFACEMEMPX
C
C     TRACE DES AXES
      CALL TRAXE3
C
C     LES SOUS-TRIANGLES
C     ------------------
      DO 100 NT=NBTRA0,NBTRA1
C
C        TRACE DU TRIANGLE NT
         IF( NUSTRA(1,NT) .LE. 0 ) GOTO 100
C
C        NO DES SOMMETS DU TRIANGLE NT
         NS1 = NOUVNOPT( NUSTRA(1,NT), NONOPT )
         NS2 = NOUVNOPT( NUSTRA(2,NT), NONOPT )
         NS3 = NOUVNOPT( NUSTRA(3,NT), NONOPT )
C
C        QUALITE DU TRIANGLE NT
         CALL QUATRID( NOSOTR, XYZPT, QUALIT )
C
C        LA COULEUR DE LA FACE VISUALISE LA QUALITE
         NCF = N1COUL + 9 - INT( 10 * ( 1 - QUALIT ) )
         NCF = MAX( NCF, N1COUL )
C
C        LES X,Y DES 3 SOMMETS POUR LE TRACE DU TRIANGLE NT
C        NUMERO DU SOMMET DANS CE TABLEAU
         XYZ1(1) = REAL( XYZPT(1,NS1) )
         XYZ1(2) = REAL( XYZPT(2,NS1) )
         XYZ1(3) = REAL( XYZPT(3,NS1) )
C
C        NUMERO DU SOMMET 2 DANS CE TABLEAU
         XYZ2(1) = REAL( XYZPT(1,NS2) )
         XYZ2(2) = REAL( XYZPT(2,NS2) )
         XYZ2(3) = REAL( XYZPT(3,NS2) )
C
C        NUMERO DU SOMMET 3 DANS CE TABLEAU
         XYZ3(1) = REAL( XYZPT(1,NS3) )
         XYZ3(2) = REAL( XYZPT(2,NS3) )
         XYZ3(3) = REAL( XYZPT(3,NS3) )
C
C        LE TRACE DE LA FACE ET DES 3 ARETES DU TRIANGLE NT
         CALL FACE3D( NCF, NCBLAN, 3, XYZTR )
C
C        TRACE DU NUMERO DES SOMMETS DU TRIANGLE NT
         CALL ENTIER3D( NCNOIR, XYZ1, NS1 )
         CALL ENTIER3D( NCNOIR, XYZ2, NS2 )
         CALL ENTIER3D( NCNOIR, XYZ3, NS3 )
C
C        LE BARYCENTRE DE NT
         XYZB(1) = (XYZ1(1) + XYZ2(1) + XYZ3(1) ) / 3
         XYZB(2) = (XYZ1(2) + XYZ2(2) + XYZ3(2) ) / 3
         XYZB(3) = (XYZ1(3) + XYZ2(3) + XYZ3(3) ) / 3
         CALL ENTIER3D( NCMAGE, XYZB, NT )
C
 100  CONTINUE
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     REPRISE DE L'ORBITE
      IF( LORBITE .GT. 0 ) THEN
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .NE. 0 ) GOTO 10
      ENDIF
C
 9900 RETURN
      END
