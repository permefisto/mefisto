      SUBROUTINE TRAREF( NDIMES, NBPOI, XYZPOI,  NBELFI, NUPGEL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES D'UN TYPE D'ELEMENTS FINIS
C -----
C
C ENTREES :
C ---------
C NDIMES : ESPACE DE TRAVAIL 1 ou 2 ou 3
C NBPOI  : NOMBRE DE POINTS DE XYZPOI
C XYZPOI : LES 3 COORDONNEES DES POINTS DE L'OBJET
C NBELFI : NOMBRE D'ELEMENTS FINIS DU TYPE A TRAITER
C NUPGEL : LE NUMERO DES POINTS DE CHAQUE ELEMENT FINI DE CE TYPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      REAL           XYZPOI(3,NBPOI)
      INTEGER        NUPGEL(NBELFI,*)
C
      IF( NDIMES .EQ. 1 ) THEN
C
C        TRACE 1D DE L'ARETE DES EF
C        --------------------------
         RMX = XYZPOI(1,1)
         RMI = XYZPOI(1,1)
         DO 5 K=2,NBPOI
            IF( XYZPOI(1,K) .LT. RMI ) RMI = XYZPOI(1,K)
            IF( XYZPOI(1,K) .GT. RMX ) RMX = XYZPOI(1,K)
 5       CONTINUE
         Y = XYZPOI(2,1) - (RMX-RMI) / 20
C
         DO 10 K=1,NBELFI
C
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NS1 = NUPGEL( K, NOSOAR(1,1) )
            NS2 = NUPGEL( K, NOSOAR(2,1) )
C
C           TRACE DE L'ARETE DROITE  NS1-NS2
            CALL TRAIT2D( NCOUAF, XYZPOI(1,NS1), XYZPOI(2,NS1),
     %                            XYZPOI(1,NS2), XYZPOI(2,NS2) )
            CALL TRAIT2D( NCOUAF, XYZPOI(1,NS1), Y,
     %                            XYZPOI(1,NS2), Y )
C
C           TRACE DU SYMBOLE DES 2 SOMMETS
            CALL SYMBOLE2D( NCGRIS, XYZPOI(1,NS1), XYZPOI(2,NS1), 'I' )
            CALL SYMBOLE2D( NCGRIS, XYZPOI(1,NS2), XYZPOI(2,NS2), 'I' )
            CALL SYMBOLE2D( NCGRIS, XYZPOI(1,NS1), Y,             'I' )
            CALL SYMBOLE2D( NCGRIS, XYZPOI(1,NS2), Y,             'I' )
C
 10      CONTINUE
C
      ELSE IF( NDIMES .EQ. 2 ) THEN
C
C        TRACE 2D DES ARETES DES EF
C        --------------------------
         DO 30 K=1,NBELFI
            DO 20 I=1,NARET
C
C              LE NUMERO DES 2 SOMMETS DE L'ARETE
               NS1 = NUPGEL( K, NOSOAR(1,I) )
               NS2 = NUPGEL( K, NOSOAR(2,I) )
C
               IF( NBNOAR(I) .LE. 0 ) THEN
C
C                 TRACE DE L'ARETE DROITE  NS1-NS2
                  CALL TRAIT2D( NCOUAF, XYZPOI(1,NS1), XYZPOI(2,NS1),
     %                                  XYZPOI(1,NS2), XYZPOI(2,NS2) )
               ELSE
C
C                 LE NUMERO DU MILIEU DE L'ARETE
                  NM = NUPGEL( K, NOPOAR(1,I) )
C
C                 TRACE DU SOMMET1 AU MILIEU
                  CALL TRAIT2D( NCOUAF, XYZPOI(1,NS1), XYZPOI(2,NS1),
     %                                  XYZPOI(1,NM ), XYZPOI(2,NM ) )
C
C                 TRACE DU MILIEU AU SOMMET 2
                  CALL TRAIT2D( NCOUAF, XYZPOI(1,NM ), XYZPOI(2,NM ),
     %                                  XYZPOI(1,NS2), XYZPOI(2,NS2) )
               ENDIF
 20         CONTINUE
 30      CONTINUE
C
      ELSE IF( NDIMES .EQ. 3 ) THEN
C
C        TRACE 3D DES ARETES DES EF
C        --------------------------
         DO 90 K=1,NBELFI
            DO 80 I=1,NARET
C
C              LE NUMERO DES 2 SOMMETS DE L'ARETE
               NS1 = NUPGEL( K, NOSOAR(1,I) )
               NS2 = NUPGEL( K, NOSOAR(2,I) )
C
               IF( NBNOAR(I) .LE. 0 ) THEN
C
C                 TRACE DE L'ARETE DROITE  NS1-NS2
                  CALL TRAIT3D( NCOUAF, XYZPOI(1,NS1), XYZPOI(1,NS2) )
C
               ELSE
C
C                 LE NUMERO DU MILIEU DE L'ARETE
                  NM = NUPGEL( K, NOPOAR(1,I) )
C
C                 TRACE DU SOMMET 1 AU MILIEU
                  CALL TRAIT3D( NCOUAF, XYZPOI(1,NS1), XYZPOI(1,NM) )
C
C                 TRACE DU MILIEU AU SOMMET 2
                  CALL TRAIT3D( NCOUAF, XYZPOI(1,NM ), XYZPOI(1,NS2) )
               ENDIF
 80         CONTINUE
 90      CONTINUE
      ENDIF
      END
