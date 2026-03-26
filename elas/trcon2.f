      SUBROUTINE TRCON2( NOFORE , NBELFI , NBPIEF , NDIMES , NBCAS ,
     %                   COPIEF , CONPRI , DIRPRI ,
     %                   XYZPOI , NUPGEL ,
     %                   NCAS   , CMPCON )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES CONTRAINTES PRINCIPALES D'UN ENSEMBLE D'EF en 2D
C -----                   VERSION xvue
C
C ENTREES :
C ---------
C NOFORE : 0 SI PAS DE FONCTION UTILISATEUR 'REGION', >0 SINON
C NBELFI : NOMBRE D'ELEMENTS
C NBPIEF : NOMBRE DE POINTS DE CALCUL DES CONTRAINTES PAR ELEMENT
C NDIMES : =2 DIMENSION DE L'ESPACE DE TRAVAIL
C NBCAS  : NOMBRE DE CAS TRAITES
C COPIEF : LES NDIMES COORDONNEES DES POINTS DE CALCUL DES CONTRAINTES
C CONPRI : LES NDIMES CONTRAINTES PRINCIPALES
C DIRPRI : LES NDIMES DIRECTIONS PRINCIPALES
C XYZPOI : LES 3 COORDONNEES DES POINTS DE L'OBJET
C NUPGEL : LE NUMERO DES POINTS DE CHAQUE ELEMENT
C NCAS   : LE NUMERO DU CAS A TRACER
C CMPCON : LE NOMBRE DE CM POUR L'UNITE DE CONTRAINTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1994
C23456---------------------------------------------------------------012
      PARAMETER        (LIGCON=0, LIGTIR=1 )
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/ctemps.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      REAL              COPIEF(1:NBELFI,1:NBPIEF,1:NDIMES)
      REAL              XYZPOI(3,*)
      INTEGER           NUPGEL(NBELFI,*)
      DOUBLE PRECISION  CONPRI(1:NBELFI,1:NBPIEF,1:NDIMES,1:NBCAS)
      REAL     DIRPRI(1:NBELFI,1:NBPIEF,1:NDIMES,1:NDIMES,1:NBCAS)
      DOUBLE PRECISION  VMAX,A(4)
      INTRINSIC         REAL
C
C     TRACE DES ARETES DES EF
C     ------------------------
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( NTLAFR )
C
      DO 30 K=1,NBELFI
         DO 20 I=1,NARET
            NS1 = NUPGEL( K , NOSOAR(1,I) )
            NS2 = NUPGEL( K , NOSOAR(2,I) )
            CALL TRAIT2D( NCOUAF, XYZPOI(1,NS1) , XYZPOI(2,NS1),
     %                            XYZPOI(1,NS2) , XYZPOI(2,NS2) )
 20      CONTINUE
 30   CONTINUE
C
C     TRACE DES CONTRAINTES
C     ---------------------
      CALL XVEPAISSEUR( 2 )
      CALL XVTYPETRAIT( LIGCON )
C
      DO 100 K=1,NBELFI
         DO 90 L=1,NBPIEF
C
C           IMPRESSION REDUITE DES CONTRAINTES
            IF( NOFORE .GT. 0 ) THEN
C
C              LES 4 PARAMETRES D'APPEL DE LA FONCTION 'CHOIX'
C              LE TEMPS EN 1-ER PARAMETRE
               A(1) = TEMPS
C              PUIS LES 3 COORDONNEES X Y Z DU NOEUD
               A(2) = COPIEF(K,L,1)
               A(3) = COPIEF(K,L,2)
               A(4) = COPIEF(K,L,3)
C              FONCTION CHOIX(TEMPS,X,Y,Z)
               CALL FONVAL( NOFORE, 4, A, NCODEV, VMAX )
               IF( NINT(VMAX) .NE. 0 ) THEN
                  J=1
               ELSE
                  J=0
               ENDIF
            ELSE
               J=1
            ENDIF
            IF( J .EQ. 0 ) GOTO 90
C
C           COORDONNEES DU POINT CENTRE DES CONTRAINTES A TRACER
            XF = COPIEF(K,L,1)
            YF = COPIEF(K,L,2)
C
C           LA PREMIERE CONTRAINTE PRINCIPALE
            VP1 = REAL( CONPRI(K,L,1,NCAS) )
C           SA DIRECTION PRINCIPALE
            XV1 = DIRPRI(K,L,1,1,NCAS)
            YV1 = DIRPRI(K,L,2,1,NCAS)
C
C           LONGUEUR CM DE LA FLECHE DE VP1 ( >VP2 )
            COXF  = VP1 * XV1 * CMPCON
            COYF  = VP1 * YV1 * CMPCON
            CONCM = SQRT( COXF * COXF + COYF * COYF )
C
C           LES 2 FLECHES DE LA 1-ERE CONTRAINTE PRINCIPALE
C           LE SIGNE DE CONCM DONNE LE SENS DE LA FLECHE
            IF( VP1 .LT. 0.D0 ) CONCM = - CONCM
            CALL T2FLEC( NCOUFL, XF,YF,CONCM, COXF, COYF )
            CALL T2FLEC( NCOUFL, XF,YF,CONCM,-COXF,-COYF )
C
C           LA SECONDE CONTRAINTE PRINCIPALE
            VP2   = REAL( CONPRI(K,L,2,NCAS) )
C           SA DIRECTION PRINCIPALE
            XV1   = DIRPRI(K,L,1,2,NCAS)
            YV1   = DIRPRI(K,L,2,2,NCAS)
C           LA LONGUEUR DE LA FLECHE DE VP2
            COXF  = VP2 * XV1 * CMPCON
            COYF  = VP2 * YV1 * CMPCON
            CONCM = SQRT( COXF * COXF + COYF * COYF )
C
C           LES 2 FLECHES DE LA 2-EME CONTRAINTE PRINCIPALE
C           LE SIGNE DE CONCM DONNE LE SENS DE LA FLECHE
            IF( VP2 .LT. 0.D0 ) CONCM = - CONCM
            CALL T2FLEC( NCOUFL, XF,YF,CONCM, COXF, COYF )
            CALL T2FLEC( NCOUFL, XF,YF,CONCM,-COXF,-COYF )
 90     CONTINUE
 100  CONTINUE
C
      CALL XVEPAISSEUR( 1 )
      CALL XVCOULEUR( NCBLAN )
      END
