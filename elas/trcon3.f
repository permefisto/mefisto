       SUBROUTINE TRCON3( NOFORE, NBELFI, NBPIEF, NDIMES, NBCAS,
     %                    COPIEF, CONPRI, DIRPRI, XYZPOI,
     %                    NCAS,   CMPCON, KNOMOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES CONTRAINTES PRINCIPALES D'UN ENSEMBLE D'EF 3D
C -----                   VERSION xvue
C
C ENTREES :
C ---------
C NOFORE : 0 SI PAS DE FONCTION UTILISATEUR 'REGION', >0 SINON
C NBELFI : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NBPIEF : NOMBRE DE POINTS DE CALCUL DES CONTRAINTES PAR ELEMENT
C NDIMES : ESPACE DE TRAVAIL 2 OU 3
C NBCAS  : NOMBRE DE CAS TRAITES
C COPIEF : LES NDIMES COORDONNEES DES POINTS DE CALCUL DES CONTRAINTES
C CONPRI : LES NDIMES CONTRAINTES PRINCIPALES
C DIRPRI : LES NDIMES DIRECTIONS PRINCIPALES
C XYZPOI : LES 3 COORDONNEES DES POINTS DE L'OBJET
C NCAS   : LE NUMERO DU CAS A TRACER
C LATOPO : LE TABLEAU TOPOLOGIE
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C NCAS   : NUMERO DU CAS A TRAITER
C CMPCON : NOMBRE DE CM POUR L'UNITE DE CONTRAINTE
C KNOMOB : NOM DE L'OBJET A TRAITER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1998
C23456---------------------------------------------------------------012
      PARAMETER        (LIGCON=0, LIGTIR=1 )
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/ctemps.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___face.inc"
      CHARACTER*(*)     KNOMOB
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      REAL              CXYZCM(3)
      EQUIVALENCE      (CXYZCM(1),CXCM), (CXYZCM(2),CYCM),
     %                 (CXYZCM(3),CZCM)
      REAL               XYZ(3)
      EQUIVALENCE      (XYZ(1),X),(XYZ(2),Y),(XYZ(3),Z)
      REAL              COPIEF(1:NBELFI,1:NBPIEF,1:NDIMES)
      REAL              XYZPOI(3,*)
      DOUBLE PRECISION  CONPRI(1:NBELFI,1:NBPIEF,1:NDIMES,1:NBCAS)
      REAL     DIRPRI(1:NBELFI,1:NBPIEF,1:NDIMES,1:NDIMES,1:NBCAS)
      DOUBLE PRECISION  VMAX,A(4)
      INTRINSIC         REAL
C
C     CREATION OU REDECOUVERTE DU TMS  OBJET>>>FACE
C     ---------------------------------------------
      CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
C
C     CREATION OU REDECOUVERTE DU TMS  OBJET>>>ARETEFR
C     DES ARETES DES FACES FRONTALIERES DE L'OBJET
C     ------------------------------------------------
      CALL HACHAF( KNOMOB, 0, NTFAOB, MNFAOB,
     %             NTAFOB, MNAFOB, I )
C
C     LE TRACE DES ARETES FRONTALIERES
C     --------------------------------
C     LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
      MOARFR = MCN( MNAFOB + WOARFR )
C     LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
      L1ARFR = MCN( MNAFOB + W1ARFR )
C     LIGNES NON EPAISSIES
      CALL XVEPAISSEUR( 1 )
C     LIGNES TIRETEES
      CALL XVTYPETRAIT( NTLAFR )
      CALL TRARFR( NCOUAF, MOARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %             XYZPOI )
C
C     TRACE DES CONTRAINTES PRINCIPALES DE L'OBJET 3D
C     -----------------------------------------------
C     TRACE AVEC 2 EPAISSEURS
      CALL XVEPAISSEUR( 2 )
C     TRACE AVEC LIGNE CONTINUE
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
C          COORDONNEES DU POINT CENTRE DES CONTRAINTES A TRACER
          XYZ(1)=COPIEF(K,L,1)
            XYZ(2)=COPIEF(K,L,2)
            XYZ(3)=COPIEF(K,L,3)
C
C           LA PREMIERE CONTRAINTE PRINCIPALE
            VP1 = REAL( CONPRI(K,L,1,NCAS) )
C           SA DIRECTION PRINCIPALE
          XV1=DIRPRI(K,L,1,1,NCAS)
            YV1=DIRPRI(K,L,2,1,NCAS)
            ZV1=DIRPRI(K,L,3,1,NCAS)
C          LONGUEUR CM DE LA FLECHE DE VP1
            CXCM  = VP1*XV1*CMPCON
            CYCM  = VP1*YV1*CMPCON
            CZCM  = VP1*ZV1*CMPCON
            CONCM = SQRT( CXCM**2 + CYCM**2 + CZCM**2 )
C          LES 3 FLECHES DE LA 1-ERE CONTRAINTE PRINCIPALE
C          LE SIGNE DE CONCM DONNE LE SENS DE LA FLECHE
          IF (VP1 .LT. 0.D0 ) CONCM=-CONCM
            CALL T3FLEC( NCOUFL, XYZ, CONCM, CXYZCM )
C           LA FLECHE OPPOSEE
            CXCM=-CXCM
            CYCM=-CYCM
            CZCM=-CZCM
            CALL T3FLEC( NCOUFL, XYZ, CONCM, CXYZCM )
C
C           LA DEUXIEME CONTRAINTE PRINCIPALE
            VP2 = REAL( CONPRI(K,L,2,NCAS) )
C           SA DIRECTION PRINCIPLALE
            XV1=DIRPRI(K,L,1,2,NCAS)
            YV1=DIRPRI(K,L,2,2,NCAS)
            ZV1=DIRPRI(K,L,3,2,NCAS)
C           LONGUEUR CM DE LA FLECHE DE VP2
            CXCM  = VP2*XV1*CMPCON
            CYCM  = VP2*YV1*CMPCON
            CZCM  = VP2*ZV1*CMPCON
            CONCM = SQRT( CXCM**2 + CYCM**2 + CZCM**2 )
C          LES 3 FLECHES DE LA 2-ERE CONTRAINTE PRINCIPALE
C          LE SIGNE DE CONCM DONNE LE SENS DE LA FLECHE
          IF (VP2 .LT. 0.D0 ) CONCM=-CONCM
            CALL T3FLEC( NCOUFL, XYZ , CONCM , CXYZCM )
C           LA FLECHE OPPOSEE
            CXCM=-CXCM
            CYCM=-CYCM
            CZCM=-CZCM
            CALL T3FLEC( NCOUFL, XYZ, CONCM, CXYZCM )
C
C           LA TROISIEME CONTRAINTE PRINCIPALE
            VP3 = REAL( CONPRI(K,L,3,NCAS) )
C           SA DIRECTION PRINCIPLALE
            XV1=DIRPRI(K,L,1,3,NCAS)
            YV1=DIRPRI(K,L,2,3,NCAS)
            ZV1=DIRPRI(K,L,3,3,NCAS)
C           LONGUEUR CM DE LA FLECHE DE VP3
            CXCM  = VP3*XV1*CMPCON
            CYCM  = VP3*YV1*CMPCON
            CZCM  = VP3*ZV1*CMPCON
            CONCM = SQRT( CXCM**2 + CYCM**2 + CZCM**2 )
C          LES 3 FLECHES DE LA 3-ERE CONTRAINTE PRINCIPALE
C          LE SIGNE DE CONCM DONNE LE SENS DE LA FLECHE
          IF (VP3 .LT. 0.D0 ) CONCM=-CONCM
            CALL T3FLEC( NCOUFL, XYZ , CONCM , CXYZCM )
C           LA FLECHE OPPOSEE
            CXCM=-CXCM
            CYCM=-CYCM
            CZCM=-CZCM
            CALL T3FLEC( NCOUFL, XYZ, CONCM, CXYZCM )
C
 90     CONTINUE
 100  CONTINUE
C
C     RETOUR A 1 EPAISSEUR DE TRACE DES LIGNES
      CALL XVEPAISSEUR( 1 )
      CALL XVCOULEUR( NCBLAN )
C
      RETURN
      END
