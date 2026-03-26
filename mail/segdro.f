      SUBROUTINE SEGDRO( NTLXLI, XYZ1,   XYZ2,   NBARLI, RAIGEO,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREER LES ARETES D'UN SEGMENT DE DROITE D'EXTREMITES XYZ1 XYZ2
C -----
C
C ENTREES:
C --------
C NTLXLI : NUMERO DU TABLEAU TS DU LEXIQUE DE LA LIGNE
C XYZ1   : PREMIER POINT DU SEGMENT DE DROITE
C XYZ2   : DERNIER POINT DU SEGMENT DE DROITE
C NBARLI : NOMBRE D'ARETES DE LA LIGNE SI NOFOTI>0
C RAIGEO : RAISON DE LA PROGRESSION GEOMETRIQUE DES LONGUEURS DES ARETES
C
C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C          CF ~/td/d/a___nsef
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C          CF ~/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          1 SI NOMBRE DE SOMMETS INCORRECT <2 OU TAILLE_IDEALE INCALCULABLE
C          2 POINT INITIAL OU FINAL NON INITIALISE
C          3 POINT INITIAL ET FINAL CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1996
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/darete.inc"
      include"./incl/langue.inc"
C
C     LE SUPER-TABLEAU OU TOUS LES TMC ET TMS OUVERTS SONT STOCKES
C     ICI LE TABLEAU DMCN REEL DOUBLE PRECISION N'EST PAS UTILE
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE     ( MCN(1), RMCN(1) )
C
C     TABLEAUX AUXILIAIRES
      DOUBLE PRECISION  RAIGEO, H, PAS(3), XYZ(3), XYZ1(3), XYZ2(3),
     %                  XYZ0(3), XYZM(3), HHH
      DOUBLE PRECISION  DTAILL, DTAIL0, DTAIL1, DTAIL2, D
      INTRINSIC         NINT, SQRT
C
C     EN CAS D'ERREUR RENCONTREE
      NTXYZS = 0
      MNXYZS = 0
      NTNSEF = 0
      MNNSEF = 0
C
C     CALCUL DE LA DISTANCE ENTRE CES 2 POINTS
      H = SQRT( (XYZ2(1)-XYZ1(1))**2
     %        + (XYZ2(2)-XYZ1(2))**2
     %        + (XYZ2(3)-XYZ1(3))**2 )
      IF( H .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LA LIGNE EST REDUITE A UN POINT'
         ELSE
            KERR(1) = 'The LINE is REDUCED TO A POINT'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE(XYZ)' DES ARETES
C     ===============================================================
      IF( NOFOTI .GT. 0 ) THEN
C
C        LA FONCTION 'TAILLE_IDEALE' EXISTE
C        DECLARATION DU TABLEAU DES TAILLES DES ARETES
         MXTAIL = 512
         CALL TNMCDC( 'REEL', MXTAIL, MNTAIL )
C        CALCUL DU NOMBRE D'ARETES DE LA LIGNE POUR DECLARER LES TABLEAUX
         DTAIL0 = 1D111
         NBARLI = 1
         HHH    = H
         XYZ(1) = XYZ1(1)
         XYZ(2) = XYZ1(2)
         XYZ(3) = XYZ1(3)
C        LES 3 PARAMETRES D'APPEL DE LA FONCTION 'TAILLE_IDEALE'
C        AU SOMMET XYZ DU SEGMENT SOMMET1 SOMMET2
 5       CALL FONVAL( NOFOTI, 3, XYZ,  NCODEV, DTAILL )
         IF( NCODEV .EQ. 0 ) THEN
C
C           NCODEV  : 0 DTAILL N'EST PAS INITIALISEE EN SORTIE
C                     1 DTAILL   EST     INITIALISEE EN SORTIE
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='FONCTION TAILLE_IDEALE INCALCULABLE AU POINT'
            ELSE
             KERR(1)='EDGE_LENGTH(X,Y,Z) FUNCTION UNCOMPUTABLE at POINT'
            ENDIF
            KERR(2) = ' '
            WRITE(KERR(2)(1:15), '(E15.7)') XYZ(1)
            WRITE(KERR(2)(18:32),'(E15.7)') XYZ(2)
            WRITE(KERR(2)(35:49),'(E15.7)') XYZ(3)
            CALL LEREUR
            IERR = 1
            GOTO 9999
C
         ELSE
C
C           TAILLE AUTOUR DU SOMMET EST INITIALISEE
C           ESSAI DE CREER UNE ARETE S1P DE CETTE TAILLE
            DTAILL = ABS( DTAILL )
C           MIN ajoute pour corriger un arret trop rapide du calcul
C           des sommets   8/2/2008
C           DTAIL0 est la taille de la derniere arete calculee
            IF( MIN(DTAIL0,DTAILL) .LT. HHH*0.65 ) THEN
C              CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
               IF( NBARLI .GE. MXTAIL ) THEN
C                 SATURATION DU TABLEAU TAIL => IL EST AUGMENTE
                  CALL TNMCAU( 'REEL', MXTAIL, MXTAIL+512,
     %                         NBARLI-1, MNTAIL )
                  MXTAIL = MXTAIL + 512
               ENDIF
               NBARLI = NBARLI + 1
               IPAS   = 0
C
C              LE POINT A LA DISTANCE DTAILL DU PRECEDENT
 7             XYZ0(1) = XYZ(1) + DTAILL * (XYZ2(1)-XYZ1(1)) / H
               XYZ0(2) = XYZ(2) + DTAILL * (XYZ2(2)-XYZ1(2)) / H
               XYZ0(3) = XYZ(3) + DTAILL * (XYZ2(3)-XYZ1(3)) / H
               IF( IPAS .EQ. 0 ) THEN
C                 LA TAILLE IDEALE EN LA SECONDE EXTREMITE
                  CALL FONVAL( NOFOTI, 3, XYZ0,  NCODEV, DTAIL2 )
C                 LE MINIMUM DES 2 TAILLES IDEALES
                  DTAIL2 = ABS( DTAIL2 )
C                 LA TAILLE IDEALE AU MILIEU DE L'ARC
                  XYZM(1) = ( XYZ(1) + XYZ0(1) ) * 0.5
                  XYZM(2) = ( XYZ(2) + XYZ0(2) ) * 0.5
                  XYZM(3) = ( XYZ(3) + XYZ0(3) ) * 0.5
                  CALL FONVAL( NOFOTI, 3, XYZM,  NCODEV, DTAIL1 )
C                 LE MINIMUM DES 2 TAILLES IDEALES
                  DTAIL1 = ABS( DTAIL1 )
C
C                 STRATEGIE POUR LE CALCUL DE LA TAILLE FINALE
                  IF( DTAIL1 .LT. DTAILL .AND. DTAIL1 .LT. DTAIL2 ) THEN
C                    LA TAILLE AU MILIEU EST INFERIEURE A CELLES DES EXTREMITES
C                    CALCUL A NOUVEAU AVEC CETTE NOUVELLE TAILLE CAR MINIMUM LOC
                     DTAILL = DTAIL1
                     GOTO 7
                  ENDIF
C                 LA TAILLE MINIMALE EST CELLE D'UNE DES 2 EXTREMITES
C                 LA PLUS FAIBLE EST GARDEE
                  IF( DTAIL2 .LT. DTAILL ) THEN
C                    PAS DE NOUVEAU PASSAGE
                     IPAS   = 1
                     DTAILL = DTAIL2
                     GOTO 7
                  ENDIF
               ENDIF
C
C              LE POINT DEFINITIF
               XYZ(1) = XYZ0(1)
               XYZ(2) = XYZ0(2)
               XYZ(3) = XYZ0(3)
C              LA TAILLE DE L'ARETE NBARLI-1
               RMCN(MNTAIL-2+NBARLI) = REAL( DTAILL )
               DTAIL0 = DTAILL
C              LA LONGUEUR RESTANTE A MAILLER
               HHH = HHH - DTAILL
               GOTO 5
            ENDIF
C
         ENDIF
C
      ELSE IF( DARETE .GT. 0D0 ) THEN
C
C        LA LONGUEUR PAR DEFAUT DES ARETES EST PRISE EN COMPTE
C        -----------------------------------------------------
         D = H / DARETE
         NBARLI = NINT( D )
         IF( NBARLI .LE. 0 ) NBARLI = 1
C
      ENDIF
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE LA LIGNE
C     -----------------------------------------------
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBSOLI = NBARLI + 1
      CALL LXTNDC( NTLXLI, 'XYZSOMMET', 'ENTIER', WYZSOM + 3*NBSOLI )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET',  NTXYZS , MNXYZS )
C
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      MCN( MNXYZS + WNBSOM ) = NBSOLI
C
      IF( NOFOTI .LE. 0 ) THEN
C
C        PAS DE FONCTION TAILLE_IDEALE
C        *****************************
C        ADRESSE-1 DU DEBUT DES COORDONNEES DU 2-EME SOMMET DE LA LIGNE
         MNS = MNXYZS + WYZSOM + 2
         IF( ABS(RAIGEO-1) .LT. 1.E-3 ) THEN
C
C           SOMMETS EQUIDISTANTS
C           ====================
C           N NOMBRE DE POINTS A CALCULER
            N = NBSOLI - 2
C           LE PAS SELON LA COORDONNEE
            DO 8 J=1,3
               PAS(J) = ( XYZ2(J) - XYZ1(J) ) / ( N + 1 )
 8          CONTINUE
C
C           GENERATION DES NBSOLI - 2 SOMMETS INTERMEDIAIRES
            DO 20 I=1,N
               DO 10 J=1,3
                  RMCN( MNS + J ) = REAL( XYZ1(J) + I * PAS(J) )
 10            CONTINUE
               MNS = MNS + 3
 20         CONTINUE
C
         ELSE
C
C           SOMMETS EN PROGRESSION GEOMETRIQUE
C           ==================================
            H  = ( 1D0 - RAIGEO ) / ( 1D0 - RAIGEO ** (NBSOLI-1) )
            DO 25 J=1,3
               PAS(J) = ( XYZ2(J) - XYZ1(J) ) * H
               XYZ(J) = XYZ1(J)
 25         CONTINUE
C           XYZ EST LE POINT DE DEPART DU SEGMENT
            DO 40 I=2,NBSOLI-1
               DO 30 J=1,3
                  XYZ(J) = XYZ(J) + PAS(J) * ( RAIGEO ** (I-2) )
                  RMCN( MNS + J ) = REAL( XYZ( J ) )
 30            CONTINUE
               MNS = MNS + 3
 40         CONTINUE
C
         ENDIF
C
      ELSE
C
C        LA FONCTION 'TAILLE_IDEALE' EXISTE
C        **********************************
C        ADRESSE DE LA PREMIERE COORDONNEE DU PREMIER SOMMET DE LA LIGNE
         MNS = MNXYZS + WYZSOM
         RMCN(MNS  ) = REAL( XYZ1(1) )
         RMCN(MNS+1) = REAL( XYZ1(2) )
         RMCN(MNS+2) = REAL( XYZ1(3) )
C
C        LES AUTRES SOMMETS INTERMEDIAIRES
         DO 60 I=1,NBARLI-1
C           LA TAILLE DE L'ARETE I
            DTAILL = RMCN(MNTAIL-1+I)
C           LES 3 COORDONNEES DU SOMMET I+1
            MNS = MNS + 3
            RMCN(MNS  ) = REAL(RMCN(MNS-3) +DTAILL *(XYZ2(1)-XYZ1(1))/H)
            RMCN(MNS+1) = REAL(RMCN(MNS-2) +DTAILL *(XYZ2(2)-XYZ1(2))/H)
            RMCN(MNS+2) = REAL(RMCN(MNS-1) +DTAILL *(XYZ2(3)-XYZ1(3))/H)
 60      CONTINUE
C
C        DESTRUCTION DU TABLEAU DES TAILLES
         CALL TNMCDS( 'REEL', MXTAIL, MNTAIL )
      ENDIF
C
C     AJOUT DES COORDONNEES DU POINT FINAL
      MNS = MNXYZS + WYZSOM + 3 * NBSOLI - 3
      RMCN( MNS     ) = REAL( XYZ2( 1 ) )
      RMCN( MNS + 1 ) = REAL( XYZ2( 2 ) )
      RMCN( MNS + 2 ) = REAL( XYZ2( 3 ) )
C
C     AJOUT DES COORDONNEES DU POINT INITIAL
      MNS = MNXYZS + WYZSOM
      RMCN( MNS     ) = REAL( XYZ1( 1 ) )
      RMCN( MNS + 1 ) = REAL( XYZ1( 2 ) )
      RMCN( MNS + 2 ) = REAL( XYZ1( 3 ) )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNXYZS) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + WBCOOR ) = 3
      MCN( MNXYZS + WNBTGS ) = 0
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' DE TYPE LIGNE STRUCTUREE
C     -------------------------------------------------------
      CALL LXTNDC( NTLXLI, 'NSEF', 'ENTIER', 1 + WBARSE )
      CALL LXTSOU( NTLXLI, 'NSEF',  NTNSEF , MNNSEF )
C
C     LE TYPE DE L'OBJET : ICI LIGNE
      MCN( MNNSEF + WUTYOB ) = 2
C     LE TYPE DU MAILLAGE : ICI SEGMENT STRUCTURE
      MCN( MNNSEF + WUTYMA ) = 2
C     LE NOMBRE DE SOMMETS ET TANGENTES PAR ARETE
      MCN( MNNSEF + WBSOEF ) = 2
      MCN( MNNSEF + WBTGEF ) = 0
C     PAS D'EF A TG
      MCN( MNNSEF + WBEFTG ) = 0
      MCN( MNNSEF + WBEFAP ) = 0
C     LE NOMBRE D'ARETES DU SEGMENT STRUCTURE
      MCN( MNNSEF + WBEFOB ) = NBARLI
C     LE NOMBRE D'ARETES DU SEGMENT STRUCTURE
      MCN( MNNSEF + WBARSE ) = NBARLI
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNNSEF + WUTFMA ) = 0
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNNSEF) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C
C     PAS D'ERREUR RENCONTREE
      IERR = 0
C
C     ERREUR
C     ======
 9999 RETURN
      END
