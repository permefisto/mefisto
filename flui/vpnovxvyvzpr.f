      SUBROUTINE VPNOVXVYVZPR( NDIM,   NUTYEL, NBELEM, NBNOEL, NONOEF,
     %                         NDDLNO, NTDL,   VXYZPN,
     %                         NBNOVI, VX, VY, VZ,     NBNOPR, PR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     DECOUPER LE VECTEUR VITESSE+PRESSION EN 4 TABLEAUX
C -----     INDEPENDANTS VITESSE en X, puis en Y, puis en Z puis en PRESSION
C           EN TOUS LES NOEUDS DU MAILLAGE POUR LE TRIANGLE ou TETRAEDRE
C           DE TAYLOR HOOD ou BREZZI-FORTIN

C ATTENTION: EF BREZZI-FORTIN PRESSION AUX SOMMETS DES EF
C            EF TAYLOR-HOOD   PRESSION AUX SOMMETS et MILIEUX des ARETES DES EF
C                             PRESSION(MILIEU)=(PRESSION(ST1)+PRESSION(ST2))/2
C                             ST1 et ST2 LES 2 SOMMETS DE L'ARETE

C REMARQUE: SP IDENTIQUE a vxvyvzpr.f MAIS avec NBVECT=1 => MOINS D'OPERATIONS

C ENTREES :
C ---------
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C NUTYEL : NO DU TYPE D'EF
C NBELEM : NOMBRE DES EF
C NBNOEL : NOMBRE DE NOEUD d'UN EF
C NONOEF : NO DES NBNOEL NOEUDS DES NBELEM EF
C NDDLNO   : NDDLNO(I)=NO DU DERNIER DL DU NOEUD VITESSE I NDDLNO(0)=0
C NTDL   : NOMBRE TOTAL DE DL (VITESSES+PRESSIONS)
C VXYZPN : TABLEAU VITESSEPRESSION SOLUTION CONNUE DL PAR NOEUDS
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C         (POUR TAYLOR HOOD CE SONT LES SOMMETS ET MILIEUX DES ARETES)
C         (POUR BREZZI FORTIN CE SONT LES SOMMETS ET BARYCENTRES DES EF
C          ET IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBSOM=NBNOPR
C          LES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF)
C NBNOPR : NOMBRE DE NOEUDS SUPPORT DE LA PRESSION
C         (POUR BREZZI-FORTIN CE SONT LES SOMMETS
C          POUR TAYLOR-HOOD   CE SONT LES SOMMETS+MILIEUX CAR
C          LA VALEUR AU MILIEU DES ARETES EST LA DEMI SOMME DES
C          VALEURS AUX SOMMETS POUR HOMOGENEITE AVEC LES VITESSES)

C SORTIES:
C --------
C VX     : LA VITESSE EN X EN CHAQUE NOEUD VITESSE
C VY     : LA VITESSE EN Y EN CHAQUE NOEUD VITESSE
C VZ     : LA VITESSE EN Z EN CHAQUE NOEUD VITESSE
C PR     : LA PRESSION     EN CHAQUE NOEUD DU MAILLAGE
C          ATTENTION:  PR A UNE VALEUR
C          POUR TAYLOR HOOD   AUX SOMMETS ET MILIEUX DES ARETES DES EF
C          POUR BREZZI FORTIN AUX SOMMETS DES EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Laboratoire J-L.LIONS UPMC Paris   Juin 2007
C MODIFS : ALAIN PERRONNET  Saint Pierre du Perray          Fevrier 2021
C23456---------------------------------------------------------------012
      INTEGER           NDDLNO(0:NBNOVI)
      INTEGER           NONOEF(NBELEM,NBNOEL)
      DOUBLE PRECISION  VXYZPN(NTDL)
      DOUBLE PRECISION  VX(NBNOVI),
     %                  VY(NBNOVI),
     %                  VZ(NBNOVI)
      DOUBLE PRECISION  PR(NBNOPR)

C     PARTAGE DES VITESSES EN VX, VY, VZ et de la PRESSION AUX SOMMETS
      NDIM1 = NDIM + 1

      IF( NUTYEL .EQ. 15 .OR. NUTYEL .EQ. 20 ) THEN

C        TAYLOR-HOOD
C        ===========
         DO I=1,NBNOVI

C           LE NUMERO DU DL AVANT LE NOEUD I
            NODGLI = NDDLNO(I-1)

C           LE NOMBRE DE DL DU NOEUD I
            NDL = NDDLNO(I) - NODGLI

C              COMPOSANTE 1 DE LA VITESSE
               VX( I ) = VXYZPN( NODGLI+1 )

C              COMPOSANTE 2 DE LA VITESSE
               VY( I ) = VXYZPN( NODGLI+2 )

               IF( NDIM .EQ. 3 ) THEN
C                 COMPOSANTE 3 DE LA VITESSE
                  VZ( I ) = VXYZPN( NODGLI+3 )
               ENDIF

               IF( NDL .EQ. NDIM1 ) THEN
C                 LA PRESSION EN CE NOEUD=SOMMET
                  PR( I ) = VXYZPN( NODGLI+NDIM1 )
               ENDIF

         ENDDO

C        PARCOURS DES EF POUR METTRE LA PRESSION AU MILIEU DES
C        ARETES COMME MOYENNE DES PRESSIONS DE SES 2 SOMMETS
         IF( NDIM .EQ. 2 ) THEN

C           TRIANGLE P2 P1 de TAYLOR-HOOD
C           -----------------------------
            DO NEF=1,NBELEM

C              BOUCLE SUR LES 3 ARETES DU TRIANGLE NEF
               DO I=1,3

C                 NO DU NOEUD MILIEU DE L'ARETE I
                  NM = NONOEF(NEF,I+3)

C                 NO DES 2 SOMMETS DE L'ARETE I
                  N1 = NONOEF(NEF,I)
                  IF( I .LT. 3 ) THEN
                     J = I+1
                  ELSE
                     J = 1
                  ENDIF
                  N2 = NONOEF(NEF,J)
                  PR( NM ) = ( PR( N1 ) + PR( N2 ) ) * 0.5D0

               ENDDO
            ENDDO

         ELSE

C           TETRAEDRE P2 P1 de TAYLOR-HOOD
C           ------------------------------
            DO NEF=1,NBELEM

C              NO DU SOMMET 4 DU TETRAEDRE
               N4 = NONOEF(NEF,4)

C              BOUCLE SUR LES 6 ARETES DU TETRAEDRE NEF
               DO I=1,3

C                 NO DU NOEUD MILIEU DE L'ARETE I DU BAS
                  NMB = NONOEF(NEF,I+4)

C                 NO DU NOEUD MILIEU DE L'ARETE I+3 DU HAUT
                  NMH = NONOEF(NEF,I+7)

C                 NO DES 2 SOMMETS DE L'ARETE I DU BAS
                  N1 = NONOEF(NEF,I)
                  IF( I .LT. 3 ) THEN
                     J = I+1
                  ELSE
                     J = 1
                  ENDIF
                  N2 = NONOEF(NEF,J)
                  PR( NMB ) = ( PR( N1 ) + PR( N2 ) ) *0.5D0
                  PR( NMH ) = ( PR( N1 ) + PR( N4 ) ) *0.5D0

               ENDDO
            ENDDO

         ENDIF

      ELSE IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN

C        EF TRIANGLE ou TETRAEDRE BREZZI-FORTIN
C        --------------------------------------
C        NBNOVI = SOMMETS + NOEUDS BARYCENTRE DE CHAQUE EF EN VITESSE
C        ATTENTION: IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBNOPR
C                   LES AUTRES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF
         DO I=1,NBNOVI

C           LE NUMERO DU DL AVANT LE NOEUD I
            NODGLI = NDDLNO(I-1)

C           LE NOMBRE DE DL PRESSION AU NOEUD I
            NDLP = NDDLNO(I) - NODGLI - NDIM

C           COMPOSANTE 1 DE LA VITESSE
            VX( I ) = VXYZPN( NODGLI+1 )

C           COMPOSANTE 2 DE LA VITESSE
            VY( I ) = VXYZPN( NODGLI+2 )

            IF( NDIM .EQ. 3 ) THEN
C              COMPOSANTE 3 DE LA VITESSE
               VZ( I ) = VXYZPN( NODGLI+3 )
            ENDIF

            IF( NDLP .GT. 0 ) THEN
C              LA PRESSION EN CE NOEUD=SOMMET DE NO AVANT LES BARYCENTRES
               PR( I ) = VXYZPN( NODGLI+NDIM1 )
            ENDIF

         ENDDO

      ENDIF

      RETURN
      END
