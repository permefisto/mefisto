      SUBROUTINE VECPRES( NDIM,   NUTYEL, NBELEM, NBNOEF, NUNDEL,
     %                    NBNOVI, NDDL,
     %                    NTDL,   NBVECT, VXYZPN,
     %                    NBNOPR, PR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECUPERER DANS UN TMC LES NBVECT VECTEURS PRESSION P1
C ----- EN TOUS LES SOMMETS DU MAILLAGE POUR LE TRIANGLE ou TETRAEDRE
C       DE TAYLOR HOOD ou BREZZI-FORTIN
C
C ENTREES :
C ---------
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C NUTYEL : NO DU TYPE D'EF
C NBELEM : NOMBRE DES EF
C NBNOEF : NOMBRE DE NOEUD d'UN EF
C NUNDEL : NO DES NBNOEF NOEUDS DES NBELEM EF
C NDDL   : NDDL(I)=NO DU DERNIER DL DU NOEUD VITESSE I NDDL(0)=0
C NTDL   : NOMBRE TOTAL DE DL (VITESSES+PRESSIONS)
C NBVECT : NOMBRE TOTAL DE VECTEURS VITESSE+PRESSION
C VXYZPN : TABLEAU VITESSEPRESSION SOLUTION CONNUE DL PAR NOEUDS
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C         (POUR TAYLOR HOOD CE SONT LES SOMMETS ET MILIEUX DES ARETES)
C         (POUR BREZZI FORTIN CE SONT LES SOMMETS ET BARYCENTRES DES EF
C          ET IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBSOM=NBNOPR
C          ET LES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF)
C NBNOPR : NOMBRE DE SOMMETS SUPPORTS DE LA PRESSION
C         (POUR TAYLOR HOOD   CE SONT LES SOMMETS ET MILIEUX DES ARETES)
C         (POUR BREZZI FORTIN CE SONT LES SOMMETS)
C
C SORTIES:
C --------
C PR     : PRESSION EN CHAQUE NOEUD DU MAILLAGE ET NBVECT FOIS
C          ATTENTION:  PR A UNE VALEUR AUX NOEUDS C'EST A DIRE
C          POUR TAYLOR HOOD:   AUX SOMMETS ET MILIEUX DES ARETES DES EF
C                              AUX MILIEUX C'EST LA MOYENNE AUX 2 SOMMETS
C          POUR BREZZI FORTIN: AUX SOMMETS DES EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Janvier 2011
C23456---------------------------------------------------------------012
      INTEGER           NDDL(0:NBNOVI)
      INTEGER           NUNDEL(NBELEM,NBNOEF)
      DOUBLE PRECISION  VXYZPN(NTDL,NBVECT)
      DOUBLE PRECISION  PR(NBNOPR,NBVECT)
C
C     PARTAGE DES VITESSES EN VX, VY, VZ et de la PRESSION AUX SOMMETS
      IF( NUTYEL .EQ. 15 .OR. NUTYEL .EQ. 20 ) THEN
C
C        TAYLOR-HOOD
C        ===========
         DO I=1,NBNOVI
C
C           LE NUMERO DU DL AVANT LE NOEUD I
            NODGLI = NDDL(I-1)
C
C           LE NOMBRE DE DL DU NOEUD I
            NDL = NDDL(I) - NODGLI
C
            IF( NDL .EQ. NDIM+1 ) THEN
               DO K=1,NBVECT
C                 LA PRESSION EN CE NOEUD=SOMMET
                  PR( I, K ) = VXYZPN( NODGLI+NDIM+1, K )
               ENDDO
            ENDIF
C
         ENDDO
C
C        PARCOURS DES EF POUR METTRE LA PRESSION AU MILIEU DES
C        ARETES COMME MOYENNE DES PRESSIONS DE SES 2 SOMMETS
         IF( NDIM .EQ. 2 ) THEN
C
C           TRIANGLE P2 P1 de TAYLOR-HOOD
C           -----------------------------
            DO NEF=1,NBELEM
C
C              BOUCLE SUR LES 3 ARETES DU TRIANGLE NEF
               DO I=1,3
C
C                 NO DU NOEUD MILIEU DE L'ARETE I
                  NM = NUNDEL(NEF,I+3)
C
C                 NO DES 2 SOMMETS DE L'ARETE I
                  N1 = NUNDEL(NEF,I)
                  IF( I .LT. 3 ) THEN
                     J = I+1
                  ELSE
                     J = 1
                  ENDIF
                  N2 = NUNDEL(NEF,J)
C
C                 LA PRESSION EN CE NOEUD=MILIEU D'UNE ARETE
                  DO K=1,NBVECT
                     PR( NM, K ) = ( PR( N1, K ) + PR( N2, K ) ) * 0.5D0
                  ENDDO
               ENDDO
            ENDDO
C
         ELSE
C
C           TETRAEDRE P2 P1 de TAYLOR-HOOD
C           ------------------------------
            DO NEF=1,NBELEM
C
C              NO DU SOMMET 4 DU TETRAEDRE
               N4 = NUNDEL(NEF,4)
C
C              BOUCLE SUR LES 6 ARETES DU TETRAEDRE NEF
               DO I=1,3
C
C                 NO DU NOEUD MILIEU DE L'ARETE I DU BAS
                  NMB = NUNDEL(NEF,I+4)
C
C                 NO DU NOEUD MILIEU DE L'ARETE I+3 DU HAUT
                  NMH = NUNDEL(NEF,I+7)
C
C                 NO DES 2 SOMMETS DE L'ARETE I DU BAS
                  N1 = NUNDEL(NEF,I)
                  IF( I .LT. 3 ) THEN
                     J = I+1
                  ELSE
                     J = 1
                  ENDIF
                  N2 = NUNDEL(NEF,J)
C
C                 LA PRESSION EN CE NOEUD=MILIEU D'UNE ARETE
                  DO K=1,NBVECT
                     PR( NMB, K ) = ( PR( N1, K ) + PR( N2, K ) ) *0.5D0
                     PR( NMH, K ) = ( PR( N1, K ) + PR( N4, K ) ) *0.5D0
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
C
      ELSE IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN
C
C        EF TRIANGLE ou TETRAEDRE BREZZI-FORTIN
C        --------------------------------------
C        NBNOVI = SOMMETS + NOEUDS BARYCENTRE DE CHAQUE EF EN VITESSE
C        ATTENTION: IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBNOPR
C                   LES AUTRES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF
C        PAR SUITE, DANS CE QUI SUIT, SEULS LES SOMMETS DONNENT UNE PRESSION
         DO I=1,NBNOPR
C
C           LE NUMERO DU DL AVANT LE NOEUD I
            NODGLI = NDDL(I-1)
C
C           LE NOMBRE DE DL PRESSION AU NOEUD I
            NDLP = NDDL(I) - NODGLI
C
            IF( NDLP .GT. NDIM ) THEN
               DO K=1,NBVECT
C                 LA PRESSION EN CE NOEUD=SOMMET DE NO AVANT LES BARYCENTRES
                  PR( I, K ) = VXYZPN( NODGLI+NDIM+1, K )
               ENDDO
            ENDIF
         ENDDO
C
      ENDIF
C
      RETURN
      END
