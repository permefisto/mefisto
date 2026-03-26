      SUBROUTINE VITROT( NBNOVI, XYZNOE, NDDLNO, NTDLVP, NBVECT, TEMPSC,
     %                   NOVOLU, NUMIVO, NUMAVO, LTDEVO, WITPRE, VXYZPN)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LA ROTATION A LA VITESSE W POUR OBTENIR LA VITESSE
C -----    CARTESIENNE V(t,X,Y,Z) POUR TOUS LES TEMPS ET TOUS LES VECTEURS
C          Vx = Wx + omega2 * z - omega3 * y
C          Vy = Wy + omega3 * x - omega1 * z
C          Vz = Wz + omega1 * y - omega2 * x
C ATTENTION: LA VITESSE ANGULAIRE OMEGA(3) EST SUPPOSEE CONSTANTE
C
C ENTREES:
C --------
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE ET PRESSION
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C XYZNOE : XYZ DES NOEUDS
C NDDLNO : TABLEAU DES POINTEURS SUR DERNIER DL POUR CHAQUE NOEUD
C          NDDLNO(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES+PRESSIONS
C NBVECT : NOMBRE DE VECTEURS VITESSES+PRESSIONS AUX NOEUDS
C TEMPSC : TEMPS DU CALCUL POUR LES NBVECT VECTEURS
C
C NOVOLU : NUMERO DU VOLUME
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C WITPRE : WITESSE PRESSION EN TOUS TEMPS DE CALCUL ET TOUS LES NOEUDS
C
C SORTIE :
C --------
C VXYZPN : VITESSE PRESSION EN TOUS TEMPS DE CALCUL ET TOUS LES NOEUDS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Aout 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION  WITPRE( NTDLVP, NBVECT ),
     %                  VXYZPN( NTDLVP, NBVECT )
      REAL              XYZNOE( 3, NBNOVI ), TEMPSC( NBVECT )
      INTEGER           NDDLNO( 0:NBNOVI )
      INTEGER           LTDEVO(1:MXDOFL, NUMIVO:NUMAVO )
      DOUBLE PRECISION  XD, YD, ZD, OMEGA(3)
C
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) '==========================================='
      WRITE(IMPRIM,*) '     V(t,xyz) = W(t,xyz) + omega  x  xyz'
      WRITE(IMPRIM,*) '==========================================='
C
C     RECUPERATION VITESSE ANGULAIRE SUPPOSEE CONSTANTE EN XYZ
C     MAIS PAS EN TEMPS
      XD = XYZNOE( 1, 1 )
      YD = XYZNOE( 2, 1 )
      ZD = XYZNOE( 3, 1 )
C
      DO K = 1, NBVECT
C
C        LE TEMPS DE CALCUL DU VECTEUR K
         TEMPS = TEMPSC( K )
C
C        LA VITESSE ANGULAIRE Omega A CE TEMPS
         CALL REVIAN( 4, NOVOLU, XD, YD, ZD, LTDEVO(LPVIAN,NOVOLU),
     %                OMEGA )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Au temps',TEMPS,' Omega=',OMEGA
         ELSE
            WRITE(IMPRIM,*) 'At time',TEMPS,' Omega=',OMEGA
         ENDIF
C
         DO I = 1, NBNOVI
C
C           LE NO DU DL AVANT LES DL VITESSE DU NOEUD I
C           LA PRESSION EST ENSUITE SI I EST UN SOMMET
            NDL = NDDLNO(I-1)
C
C           XYZ DU NOEUD I
            XD = XYZNOE( 1, I )
            YD = XYZNOE( 2, I )
            ZD = XYZNOE( 3, I )
C
C           Vx = Wx + omega2 * z - omega3 * y
            VXYZPN( NDL+1, K ) = WITPRE( NDL+1, K )
     %                         + OMEGA(2) * ZD - OMEGA(3) * YD
C
C           Vy = Wy + omega3 * x - omega1 * z
            VXYZPN( NDL+2, K ) = WITPRE( NDL+2, K )
     %                         + OMEGA(3) * XD - OMEGA(1) * ZD
C
C           Vz = Wz + omega1 * y - omega2 * x
            VXYZPN( NDL+3, K ) = WITPRE( NDL+3, K )
     %                         + OMEGA(1) * YD - OMEGA(2) * XD
C
            IF( NDDLNO(I)-NDL .EQ. 4 ) THEN
C              PRESSION IDENTIQUE
               VXYZPN( NDL+4, K ) = WITPRE( NDL+4, K )
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
