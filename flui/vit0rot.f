      SUBROUTINE VIT0ROT( OMEGA, NBNOVI, XYZNOE, NDDLNO, NTDLVP, VXYZPN,
     %                    WXYZPN)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETRANCHER LA ROTATION A LA VITESSE CARTESIENNE V(t0,X,Y,Z)
C -----    POUR OBTENIR LA VITESSE  W(t0,X,Y,Z) = V(t0,X,Y,Z) - Omega x XYZ
C          Wx = Vx - ( omega2 * z - omega3 * y )
C          Wy = Vy - ( omega3 * x - omega1 * z )
C          Wz = Vz - ( omega1 * y - omega2 * x )

C ATTENTION: LA VITESSE ANGULAIRE OMEGA(1:3) EST SUPPOSEE CONSTANTE
C
C ENTREES:
C --------
C OMEGA  : VITESSE ANGULAIRE DANS LES 3 DIRECTIONS X Y Z
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE ET PRESSION
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C XYZNOE : XYZ DES NOEUDS
C NDDLNO : TABLEAU DES POINTEURS SUR DERNIER DL POUR CHAQUE NOEUD
C          NDDLNO(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES+PRESSIONS
C VXYZPN : VITESSE PRESSION EN TOUS LES NOEUDS AU TEMPS t0
C
C SORTIE :
C --------
C WXYZPN : VITESSE PRESSION EN TOUS LES NOEUDS - Omega x XYZ AU TEMPS t0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Juillet 2013
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION  OMEGA(3),
     %                  VXYZPN( NTDLVP ),
     %                  WXYZPN( NTDLVP ),
     %                  XD, YD, ZD
      REAL              XYZNOE( 3, NBNOVI )
      INTEGER           NDDLNO( 0:NBNOVI )
C
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) '==========================================='
      WRITE(IMPRIM,*) '  W(t0,xyz) = V(t0,xyz) - Omega  x  xyz'
      WRITE(IMPRIM,*) '==========================================='
C
      DO I = 1, NBNOVI
C
C        LE NO DU DL AVANT LES DL VITESSE DU NOEUD I
C        LA PRESSION EST ENSUITE SI I EST UN SOMMET
         NDL = NDDLNO(I-1)
C
C        XYZ DU NOEUD I
         XD = XYZNOE( 1, I )
         YD = XYZNOE( 2, I )
         ZD = XYZNOE( 3, I )
C
C        Vx = Wx - ( omega2 * z - omega3 * y )
         WXYZPN( NDL+1 ) = VXYZPN( NDL+1 )
     %                   - OMEGA(2) * ZD + OMEGA(3) * YD
C
C        Vy = Wy - ( omega3 * x - omega1 * z )
         WXYZPN( NDL+2 ) = VXYZPN( NDL+2 )
     %                   - OMEGA(3) * XD + OMEGA(1) * ZD
C
C        Vz = Wz - ( omega1 * y - omega2 * x )
         WXYZPN( NDL+3 ) = VXYZPN( NDL+3 )
     %                   - OMEGA(1) * YD + OMEGA(2) * XD
C
         IF( NDDLNO(I)-NDL .EQ. 4 ) THEN
C           PRESSION IDENTIQUE
            WXYZPN( NDL+4 ) = VXYZPN( NDL+4 )
         ENDIF

      ENDDO
C
      RETURN
      END
