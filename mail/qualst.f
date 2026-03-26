      SUBROUTINE QUALST( NS,     XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                   VOLUMT, QUALIT,
     %                   VOLUNS, QUALNS, NTQMIN, NBTENS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU VOLUME ET DE LA QUALITE D'UN SOMMET D'UNE TETRAEDRISATION
C -----    EGALE AU MINIMUM DES QUALITES DES TETRAEDRES DE CE SOMMET
C
C ENTREES:
C --------
C NS     : NUMERO DU SOMMET DE QUALITE DES TETRAEDRES A CALCULER
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C MODIFIES :
C ----------
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VOLUNS : VOLUME DES TETRAEDRES DE SOMMET NS
C QUALNS : QUALITE DU SOMMET NS DE LA TETRAEDRISATION
C NTQMIN : NUMERO NSTETR DU TETRAEDRE DE SOMMET NS ET DE QUALITE MINIMALE
C NBTENS : NOMBRE DE TETRAEDRES DE SOMMET NS
C          0 SI NS NON SOMMET D'UN TETRAEDRE CHAINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1993
C....................................................................012
      INTEGER           NSTETR(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*)
      REAL              XYZSOM(3,*),
     %                  VOLUMT(*),
     %                  QUALIT(*)

      REAL              VOLUNS, QUALNS, ARMIN, ARMAX, SURFTR(4)
C
      NBTENS = 0
      VOLUNS = 0.0
      QUALNS = 100.0
C
C     POSITION DANS NOTESO DU 1-ER TETRAEDRE DE SOMMET NS
      NDT = NO1TSO( NS )

C     TANT QU'IL EXISTE UN TETRAEDRE DE SOMMET NS FAIRE
 10   IF( NDT .GT. 0 ) THEN
C
C        LE NUMERO DU TETRAEDRE DANS NSTETR
         NT = NOTESO(1,NDT)

         if( nstetr(1,nt) .le. 0 ) then
            goto 30
         endif

c        verification du chainage des tetraedres du sommet NS
C        sequence a supprimer si plus de detection d'anomalie!...
         do j=1,4
            if( nstetr(j,nt) .eq. ns ) goto 20
         enddo
         print *,'qualst: ANOMALIE ns=',ns,' NON SOMMET de NSTETR(',
     %   nt,')=',(nstetr(kk,nt),kk=1,4)
         NBTENS = 0
         GOTO 9999
C
C        CALCUL DU VOLUME ET QUALITE MINIMALE DES TETRAEDRES DE SOMMET NS
 20      CALL QUATET( XYZSOM(1,NSTETR(1,NT)),
     %                XYZSOM(1,NSTETR(2,NT)),
     %                XYZSOM(1,NSTETR(3,NT)),
     %                XYZSOM(1,NSTETR(4,NT)),
     %   ARMIN, ARMAX, SURFTR, VOLUMT(NT), QUALIT(NT) )
C
         VOLUNS = VOLUNS + ABS( VOLUMT(NT) )
         IF( QUALIT(NT) .LT. QUALNS ) THEN
            QUALNS = QUALIT(NT)
            NTQMIN = NT
         ENDIF
C
C        UN TETRAEDRE DE PLUS CONTIENT LE SOMMET NS
         NBTENS = NBTENS + 1
C
C        LE TETRAEDRE SUIVANT
 30      NDT = NOTESO(2,NDT)
         GOTO 10
      ENDIF
C
 9999 RETURN
      END
