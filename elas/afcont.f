      SUBROUTINE AFCONT( NOFORE, NOAXIS, NOMELE,
     %                   NBELFI, NBPIEF, NDIMES ,NBCAS,
     %                   COPIEF, CONPRI, DIRPRI,
     %                   NCAS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES CONTRAINTES PRINCIPALES D'UN ENSEMBLE D'EF
C -----
C
C ENTREES :
C ---------
C NOFORE : 0 SI PAS DE FONCTION UTILISATEUR 'REGION', >0 SINON
C NOAXIS : 1 SI PB 2D AXISYMETRIQUE, 0 SINON
C NOMELE : NOM DU TYPE DES EF TRAITES ICI
C NBELFI : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NBPIEF : NOMBRE DE POINTS DE CALCUL DES CONTRAINTES PAR ELEMENT
C NDIMES : ESPACE DE TRAVAIL 2 OU 3
C NBCAS  : NOMBRE DE CAS TRAITES
C COPIEF : LES NDIMES COORDONNEES DES POINTS DE CALCUL DES CONTRAINTES
C CONPRI : LES NDIMES CONTRAINTES PRINCIPALES
C DIRPRI : LES NDIMES DIRECTIONS PRINCIPALES
C NCAS   : LE NUMERO DU CAS A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JUILLET 2002
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ctemps.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL              COPIEF(1:NBELFI,1:NBPIEF,1:NDIMES)
      DOUBLE PRECISION  CONPRI(1:NBELFI,1:NBPIEF,1:NDIMES,1:NBCAS)
      REAL     DIRPRI(1:NBELFI,1:NBPIEF,1:NDIMES,1:NDIMES,1:NBCAS)
      DOUBLE PRECISION  VMAX,A(4),CONMIN,CONMAX,DINFO
      CHARACTER*4       NOMELE(2)
      INTRINSIC         REAL
C
      WRITE(IMPRIM,*) 'CAS',NCAS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'CONTRAINTES des EF de TYPE ',NOMELE(1),
     %                   ' ',NOMELE(2)
         WRITE(IMPRIM,*) 'CP pour CONTRAINTE PRINCIPALE'
         WRITE(IMPRIM,*) 'DP pour DIRECTION PRINCIPALE'
      ELSE
         WRITE(IMPRIM,*) 'MAIN STRESSES and DIRECTIONS of FE of TYPE ',
     %                    NOMELE(1),' ',NOMELE(2)
         WRITE(IMPRIM,*) 'MS for MAIN STRESS'
         WRITE(IMPRIM,*) 'MD for MAIN DIRECTION'
      ENDIF
C
      CONMIN =  DINFO( 'GRAND' )
      CONMAX = -CONMIN
      ZMIN = 0
      ZMAX = 0
      DO 100 K=1,NBELFI
         DO 50 L=1,NBPIEF
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
            IF( VP1 .GT. CONMAX ) THEN
               CONMAX = VP1
               XMAX   = XF
               YMAX   = YF
            ENDIF
            IF( VP1 .LT. CONMIN ) THEN
               CONMIN = VP1
               XMIN   = XF
               YMIN   = YF
            ENDIF
C
C           LA SECONDE CONTRAINTE PRINCIPALE
            VP2 = REAL( CONPRI(K,L,2,NCAS) )
C           SA DIRECTION PRINCIPALE
            XV2 = DIRPRI(K,L,1,2,NCAS)
            YV2 = DIRPRI(K,L,2,2,NCAS)
C
            IF( VP2 .GT. CONMAX ) THEN
               CONMAX = VP2
               XMAX   = XF
               YMAX   = YF
            ENDIF
            IF( VP2 .LT. CONMIN ) THEN
               CONMIN = VP2
               XMIN   = XF
               YMIN   = YF
            ENDIF
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
            IF( J .NE. 0 ) THEN
               IF( NDIMES .EQ. 2 ) THEN
                  IF( NOAXIS .NE. 0 ) THEN
C
C                    PB AXISYMETRIQUE
                     IF( LANGAG .EQ. 0 ) THEN
                        WRITE(IMPRIM,10010) K,L,XF,YF,NCAS,
     %                                      VP1,VP2,XV1,YV1,XV2,YV2
                     ELSE
                        WRITE(IMPRIM,20010) K,L,XF,YF,NCAS,
     %                                      VP1,VP2,XV1,YV1,XV2,YV2
                     ENDIF
10010                FORMAT(' EF',I5,' : PT',I3,' R=',G12.4,' Z=',G12.4,
     %                      ' CAS',I2,': CP1',G12.4,' CP2',G12.4,
     %                      '  DP1',2F7.3,' DP2',2F7.3)
20010                FORMAT(' FE',I5,' : PT',I3,' R=',G12.4,' Z=',G12.4,
     %                      ' CASE',I2,': MS1',G12.4,' MS2',G12.4,
     %                      '  MD1',2F7.3,' MD2',2F7.3)
C
                  ELSE
C
C                    PB 2D
                     IF( LANGAG .EQ. 0 ) THEN
                        WRITE(IMPRIM,10020) K,L,XF,YF,NCAS,
     %                                      VP1,VP2,XV1,YV1,XV2,YV2
                     ELSE
                        WRITE(IMPRIM,20020) K,L,XF,YF,NCAS,
     %                                      VP1,VP2,XV1,YV1,XV2,YV2
                     ENDIF
10020                FORMAT(' EF',I5,' : PT',I3,' X=',G12.4,' Y=',G12.4,
     %                      ' CAS',I2,': CP1',G12.4,' CP2',G12.4,
     %                      '  DP1',2F7.3,' DP2',2F7.3)
20020                FORMAT(' FE',I5,' : PT',I3,' X=',G12.4,' Y=',G12.4,
     %                      ' CASE',I2,': MS1',G12.4,' MS2',G12.4,
     %                      '  MD1',2F7.3,' MD2',2F7.3)
C
                  ENDIF
               ELSE
C
C                 PB 3D
                  ZF  = COPIEF(K,L,3)
                  VP3 = REAL( CONPRI(K,L,3,NCAS) )
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10030) K,L,XF,YF,ZF,
     %                                   NCAS,VP1,VP2,VP3,
     %                 ((DIRPRI(K,L,M,N,NCAS),M=1,3),N=1,3)
                  ELSE
                     WRITE(IMPRIM,20030) K,L,XF,YF,ZF,
     %                                   NCAS,VP1,VP2,VP3,
     %                 ((DIRPRI(K,L,M,N,NCAS),M=1,3),N=1,3)
                  ENDIF
10030             FORMAT(' EF',I5,': PT',I3,' X=',G12.4,' Y=',G12.4,
     %                   ' Z=',G12.4/' CAS',I2,
     %                   ': CP1',G12.4,' CP2',G12.4,' CP3',G12.4,
     %                   ' DP1',3F7.3,' DP2',3F7.3,' DP3',3F7.3)
20030             FORMAT(' FE',I6,': PT',I3,' X=',G12.4,' Y=',G12.4,
     %                   ' Z=',G12.4/' CASE',I2,
     %                   ': MS1',G12.4,' MS2',G12.4,' MS3',G12.4,
     %                   ' MD1',3F7.3,' MD2',3F7.3,' MD3',3F7.3)
C
                  IF( VP3 .GT. CONMAX ) THEN
                     CONMAX = VP3
                     XMAX   = XF
                     YMAX   = YF
                     ZMAX   = ZF
                  ENDIF
                  IF( VP3 .LT. CONMIN ) THEN
                     CONMIN = VP3
                     XMIN   = XF
                     YMIN   = YF
                     ZMIN   = ZF
                  ENDIF
               ENDIF
            ENDIF
 50      CONTINUE
 100  CONTINUE
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10110) NOMELE,CONMIN,XMIN,YMIN,ZMIN
         WRITE(IMPRIM,10120) NOMELE,CONMAX,XMAX,YMAX,ZMAX
      ELSE
         WRITE(IMPRIM,20110) NOMELE,CONMIN,XMIN,YMIN,ZMIN
         WRITE(IMPRIM,20120) NOMELE,CONMAX,XMAX,YMAX,ZMAX
      ENDIF
10110 FORMAT(/' La CONTRAINTE PRINCIPALE minimale DES EF ',A4,' ',A4,
     %' VAUT',G15.7,' AU POINT X=',G13.5, '  Y=',G13.5,'  Z=',G13.5)
20110 FORMAT(/' The minimum of MAIN STRESSES of FE ',A4,' ',A4,
     %' =',G15.7,' at POINT X=',G13.5, '  Y=',G13.5,'  Z=',G13.5)
C
10120 FORMAT(' La CONTRAINTE PRINCIPALE MAXIMALE DES EF ',A4,' ',A4,
     %' VAUT',G15.7,' AU POINT X=',G13.5, '  Y=',G13.5,'  Z=',G13.5/)
20120 FORMAT(' The MAXIMUM of MAIN STRESSES of FE ',A4,' ',A4,
     %' =',G15.7,' at POINT X=',G13.5, '  Y=',G13.5,'  Z=',G13.5/)
      RETURN
      END
