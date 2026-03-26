      SUBROUTINE TETSPH( NCOGEL, NSEFJ, TE1SOB, PLUTET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PERMUTATION DES NUMEROS DE SOMMETS DE L'EF POUR QUE LE PLUS
C -----    PETIT NUMERO DE SOMMET SOIT LE PREMIER SOMMET DE L'EF
C
C ENTREE :
C --------
C NCOGEL : CODE GEOMETRIQUE DE L'EF
C          5:TETRAEDRE, 6:PENTAEDRE, 7:HEXAEDRE, 9:PYRAMIDE A BASE CARREE
C NSEFJ  : SOMMETS D'UN EF
C
C SORTIE :
C --------
C TE1SOB : TETRAEDRISATION D'UN EF
C PLUTET : NOMBRE DE TETRAEDRES DANS L'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: F.CHOUKROUN O.RICOU UPMC ANALYSE NUMERIQUE PARIS JANVIER 1990
C MODIFS : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C--------------------------------------------------------------------012
      INTEGER NSEFJ(8),    TE1SOB(4,6), PLUTET, KONFIG(8), MODECO(3)
      INTEGER FACE(4,3),   TETASS(24)
      INTEGER CNFPYR(5,5), CNFPEN(6,6), CNFHEX(8,8)
C
      DATA (CNFPYR(I,1),I=1,5) /1,2,3,4, 5/
      DATA (CNFPYR(I,2),I=1,5) /2,3,4,1, 5/
      DATA (CNFPYR(I,3),I=1,5) /3,4,1,2, 5/
      DATA (CNFPYR(I,4),I=1,5) /4,1,2,3, 5/
      DATA (CNFPYR(I,5),I=1,5) /1,2,3,4, 5/
C     ATTENTION: LE 5 EME SOMMET NE PEUT ETRE PERMUTE POUR UNE PYRAMIDE!
C
      DATA (CNFPEN(I,1),I=1,6) /1,2,3, 4,5,6/
      DATA (CNFPEN(I,2),I=1,6) /2,3,1, 5,6,4/
      DATA (CNFPEN(I,3),I=1,6) /3,1,2, 6,4,5/
      DATA (CNFPEN(I,4),I=1,6) /4,6,5, 1,3,2/
      DATA (CNFPEN(I,5),I=1,6) /5,4,6, 2,1,3/
      DATA (CNFPEN(I,6),I=1,6) /6,5,4, 3,2,1/
C
      DATA (CNFHEX(I,1),I=1,8) /1,2,3,4, 5,6,7,8/
      DATA (CNFHEX(I,2),I=1,8) /2,3,4,1, 6,7,8,5/
      DATA (CNFHEX(I,3),I=1,8) /3,4,1,2, 7,8,5,6/
      DATA (CNFHEX(I,4),I=1,8) /4,1,2,3, 8,5,6,7/
      DATA (CNFHEX(I,5),I=1,8) /5,8,7,6, 1,4,3,2/
      DATA (CNFHEX(I,6),I=1,8) /6,5,8,7, 2,1,4,3/
      DATA (CNFHEX(I,7),I=1,8) /7,6,5,8, 3,2,1,4/
      DATA (CNFHEX(I,8),I=1,8) /8,7,6,5, 4,3,2,1/
C
       IF( NCOGEL .EQ. 9 ) THEN
C         PYRAMIDE
          IQUAD = 1
          NBS   = 5
       ELSE IF( NCOGEL .EQ. 6 ) THEN
C         PENTAEDRE
          IQUAD = 1
          NBS   = 6
       ELSE
C         HEXAEDRE
          IQUAD = 3
          NBS   = 8
       ENDIF
C
C ---- RECHERCHE DU PLUS PETIT SOMMET DE L'EF ----
       INF    = NSEFJ(1)
       INDINF = 1
       DO 20 I=2,NBS
          IX=NSEFJ(I)
          IF( IX .LT. INF ) THEN
             INF    = IX
             INDINF = I
c
             if( ncogel .eq. 9 .and. indinf .eq. 5 ) then
                print *,'Attention pyramide de sommet min en 5'
                print *,'Pyramide:',(nsefj(k),k=1,nbs)
             endif
c
          ENDIF
 20    CONTINUE
C
C ---- KONFIG=No des Sommets le plus petit en tete -----
       KONFIG(6)=0
       KONFIG(7)=0
       KONFIG(8)=0
       IF( NCOGEL .EQ. 9 ) THEN
C         PYRAMIDE
          DO 25 I=1,5
             KONFIG(I)=NSEFJ( CNFPYR(I,INDINF) )
 25       CONTINUE
       ELSE IF( NCOGEL .EQ. 6 ) THEN
C         PENTAEDRE
          DO 30 I=1,6
             KONFIG(I)=NSEFJ( CNFPEN(I,INDINF) )
 30       CONTINUE
       ELSE
C         HEXAEDRE
          DO 40 I=1,8
             KONFIG(I)=NSEFJ( CNFHEX(I,INDINF) )
 40       CONTINUE
       ENDIF
C
C----- DECOUPE DES 1 ou 3 FACES RESTANTES SANS LE SOMMET MIN -----
       CALL NOFAPH( NCOGEL, KONFIG, FACE )
       DO 50 I=2,3
          MODECO(I)=0
 50    CONTINUE
C      RECHERCHE DU SOMMET MIN (IDINF) DE LA FACE SANS SOMMET MIN DE L'EF
       DO 60 I=1,IQUAD
          INFFCE=FACE(1,I)
          IDINF=1
          DO 70 J=2,4
             IF ( INFFCE .GT. FACE(J,I) ) THEN
                INFFCE = FACE(J,I)
                IDINF = J
             ENDIF
 70       CONTINUE
          IF( IDINF.EQ.2 .OR. IDINF.EQ.4 ) THEN
C            DECOUPE SELON LA DIAGONALE 24
             MODECO(I)=2
          ELSE
C            DECOUPE SELON LA DIAGONALE 13
             MODECO(I)=1
          ENDIF
 60    CONTINUE
C
C----- DECOUPE DE L'EF EN TETRA A PARTIR DE LA COUPE DE SES FACES --
C----- ET RANGEMENT DU RESULTAT DANS TE1SOB
       CALL DECOBJ( NCOGEL, MODECO, TETASS, PLUTET )
       K=0
       DO 90 I=1,PLUTET
          DO 80 J=1,4
             K=K+1
             TE1SOB(J,I)=KONFIG(TETASS(K))
 80       CONTINUE
 90    CONTINUE
       RETURN
       END
