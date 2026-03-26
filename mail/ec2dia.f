        SUBROUTINE EC2DIA( COSMAXPL, NT1, I, NOTRIA, NBSOM, XYZSOM,
     %                     NT2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER LE TRIANGLE ADJACENT A L'ARETE I DE NT1
C -----    REDECOUPER LE QUADRANGLE AINSI FORME SI CES 2 TRIANGLES
C          SONT COPLANAIRES ET DE MEILLEURES QUALITES

C          ATTENTION: DANS R3 LES 4 SOMMETS DES 2 TRIANGLES
C                     FORMENT UN TETRAEDRE QUI DISPARAIT OU
C                     APPARAIT SELON LE CHANGEMENT DE DIAGONALE

C ENTREES:
C --------
C COSMAXPL: COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 PLANS AU
C           DESSOUS DUQUEL LES 2 TRIANGLES SONT JUGES NON COPLANAIRES
C NT1     : NUMERO DANS NOTRIA DU TRIANGLE A TRAITER
C I       : NUMERO DE L'ARETE A TRAITER
C NBSOM   : NOMBRE DE SOMMETS DU MAILLAGE

C MODIFIES:
C ---------
C NOTRIA  : NUMEROS DES 3 SOMMETS ET TRIANGLES ADJACENTS
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                                  ADJACENT PAR L'ARETE i
C XYZSOM : 3 COORDONNEES DES SOMMETS

C SORTIES:
C --------
C NT2    : >0 NUMERO DU TRIANGLE ADJACENT SI L'ECHANGE A LIEU
C          =0 SI PAS D'ECHANGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    OCTOBRE 1993
C MODIFS : ALAIN PERRONNET  Saint Pierre du Perray          Mars    2020
C....................................................................012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL              XYZSOM(3,NBSOM), PTPROJ(3)
      DOUBLE PRECISION  CBPTTR(3)
      INTEGER           NOTRIA(6,*)
      INTEGER           NOSOTR(3)

      IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 9000

C     LE TRIANGLE ADJACENT A L'ARETE I
      NT2 = NOTRIA( 3+I, NT1 )
      IF( NT2 .GT. 0 ) THEN

         IF( NOTRIA(1,NT2) .LE. 0 ) GOTO 9000


ccc         if( ( nt1 .eq. 1829 .and. nt2 .eq. 1870 ) .OR.
ccc     %       ( nt2 .eq. 1829 .and. nt1 .eq. 1870 ) ) THEN
ccc            print*,'ec2dia: echange des triangles',nt1,nt2
ccc            print*,'ec2dia: notria(',nt1,')=',(notria(k,nt1),k=1,6)
ccc            print*,'ec2dia: notria(',nt2,')=',(notria(k,nt2),k=1,6)
ccc            print*
ccc         endif


C        L'ARETE I EST INTERNE PARTAGEE PAR NT1 ET NT2
         NS1 = NOTRIA( I, NT1 )
         IF( I .LT. 3 ) THEN
            I1 = I + 1
         ELSE
            I1 = 1
         ENDIF

         NS2 = NOTRIA( I1, NT1 )
         IF( I .GT. 1 ) THEN
            I2 = I - 1
         ELSE
            I2 = 3
         ENDIF

C        LE 3-EME SOMMET DU TRIANGLE NT1
         NS3 = NOTRIA( I2, NT1 )

C        LE NUMERO DES ARETES DANS NT2
         IF( NOTRIA(1,NT2) .EQ. NS2 .AND.
     %       NOTRIA(2,NT2) .EQ. NS1 ) THEN
            NA  = 1
            NA1 = 2
            NA2 = 3
            NS4 = NOTRIA(NA2,NT2)
            NSENS = 1
         ELSE IF( NOTRIA(2,NT2) .EQ. NS2 .AND.
     %            NOTRIA(3,NT2) .EQ. NS1 ) THEN
            NA  = 2
            NA1 = 3
            NA2 = 1
            NS4 = NOTRIA(NA2,NT2)
            NSENS = 1
         ELSE IF( NOTRIA(3,NT2) .EQ. NS2 .AND.
     %            NOTRIA(1,NT2) .EQ. NS1 ) THEN
            NA  = 3
            NA1 = 1
            NA2 = 2
            NS4 = NOTRIA(NA2,NT2)
            NSENS = 1
         ELSE IF( NOTRIA(1,NT2) .EQ. NS1 .AND.
     %            NOTRIA(2,NT2) .EQ. NS2 ) THEN
            NA  = 1
            NA1 = 3
            NA2 = 2
            NS4 = NOTRIA(NA1,NT2)
            NSENS = -1
         ELSE IF( NOTRIA(2,NT2) .EQ. NS1 .AND.
     %            NOTRIA(3,NT2) .EQ. NS2 ) THEN
            NA  = 2
            NA1 = 1
            NA2 = 3
            NS4 = NOTRIA(NA1,NT2)
            NSENS = -1
         ELSE IF( NOTRIA(3,NT2) .EQ. NS1 .AND.
     %            NOTRIA(1,NT2) .EQ. NS2 ) THEN
            NA  = 3
            NA1 = 2
            NA2 = 1
            NS4 = NOTRIA(NA1,NT2)
            NSENS = -1
         ELSE
            WRITE(IMPRIM,*) 'ec2dia: ANOMALIE A TRAITER'
            GOTO 9000
         ENDIF

C        LE COSINUS DE L'ANGLE DES NORMALES AUX TRIANGLES NT1 NT2
         CALL COS2TR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %                XYZSOM(1,NS2), XYZSOM(1,NS1), XYZSOM(1,NS4),
     %                COS2PL, IERR1, IERR2 )

         IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) GOTO 20
C        SI UN TRIANGLE EST SANS NORMALE
C        ALORS L'ECHANGE EST FORCE INDEPENDAMMENT DE COSMAXPL

         IF( COS2PL .LT. COSMAXPL ) THEN
C           ANGLE DIEDRE DES 2 TRIANGLES TROP GRAND
            GOTO 9000
         ENDIF

C        LES 2 TRIANGLES FORMENT ILS UN QUADRANGLE 1423 CONVEXE?
C        -------------------------------------------------------
         CALL QUADCXR( XYZSOM(1,NS1), XYZSOM(1,NS4),
     %                 XYZSOM(1,NS2), XYZSOM(1,NS3),  NONOUI )
         IF( NONOUI .EQ. 0 ) GOTO 9000

C        LE POINT PROJECTION DE NS4 SUR LE PLAN DU TRIANGLE NT1
C        EST IL INTERNE AU TRIANGLE NT1?
C        -------------------------------------------------------
C        POINT DE PROJECTION DE NS4 SUR LE PLAN DU TRIANGLE NT1
         CALL PRPTPL( XYZSOM(1,NS4),
     %                XYZSOM(1,NOTRIA(1,NT1)),
     %                XYZSOM(1,NOTRIA(2,NT1)),
     %                XYZSOM(1,NOTRIA(3,NT1)),
     %                PTPROJ, IER )
C        PTPROJ(3) : XYZ DU POINT PROJETE SUR LE PLAN de NT1
         IF( IER .NE. 0 ) GOTO 9000

C        3 COORDONNEES BARYCENTRIQUES DU POINT PTPROJ DANS LE TRIANGLE NT1
         CALL CBPTTRR( XYZSOM(1,NOTRIA(1,NT1)),
     %                 XYZSOM(1,NOTRIA(2,NT1)),
     %                 XYZSOM(1,NOTRIA(3,NT1)),
     %                 PTPROJ,  CBPTTR )

C        PTPROJ EST IL INTERNE OU PROCHE DU TRIANGLE NT1?
         IF( ABS(CBPTTR(1))+ABS(CBPTTR(2))+ABS(CBPTTR(3)) .LE. 1.01D0 )
     %      THEN
C           OUI: PAS D'ECHANGE
            GOTO 9000
         ENDIF

C        LE POINT PROJECTION DE NS3 SUR LE PLAN DU TRIANGLE NT2
C        EST IL INTERNE AU TRIANGLE NT2?
C        -------------------------------------------------------
C        POINT DE PROJECTION DE NS4 SUR LE PLAN DU TRIANGLE NT2
         CALL PRPTPL( XYZSOM(1,NS3),
     %                XYZSOM(1,NOTRIA(1,NT2)),
     %                XYZSOM(1,NOTRIA(2,NT2)),
     %                XYZSOM(1,NOTRIA(3,NT2)),
     %                PTPROJ, IER )
C        PTPROJ(3) : XYZ DU POINT PROJETE SUR LE PLAN de NT2
         IF( IER .NE. 0 ) GOTO 9000

C        3 COORDONNEES BARYCENTRIQUES DU POINT PTPROJ DANS LE TRIANGLE NT2
         CALL CBPTTRR( XYZSOM(1,NOTRIA(1,NT2)),
     %                 XYZSOM(1,NOTRIA(2,NT2)),
     %                 XYZSOM(1,NOTRIA(3,NT2)),
     %                 PTPROJ,  CBPTTR )

C        PTPROJ EST IL INTERNE OU PROCHE DU TRIANGLE NT2?
         IF( ABS(CBPTTR(1))+ABS(CBPTTR(2))+ABS(CBPTTR(3)) .LE. 1.01D0 )
     %      THEN  
C           OUI: PAS D'ECHANGE
            GOTO 9000
         ENDIF

C        LES 2 FACES NT1 NT2 SONT COPLANAIRES OU PRESQUE
C        QUALITE DES 2 CHOIX POSSIBLES DE TRIANGULATIONS
C        -----------------------------------------------
C        TRIANGLE 123
         CALL QUATRI( NOTRIA(1,NT1), XYZSOM, Q12 )

C        TRIANGLE 214
         CALL QUATRI( NOTRIA(1,NT2), XYZSOM, Q   )
         Q12 = MIN( Q12, Q )

C        TRIANGLE 342
         NOSOTR(1) = NS3
         NOSOTR(2) = NS4
         NOSOTR(3) = NS2
         CALL QUATRI( NOSOTR, XYZSOM, Q34 )
C        COMPARAISON DES QUALITES  MIN( Q123, Q214 ) ET Q342
         IF( Q34 .LE. Q12 ) GOTO 9000

C        TRIANGLE 431
         NOSOTR(1) = NS4
         NOSOTR(2) = NS3
         NOSOTR(3) = NS1
         CALL QUATRI( NOSOTR, XYZSOM, Q )
C        COMPARAISON DES QUALITES  MIN( Q123, Q214 ) ET Q431
         IF( Q .LE. Q12 ) GOTO 9000
         Q34 = MIN( Q34, Q )

C        COMPARAISON DU RAPPORT DES SURFACES DES 2 COUPLES DE TRIANGLES
C        --------------------------------------------------------------
         ST1 = SURTRR( XYZSOM( 1, NOTRIA(1,NT1) ),
     %                 XYZSOM( 1, NOTRIA(2,NT1) ),
     %                 XYZSOM( 1, NOTRIA(3,NT1) ) )

         ST2 = SURTRR( XYZSOM( 1, NOTRIA(1,NT2) ),
     %                 XYZSOM( 1, NOTRIA(2,NT2) ),
     %                 XYZSOM( 1, NOTRIA(3,NT2) ) )
C        LE RAPPORT DES SURFACES ENTRE 0 et 1
         IF( ST1 .LE. ST2 ) THEN
            RAPS12 = ST1 / ST2
         ELSE
            RAPS12 = ST2 / ST1
         ENDIF

         ST342 = SURTRR( XYZSOM( 1, NS3 ),
     %                   XYZSOM( 1, NS4 ),
     %                   XYZSOM( 1, NS2 ) )

         ST431 = SURTRR( XYZSOM( 1, NS4 ),
     %                   XYZSOM( 1, NS3 ),
     %                   XYZSOM( 1, NS1 ) )
         IF( ST342 .LE. ST431 ) THEN
            RAPS34 = ST342 / ST431
         ELSE
            RAPS34 = ST431 / ST342
         ENDIF

ccc      IF( RAPS34 .LT. RAPS12 .AND. RAPS34 .LT. 0.055 ) THEN
         IF( RAPS34 .LT. 0.0667 ) THEN
C           L'ECHANGE ARETE 12 PAR L'ARETE 34 ENTRAINERAIT UN TRIANGLE
C           DE TRES PETITE SURFACE PAR RAPPORT AU COUPLE NT1 NT2 EXISTANT
            GOTO 9000
         ENDIF

C        ECHANGE EFFECTIF DE LA DIAGONALE DES 2 TRIANGLES NT1 NT2
C        ========================================================
C        ICI NT1 NT2 DOIVENT ETRE REMPLACES SUR EUX MEMES
 20      NOTRIA(I, NT1) = NS4
         IF( NSENS .GT. 0 ) THEN
            NOTRIA(NA ,NT2) = NS3
         ELSE
            NOTRIA(NA2,NT2) = NS3
         ENDIF
C
C        LES TRIANGLES ADJACENTS AVANT ECHANGE
         NTT1 = NOTRIA(I1+3, NT1)
         NTT2 = NOTRIA(I2+3, NT1)
         NTT3 = NOTRIA(NA1+3,NT2)
         NTT4 = NOTRIA(NA2+3,NT2)

C        NOUVELLE REPARTITION DES TRIANGLES ADJACENTS PAR LES ARETES
         NOTRIA(I +3 ,NT1) = NTT4
         NOTRIA(I1+3 ,NT1) = NTT1
         NOTRIA(I2+3 ,NT1) = NT2

         NOTRIA(NA +3,NT2) = NTT2
         NOTRIA(NA1+3,NT2) = NTT3
         NOTRIA(NA2+3,NT2) = NT1

C        L'INVERSE POUR LES TRIANGLES NTT2 ET NTT4
C        LE NUMERO D'ARETE DE NS1-NS2 DANS NT2
         IF( (NOTRIA(1,NTT4) .EQ. NS2 .AND.
     %        NOTRIA(2,NTT4) .EQ. NS4) .OR.
     %       (NOTRIA(1,NTT4) .EQ. NS4 .AND.
     %        NOTRIA(2,NTT4) .EQ. NS2) ) THEN
            NA = 4
         ELSE IF( (NOTRIA(2,NTT4) .EQ. NS2 .AND.
     %             NOTRIA(3,NTT4) .EQ. NS4) .OR.
     %            (NOTRIA(2,NTT4) .EQ. NS4 .AND.
     %             NOTRIA(3,NTT4) .EQ. NS2) ) THEN
            NA = 5
         ELSE
            NA = 6
         ENDIF
         NOTRIA(NA,NTT4) = NT1

         IF( (NOTRIA(1,NTT2) .EQ. NS1 .AND.
     %        NOTRIA(2,NTT2) .EQ. NS3) .OR.
     %       (NOTRIA(1,NTT2) .EQ. NS3 .AND.
     %        NOTRIA(2,NTT2) .EQ. NS1) ) THEN
            NA = 4
         ELSE IF( (NOTRIA(2,NTT2) .EQ. NS1 .AND.
     %             NOTRIA(3,NTT2) .EQ. NS3) .OR.
     %            (NOTRIA(2,NTT2) .EQ. NS3 .AND.
     %             NOTRIA(3,NTT2) .EQ. NS1) ) THEN
            NA = 5
         ELSE
            NA = 6
         ENDIF
         NOTRIA(NA,NTT2) = NT2

         print*
         print*,'ec2dia: Echange diagonale des 2 TRIANGLES FINAUX avec Q
     %12=',Q12,' Q34=',Q34,' RAPS12=',RAPS12,' RAPS34=',RAPS34,
     %' COS2PL=',COS2PL
         print*,'ec2dia: Triangle',NT1,' St:',(NOTRIA(k,NT1),k=1,6)
         print*,'ec2dia: Triangle',NT2,' St:',(NOTRIA(k,NT2),k=1,6)

         GOTO 9999
      ENDIF

C     PAS D'ECHANGE DES 2 DIAGONALES
C     ==============================
 9000 NT2 = 0

 9999 RETURN
      END
