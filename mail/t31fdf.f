      SUBROUTINE T31FDF( NMOBJT, NUOBJT, MNTSMA, MNSOMM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES FIL DE FER DES FACES DE LA SURFACE
C -----
C ENTREE :
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOMM : ADRESSE MCN DU TABLEAU 'XYZSOMMET'    A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C ......................................................................
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*(*)     NMOBJT
      CHARACTER*8       NMSOMM
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOSOEL(1:64)
      REAL              XYZP(3), XYZTG1(3), XYZTG2(3)

      IF( MNTSMA.LE.0 .OR. MNSOMM.LE.0 ) RETURN

C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS   = MNSOMM + WYZSOM - 3
      NBSOM = MCN( MNSOMM + WNBSOM )
      MNTG  = MNS + 3 * NBSOM

C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNTSMA),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN

C     LA BOUCLE SUR LES EF DU MAILLAGE
C     --------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      MNSTS = MNSOMM + WYZSOM
      DO 100 N=1,NBEFOB

C        LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNTSMA, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( NOSOEL(1) .EQ. 0 ) GOTO 100

C        LE NOMBRE DE SOMMETS DE CET ELEMENT EST NCOGEL=3 OU 4
C        LE TRACE DES ARETES DE LA FACE
         IF( NUEFTG .LE. 0 ) THEN

C           EF DROIT
            MN2 = MNS + 3 * NOSOEL(NCOGEL)
            DO I=1,NCOGEL
               MN1 = MN2
               MN2 = MNS + 3 * NOSOEL(I)
               CALL TRAIT3D( NCOUAF, RMCN(MN1), RMCN(MN2) )
            ENDDO

         ELSE

C           EF A TG
            DO 20 I=1,NCOGEL
               IF( I .NE. NCOGEL ) THEN
                  I1 = I + 1
               ELSE
                  I1 = 1
               ENDIF
C              ADRESSE DES 2 SOMMETS DE L'ARETE I
               MN1 = MNS + 3 * NOSOEL(I)
               MN2 = MNS + 3 * NOSOEL(I1)

C              LE NUMERO DE LA PREMIERE TANGENTE DE L'ARETE AU SOMMET J
               NTG = NOSOEL( 3 + 2 * I )
               IF( NTG .NE. 0 ) THEN
C                 IL EXISTE UNE TANGENTE
                  MN = MNTG + 3 * ABS(NTG)
                  IF( NTG .GT. 0 ) THEN
                      XYZTG1(1) = RMCN( MN     )
                      XYZTG1(2) = RMCN( MN + 1 )
                      XYZTG1(3) = RMCN( MN + 2 )
                  ELSE
                      XYZTG1(1) = -RMCN( MN     )
                      XYZTG1(2) = -RMCN( MN + 1 )
                      XYZTG1(3) = -RMCN( MN + 2 )
                  ENDIF
               ELSE
C                 PAS DE TANGENTE: COTE DROIT
                  XYZTG1(1) = RMCN( MN2     ) - RMCN( MN1     )
                  XYZTG1(2) = RMCN( MN2 + 1 ) - RMCN( MN1 + 1 )
                  XYZTG1(3) = RMCN( MN2 + 2 ) - RMCN( MN1 + 2 )
               ENDIF

C              LE NUMERO DE LA SECONDE TANGENTE DU SOMMET J
               NTG = NOSOEL( 4 + 2 * I1 )
               IF( NTG .NE. 0 ) THEN
C                 IL EXISTE UNE TANGENTE
                  MN = MNTG + 3 * ABS(NTG)
                  IF( NTG .GT. 0 ) THEN
                      XYZTG2(1) = RMCN( MN     )
                      XYZTG2(2) = RMCN( MN + 1 )
                      XYZTG2(3) = RMCN( MN + 2 )
                  ELSE
                      XYZTG2(1) = -RMCN( MN     )
                      XYZTG2(2) = -RMCN( MN + 1 )
                      XYZTG2(3) = -RMCN( MN + 2 )
                  ENDIF
               ELSE
C                 PAS DE TANGENTE: COTE DROIT
                  XYZTG2(1) = RMCN( MN1     ) - RMCN( MN2     )
                  XYZTG2(2) = RMCN( MN1 + 1 ) - RMCN( MN2 + 1 )
                  XYZTG2(3) = RMCN( MN1 + 2 ) - RMCN( MN2 + 2 )
               ENDIF
               CALL TRAR3D( NCOUAF, PREDUA, RMCN(MN1), RMCN(MN2),
     %                      XYZTG1, XYZTG2 )
 20         ENDDO

         ENDIF

C        TRACE EVENTUEL DU NO DE L'EF
         IF( IAVNEF .NE. 0 ) THEN
            XYZP(1) = 0
            XYZP(2) = 0
            XYZP(3) = 0
            DO J=1,NCOGEL
               MN1 = MNS + 3 * NOSOEL(J)
               XYZP(1) = XYZP(1) + RMCN( MN1   )
               XYZP(2) = XYZP(2) + RMCN( MN1+1 )
               XYZP(3) = XYZP(3) + RMCN( MN1+2 )
            ENDDO
            XYZP(1) = XYZP(1) / NCOGEL
            XYZP(2) = XYZP(2) / NCOGEL
            XYZP(3) = XYZP(3) / NCOGEL
            WRITE( NMSOMM, '(I8)' ) N
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, XYZP, NMSOMM(1:L) )
         ENDIF

C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            DO J=1,NCOGEL
               CALL TRST3D( NCONSO, NOSOEL(J), RMCN(MNSTS) )
ccc               MN = MNS + 3 * NOSOEL(J)
ccc               WRITE( NMSOMM, '(I8)' ) NOSOEL(J)
ccc               CALL SANSBL( NMSOMM, L )
ccc               CALL TEXTE3D( NCONSO, RMCN(MN), NMSOMM(1:L) )
            ENDDO
         ENDIF
 100  ENDDO

C     TRACE DE LA POIGNEE ET DU NOM DE LA SURFACE
C     LES COORDONNEES DE LA POIGNEE
      XYZP(1) = 0
      XYZP(2) = 0
      XYZP(3) = 0
      DO J=1,NCOGEL
         MN1 = MNS + 3 * NOSOEL(J)
         XYZP(1) = XYZP(1) + RMCN( MN1   )
         XYZP(2) = XYZP(2) + RMCN( MN1+1 )
         XYZP(3) = XYZP(3) + RMCN( MN1+2 )
      ENDDO
      XYZP(1) = XYZP(1) / NCOGEL
      XYZP(2) = XYZP(2) / NCOGEL
      XYZP(3) = XYZP(3) / NCOGEL
      CALL ITEMS3( XYZP, NMOBJT, NUOBJT )

      RETURN
      END
