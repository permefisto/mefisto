      SUBROUTINE SE30ISA6( QUALMN, MXSOM,  NBSOM,  XYZSOM,
     %                     NT0,    NBTRIA, NOTRIA6, MODIF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     SI LE TRIANGLE NT0 A UNE TRES PETITE ARETE ALORS
C -----     SES 2 SOMMETS SONT IDENTIFIES ET NT0 EST DETRUIT
C           AINSI QUE SON TRIANGLE OPPOSE PAR CETTE ARETE

C ENTREES:
C --------
C QUALMN : QUALITE MINIMALE AU DESSOUS DE LAQUELLE LA QUALITE DOIT
C          ETRE AMELIOREE
C NT0    : NUMERO NOTRIA6 DU TRIANGLE D'ANGLE TROP GRAND

C MODIFIES:
C ---------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE   AVANT et APRES
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE AVANT et APRES
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS AVANT et APRES
C NOTRIA6 : NUMERO DES 3 SOMMETS ET 3 TRIANGLES ADJACENTS PAR LES ARETES

C SORTIES:
C --------
C NBSOM  : NOMBRE DE SOMMETS   DU MAILLAGE
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE
C MODIF  : 2 SI 2 SOMMETS D'UNE PETITE ARETE ONT ETE IDENTIFIES
C          0 SI LE TRIANGLE NT0 N'A PAS ETE MODIFIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET  Saint Pierre du Perray             Novembre 2019
C2345X7..............................................................012
      PARAMETER        (MXTRIT=128)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0
      CHARACTER*48      KNM

      INTEGER           NOTRIA6(6,*), LITRIT(MXTRIT)
      REAL              XYZSOM(3,MXSOM), XYZ(3,3), LONARE(3)

      TRACTE0 = TRACTE

C     CALCUL DE LA QUALITE Q0 DU TRIANGLE NT0
      CALL QUATRI( NOTRIA6(1,NT0), XYZSOM, Q0 )
      IF( Q0 .GE. QUALMN ) GOTO 9999

      MODIF = 0
      NT12  = 0
      NT31  = 0
      NT24  = 0
      NT43  = 0

C     CONVERSION RADIANS DEGRES
      RADEGR = 45. / ATAN(1.)

C     TRACE DU TRIANGLE NT0 ET DES 3 TRIANGLES ADJACENTS
      LITRIT(1) = NT0
      KNM = 'Debut se30isa6: Triangle                     '
      WRITE(KNM(27:33),'(I7)') NT0
      CALL SANSDBL( KNM, NC )
      CALL TRTRIAN( KNM(1:NC), XYZSOM, 6, MXTRIT,
     %              1, LITRIT, NOTRIA6 )


C     LE TRIANGLE NT0 EST DE MAUVAISE QUALITE
C     CALCUL DES ANGLES MIN MAX ET PERIMETRE DU TRIANGLE NT0
C     ------------------------------------------------------
      DO J=1,3
C        LES 3 COORDONNEES DU SOMMET J DU TRIANGLE NT0
         NS = NOTRIA6( J, NT0 )
         XYZ( 1, J ) = XYZSOM( 1, NS )
         XYZ( 2, J ) = XYZSOM( 2, NS )
         XYZ( 3, J ) = XYZSOM( 3, NS )
      ENDDO

C     LA LONGUEUR DES 3 COTES DU TRIANGLE NT0
      JMIN    = 0
      ARETMIN = 1E28
      PERIME  = 0
      DO J=1,3
         IF( J .EQ. 3 ) THEN
            J1=1
         ELSE
            J1 = J+1
         ENDIF
         LONARE(J) = SQRT( ( XYZ(1,J)-XYZ(1,J1) ) ** 2
     %                   + ( XYZ(2,J)-XYZ(2,J1) ) ** 2
     %                   + ( XYZ(3,J)-XYZ(3,J1) ) ** 2 )

         PERIME = PERIME + LONARE(J)
         IF( LONARE(J) .LT. ARETMIN ) THEN
            ARETMIN = LONARE(J)
            JMIN    = J
         ENDIF

      ENDDO

      ARPEMIN = ARETMIN / PERIME

ccc   IF( ARPEMIN .LE. 0.015 ) THEN    IDENTIFIE PAS ASSEZ
ccc   IF( ARPEMIN .LE. 0.01  ) THEN    IDENTIFIE PAS ASSEZ
ccc   IF( ARPEMIN .LE. 0.03  ) THEN    IDENTIFIE TROP
ccc   IF( ARPEMIN .LE. 0.026 ) THEN
      IF( ARPEMIN .LE. 0.028 ) THEN

C        L'ARETE JMIN PEUT ELLE ETRE REDUITE A UN SOMMET PAR
C        IDENTIFICATION DE SES 2 SOMMETS EXTREMITES?
C        ===================================================
         NS1 = NOTRIA6( JMIN, NT0 )
         IF( JMIN .EQ. 3 ) THEN
            J1=1
         ELSE
            J1 = JMIN+1
         ENDIF
         NS2 = NOTRIA6( J1, NT0 )

C        RECHERCHE DE 2 TRIANGLES OPPOSES A NT0 AYANT UNE ARETE COMMUNE
C        SI LES 2 SOMMETS NS2 NS1 SONT IDENTIFIES ALORS
C        LE TETRAEDRE AINSI FORME DEVIENT UN DOUBLE TRIANGLE
C        CE QUI REND LA SURFACE NON FERMEE. DONC A INTERDIRE
         CALL TR2ADJ1AR( NT0, NOTRIA6, NONOUI )
         IF( NONOUI .NE. 0 ) THEN
C           IL EXISTE UNE ARETE COMMUNE A 2 TRIANGLES ADJACENTS A NT0
C           => PAS D'IDENTIFICATION DE NS1 et NS2
C              DE PLUS, COMME L'EPAISSEUR DU MAILLAGE
C              EST FAIBLE NT0 N'EST PAS MODIFIE
            GOTO 9990
         ENDIF

C        DE MEME POUR LE TRIANGLE AJACENT NT1 PAR L'ARETE JMIN DE NT0
         NT1 = NOTRIA6( 3+JMIN, NT0 )
         CALL TR2ADJ1AR( NT1, NOTRIA6, NONOUI )
         IF( NONOUI .NE. 0 ) THEN
C           IL EXISTE UNE ARETE COMMUNE A 2 TRIANGLES ADJACENTS A NT1
C           => PAS D'IDENTIFICATION DE NS1 ET NS2
C              DE PLUS, COMME L'EPAISSEUR DU MAILLAGE
C              EST FAIBLE NT0 ET NT1 NE SONT PAS MODIFIES
            GOTO 9990
         ENDIF

C        OUI: PAS DE RISQUE D'ECRASEMENT. IDENTIFICATION NS2->NS1
C        --------------------------------------------------------
         PRINT*,'se30isa6: le SOMMET',NS2,' DEVIENT le SOMMET',NS1,
     %          ' les 3 ARETES:',LONARE,
     %          ' ARETE MIN=',ARETMIN,' TRES PETITE. PERIMETRE=',PERIME,
     %          ' ARETE Min/PERIMETRE=',ARPEMIN

         LITRIT(1)=NT0
         KNM = 'se30isa6 ident 2 st:                           '
         WRITE(KNM(23:29),'(I7)') NS2
         WRITE(KNM(31:37),'(I7)') NS1
         CALL SANSDBL( KNM, NC )
         CALL TRTRIAN( KNM(1:NC), XYZSOM, 6, MXTRIT,
     %                 1, LITRIT, NOTRIA6 )

         DO N=1,NBTRIA
            DO J=1,3
               IF( NOTRIA6(J,N) .EQ. NS2 ) THEN
C                 IDENTIFICATION DE NS2 A NS1
                  NOTRIA6(J,N) = NS1
               ENDIF
            ENDDO
         ENDDO

C        RECHERCHE ET SUPPRESSION DES TRIANGLES AYANT 2 MEMES SOMMETS
C        NSMIN APRES L'IDENTIFICATION NS2 -> NS1
C        MISE A JOUR DU TABLEAU XYZSOM ET NOSOEF EN RENUMEROTANT
C        LES SOMMETS ET ELIMINANT LES EF DESACTIVES
C        NEWS EST UN TABLEAU AUXILIAIRE
         MNNEWS = 0
         MXS    = 1+NBSOM
ccc         CALL TNMCDC( 'ENTIER', MXS, MNNEWS )
ccc pb         CALL MAJXYZNSE( NBSOM, XYZSOM, MCN(MNNEWS),
ccc pb     %                   3,  6, NBTRIA, NOTRIA6 )
ccc         CALL TNMCDS( 'ENTIER', MXS, MNNEWS )

         MODIF = 2
         GOTO 9999

      ENDIF


C     PAS DE MODIFICATION DU TRIANGLE NT0
C     -----------------------------------
 9990 MODIF = 0

      PRINT*,'se30isa6: TRIANGLE',NT0,
     %       ' St:',(NOTRIA6(K,NT0),K=1,3),
     %       ' Qualite=',Q0,' PAS de d''IDENTIFICATION de 2 SOMMETS'

 9999 TRACTE = TRACTE0
      RETURN
      END
