      SUBROUTINE GEFAVO( NBSOEF , NBCUVO , NOSOCU ,
     &                   MXMOFA , L1FACE , L2FACE , MNLFAC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FORMER PAR HACHAGE LE TABLEAU DES FACES D'UN VOLUME
C -----
C ENTREES:
C --------
C NBSOEF : NOMBRE DE SOMMETS ET TANGENTES PAR CUBE (8 SANS TG ET 32 AVEC TG)
C NBCUVO : NOMBRE DE CUBES DU VOLUME
C NOSOCU : LES 8 NUMEROS DES SOMMETS DES NBCUVO CUBES DU VOLUME
C MXMOFA : NOMBRE DE MOTS EN PLUS PAR FACE
C          2 PAR EXEMPLE POUR LE NUMERO DES 2 CUBES CONTENANT LA FACE
C          0 SI PAS D'UN TEL STOCKAGE
C
C SORTIES:
C --------
C L1FACE : NOMBRE DE MOTS PAR FACE DU TABLEAU LFACES
C L2FACE : NOMBRE DE FACES DU TABLEAU LFACES
C MNLFAC : ADRESSE DANS M DU TABLEAU LFACES DES FACES DU MAILLAGE
C          EN SORTIE LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C                    LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C                    LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C                    LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                                0 SI TRIANGLE
C                    LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C                    LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                                0 SI PAS DE 1-ER  CUBE
C                    LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                                0 SI PAS DE 2-EME CUBE
C                    LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE
C                                TABLEAU NUTGFA
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1988
C ......................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/ponoel.inc"
      INTEGER           NOSOCU(NBSOEF,NBCUVO)
      INTEGER           NGS(4),NGS1(4)
C
C     CALCUL DU NOMBRE DES FACES DU MAILLAGE
C     ---------------------------------------
C     MXMOFA : ICI NOMBRE MAXIMAL DE CUBES CONTENANT UNE MEME FACE
      L1FACE = 5 + MXMOFA
C
C     NOMBRE DE FACES : FORMULE A MODIFIER EVENTUELLEMENT
C     ---------------------------------------------------
      L2FACE = NBCUVO * 6
C
C     ADRESSAGE DES TABLEAUX KFACE ET DU POINTEUR DES FACES REFEREES
C     ---------------------------------------------------------------
      L1     = L1FACE * L2FACE
      MNLFAC = 0
      CALL TNMCDC( 'ENTIER' , L1 , MNLFAC )
C
C     LE TABLEAU DES FACES EST INITIALISE A ZERO
      CALL AZEROI( L1 , MCN(MNLFAC) )
C
C     LA 1-ERE FACE LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
      LIBREF = L2FACE
C
C     FORMATION DU TABLEAU DES NO DES SOMMETS DES FACES
C     ==================================================
      DO 100 N=1,NBCUVO
C
C        LE NOMBRE DE SOMMETS DU CUBE
         DO NBSOMM=8,1,-1
            IF( NOSOCU(NBSOMM,N) .GT. 0 ) GOTO 20
         ENDDO
C
C        LE NUMERO DE TYPE DE L'ELEMENT DROIT
 20      IF( NBSOMM .EQ. 4 ) THEN
            NUTYEL = 19
         ELSE IF( NBSOMM .EQ. 6 ) THEN
            NUTYEL = 21
         ELSE IF( NBSOMM .EQ. 8 ) THEN
            NUTYEL = 23
         ELSE
            GOTO 100
         ENDIF
         CALL ELTYCA( NUTYEL )
C
C        BOUCLE SUR LES FACES DE L ELEMENT
C        ----------------------------------
         DO 60 I=1,NFACE
C           LE NO GLOBAL DES POINTS SOMMETS DE LA FACE I
            NS = NBSOFA(I)
            DO 27  J=1,NS
               NGS(J) = NOSOCU( NOSOFA(J,I) , N )
   27       CONTINUE
C
C           PERMUTATION CIRCULAIRE DES SOMMETS POUR AMENER
C           LE PLUS PETIT NO EN PREMIER
            ISENS = 1
            L     = NGS(1)
            DO 29  J=2,NS
               IF( L .LE. NGS(J) ) GOTO 29
               L     = NGS(J)
               ISENS = J
   29       CONTINUE
            CALL TRTATA( NGS(ISENS) , NGS1(1) , NS-ISENS+1 )
            IF(ISENS .EQ. 1) GOTO 31
            CALL TRTATA( NGS(1) , NGS1(NS-ISENS+2) , ISENS-1 )
C
C           LA NUMEROTATION DES SOMMETS EST MISE A JOUR POUR
C           QUE LA 1-ERE ARETE DE LA FACE SOIT
C           1:NO LE PLUS FAIBLE DES SOMMETS DE LA FACE
C           2:NO MIN(SOMMET2 , SOMMET NS)
C           SI MIN =SOMMET2 FACE DIRECTE
C                   SOMMET NS FACE INDIRECTE
C           UNE FACE DIRECTE EST VUE DE L EXTERIEUR AU SOLIDE
C           SOUS LA FORME DIRECTE
C
   31       IF( NGS1(2) .LT. NGS1(NS) ) GOTO 33
C
C           FACE INDIRECTE
            ISENS    = - ISENS
            L        = NGS1(2)
            NGS1(2)  = NGS1(NS)
            NGS1(NS) = L
C
C           RECHERCHE OU ADJONCTION DE LA FACE
C           ----------------------------------
   33       CALL HACHAG(NS,NGS1,L1FACE,L2FACE,MCN(MNLFAC),5,
     &                  LIBREF,NOFAC)
            MN = MNLFAC + L1FACE * ( ABS(NOFAC) - 1 )
C
C           STOCKAGE ADRESSE DE L ELEMENT DANS NFACE(6,...;NOFAC)
C           -----------------------------------------------------
            L1 = 4
   35       L1 = L1 + 1
            IF( L1 .LT. L1FACE     ) THEN
               IF( MCN(MN + L1) .NE. 0) GOTO 35
C
C              ADRESSE +1 > 0 SI LA FACE EST DIRECTE DANS L ELEMENT
C              ADRESSE +1 < 0 SI LA FACE EST INDIRECTE DANS L ELEMENT
               MCN( MN + L1 ) = N * ISENS
C
            ELSE
C
C              ERREUR IL Y A PLUS DE MXMOFA CUBES CONTENANT CETTE FACE
               IF( NS.LT.4 ) NGS1(4)=0
               NBLGRC(NRERR) = 1
               KERR(1) = 'TROP DE CUBES CONTIENNENT CETTE FACE'
               CALL LEREUR
               WRITE(IMPRIM,10060)  N,(NOSOCU(L,N),L=1,8),
     %                             (NGS1(L),L=1,4),MXMOFA
            ENDIF
   60    CONTINUE
10060 FORMAT(' GEFAVO:DANS LE CUBE',I7,' DE SOMMETS',8I8/
     %' LA FACE DE SOMMETS',4I7,' APPARTIENT A PLUS DE',I3,' ELEMENTS'/
     %' VOIR MXMOFA DU SP GEFAVO'/)
C
C        PASSAGE AU CUBE SUIVANT
  100 CONTINUE
C
      END
