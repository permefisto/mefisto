      SUBROUTINE GEARFP( NBFPPL , NOFPPL , NOFAPE , NOFACE , MOARFP ,
     &                   L1ARFP , L2ARFP , MNARFP , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FORMER PAR HACHAGE LE TABLEAU DES ARETES DES FACES PERDUES
C -----    DE L'OBJET TRAITE PAR VOEX09
C          ( 1 ARETE DANS AU PLUS MOARFP FACES )
C ENTREES:
C --------
C NBFPPL : NOMBRE DE FACES DE LA SURFACE
C NOFPPL : NUMEROS DANS NOFAPE DES FACES COPLANAIRES PERDUES
C NOFAPE : NUMEROS DES FACES PERDUES DANS NOFACE
C NOFACE : LES 4 NUMEROS DES SOMMETS DES FACES DE L'OBJET
C MOARFP : NOMBRE MAXIMAL DE FACES ADJACENTES PAR UNE MEME ARETE
C
C SORTIES:
C --------
C L1ARFP : NOMBRE DE MOTS PAR ARFP DU TABLEAU NARFP
C L2ARFP : NOMBRE DE FACES DU TABLEAU NARFP
C MNARFP : ADRESSE DANS M DU TABLEAU NARFP DES FACES DU MAILLAGE
C          EN SORTIE NARFP(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C                    NARFP(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C                    NARFP(3,I)= CHAINAGE HACHAGE SUR LA FACE SUIVANTE
C                    NARFP(4:4+MOARFP,I)= NUMERO DANS NOFPPL DES 2 FACES
C                                         ADJACENTES PAR CETTE ARETE
C IERR   : =0 SI PAS DE PROBLEME
C          >0 NOMBRE D'ARETES APPARTENANT A PLUS DE MOARFP FACES
C          <0 SATURATION DU TABLEAU ARFP
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   JANVIER 1989
C ....................................................................12
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           NOFACE(7,*),NOFPPL(NBFPPL),NOFAPE(*)
      INTEGER           NS(2)
C
      IERR = 0
C
C     CALCUL DU NOMBRE MAJORANT DES ARETES PERDUES
C     --------------------------------------------
      L1ARFP = 3 + MOARFP
      L2ARFP = 3 * NBFPPL
C
C     ADRESSAGE DU TABLEAU ARFP
C     -------------------------
      L      = L1ARFP * L2ARFP
      MNARFP = 0
      CALL TNMCDC( 'ENTIER' , L , MNARFP )
C
C     LE TABLEAU DES FACES EST INITIALISE A ZERO
      CALL AZEROI( L , MCN(MNARFP) )
C
C     LA 1-ERE ARFP LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
      LIBREF = L2ARFP
C
C     FORMATION DU TABLEAU DES NO DES SOMMETS DES ARETES DES FACES
C     ============================================================
      DO 100 N=1,NBFPPL
C        LE NUMERO DE LA FACE PERDUE DANS NOFAPE
         NF = NOFPPL( N )
         IF( NF .LE. 0 ) GOTO 100
C        LE NUMERO DE LA FACE PERDUE DANS NOFACE
         NF = ABS( NOFAPE( NF ) )
         IF( NF .LE. 0 ) GOTO 100
C        LA FACE EXISTE
         DO 50 I=1,3
C           LES NUMEROS DES 2 SOMMETS DE L'ARETE NS1 < NS2
            IF( I .NE. 3 ) THEN
               NS(1) = NOFACE(I  ,NF)
               NS(2) = NOFACE(I+1,NF)
            ELSE
               NS(1) = NOFACE(1,NF)
               NS(2) = NOFACE(3,NF)
            ENDIF
C
C           ADJONCTION DE L'ARETE SI ELLE N'EXISTE PAS DEJA
            CALL HACHAG( 2,NS,L1ARFP,L2ARFP,MCN(MNARFP),3,
     &                   LIBREF, NAR )
C
C           AJOUT DU NUMERO DE LA FACE DE CETTE ARETE
C           NAR > 0 SI ARETE RETROUVEE
C           NAR < 0 SI ARETE AJOUTEE
            MN = MNARFP + L1ARFP * ( ABS(NAR) - 1 )
            IF( NAR .LT. 0 ) THEN
C              ARETE VUE POUR LA 1-ERE FOIS
               MCN( MN + 3 ) = N
            ELSE IF( NAR .GT. 0 ) THEN
C              LA 2-EME OU 3-EME OU ... FACE AVEC CETTE ARETE
               DO 20 J=4,L1ARFP-1
                  IF( MCN( MN + J ) .EQ. 0 ) THEN
C                    LA J-2 EME FACE EST AJOUTEE
                     MCN( MN + J ) = N
                  ELSE
C                    MOARFP FACES PAR ARETE
                     IERR = IERR + 1
                     NBLGRC(NRERR) = 1
                     KERR(1) = KERR(MXLGER)(1:10)//
     %              'ARETES COMMUNES A PLUS DE 2 FACES'
                     CALL LEREUR
                     WRITE(IMPRIM,*)' GEARFP: ARETE',NAR,
     %               MCN(MN),MCN(MN+1), ' DANS LES FACES',
     %              (ABS(NOFAPE(NOFPPL(MCN(MN+K)))),K=3,L1ARFP-1),
     %              ' PLUS LA FACE ',NF, ' NON STOCKABLE'
                  ENDIF
 20            CONTINUE
            ELSE
               NBLGRC(NRERR) = 1
               KERR(1) ='GEARFP:SATURATION DU TABLEAU ARFP'
               CALL LEREUR
               IERR = -1 000 000
            ENDIF
 50      CONTINUE
 100  CONTINUE
      END
