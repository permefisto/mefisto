      SUBROUTINE GEARFA( XYZSOM, NBSOEF, NBEFOB, NOSOEF, MXFAAR,
     %                   L1ARFA, L2ARFA, MNARFA, NBARXF, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FORMER PAR HACHAGE LE TABLEAU DES ARETES DES QTANGLES
C -----    D'UNE SURFACE
C          UNE ARETE APPARTIENT A AU PLUS MXFAAR QTANGLES
C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE DE LA SURFACE
C NBSOEF : NOMBRE DE SOMMETS ET TANGENTES PAR QTANGLE
C          (4 SANS TG ET 12 AVEC TG )
C NBEFOB : NOMBRE DE FACES QTANGLES DE LA SURFACE
C NOSOEF : LES 4 NUMEROS DES SOMMETS DES NBEFOB FACES DE LA SURFACE
C MXFAAR : NOMBRE MAXIMUM DE FACES QTANGLES POUR UNE ARETE DANS NARFA
C          0 SI PAS DE DEMANDE D'UN TEL STOCKAGE
C MNARFA : ADRESSE DANS MCN DU TABLEAU NARFA DES QTANGLES DU MAILLAGE
C          >0 DESTRUCTION DE L'ANCIEN TABLEAU AVEC LES ANCIENNES VALEURS
C             L1ARFA L2ARFA QUI SONT MISES A JOUR ENSUITE
C          =0 PAS DE DESTRUCTION
C SORTIES:
C --------
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE DE QTANGLES DU TABLEAU NARFA
C MNARFA : >0 ADRESSE DANS MCN DU TABLEAU NARFA DES QTANGLES DU MAILLAGE
C          EN SORTIE NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C                    NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C                    NARFA(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C                    NARFA(4:3+MXFAAR,I)= NO NOSOEF DU QTANGLE
C                                         CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE MXFAAR QTANGLES, 
C          LE NUMERO DE QTANGLE MXFAAR EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES QTANGLES EST INCOMPLETE
C NBARXF : NBARXF(n) NOMBRE D'ARETES APPARTENANT A n  QTANGLES(ou QT) n=1,2
C          NBARXF(3) NOMBRE D'ARETES APPARTENANT A >2 QTANGLES
C IERR   : -1 000 000 SI L2ARFA INSUFFISANT DU TABLEAU ARFA DES ARETES
C           1 SI MXFAAR TROP PETIT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    Novembre 1988
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Novembre 2015
C ......................................................................
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              XYZSOM(3,*)
      INTEGER           NOSOEF(NBSOEF,NBEFOB), NBARXF(6)
      INTEGER           NS(2)

C     DESTRUCTION DE L'ANCIEN TABLEAU DES ARETES
C     ------------------------------------------
      IF( MNARFA.GT.0 ) CALL TNMCDS('ENTIER',L1ARFA*L2ARFA,MNARFA)

C     MAJORATION DU NOMBRE DES ARETES DE LA SURFACE
C     ---------------------------------------------
      L1ARFA = 3 + MXFAAR

C     LE NOMBRE MAXIMUM D'ARETES EN L'ABSENCE DU NOMBRE DE SOMMETS
cccC     FORMULE EMPIRIQUE
ccc      IF( NBEFOB .LT. 1024 ) THEN

C     MAJORATION BESTIALE MAIS SURE
      L2ARFA = 4 * NBEFOB

ccc      ELSE
ccc         L2ARFA = NINT( 3.5 * NBEFOB )
ccc      ENDIF

C     ADRESSAGE DU TABLEAU ARFA
C     -------------------------
      L      = L1ARFA * L2ARFA
      MNARFA = 0
      CALL TNMCDC( 'ENTIER', L, MNARFA )

C     LE TABLEAU DES ARETES EST INITIALISE A ZERO
      CALL AZEROI( L, MCN(MNARFA) )

C     LA 1-ERE ARFA LIBRE A PARTIR DES DERNIERES PAR VALEUR DECROISSANTE
      LIBREA = L2ARFA

C     FORMATION DU TABLEAU DES NO DES SOMMETS DES ARETES DES QTANGLES
C     ===============================================================
      NBARERR = 0
      DO 100 NF = 1, NBEFOB

C        NF QTANGLE ACTIF?
         IF( NOSOEF(1,NF) .EQ. 0 ) GOTO 100

C        LE NOMBRE DE SOMMETS DU QTANGLE NF ACTIF
         IF( NOSOEF(4,NF) .EQ. 0 ) THEN
            NBSTFA = 3
         ELSE
            NBSTFA = 4
         ENDIF

C        BOUCLE SUR LES ARETES DU QTANGLE
C        --------------------------------
         J1 = NBSTFA
         DO 30 J=1,NBSTFA

C           L'ARETE J DU QTANGLE NF
            IF( NOSOEF(J1,NF) .LE. NOSOEF(J,NF) ) THEN
               NS(1) = NOSOEF( J1, NF )
               NS(2) = NOSOEF( J , NF )
            ELSE
               NS(1) = NOSOEF( J , NF )
               NS(2) = NOSOEF( J1, NF )
            ENDIF

C           ADJONCTION DE L'ARETE SI ELLE N'EXISTE PAS DEJA
            CALL HACHAG( 2, NS, L1ARFA, L2ARFA, MCN(MNARFA), 3,
     %                   LIBREA, NAR )
            J1 = J

C           AJOUT DU NUMERO DU QTANGLE DE CETTE ARETE
C           NAR = 0 SI LE TABLEAU MNARFA EST SATURE
            IF( NAR .EQ. 0 ) THEN
C
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) =
     %            'gearfa: SATURATION du TABLEAU ARETES des QTANGLES'
               ELSE
                  KERR(1) ='gearfa: EDGES ARRAY of SURFACE IS SATURATED'
               ENDIF
               CALL LEREUR
               IERR = -1 000 000
               RETURN

            ENDIF

C           NAR > 0 SI ARETE RETROUVEE
C           NAR < 0 SI ARETE AJOUTEE
            IF( MXFAAR .LE. 0 ) GOTO 30

C           ADRESSE DE L'ARETE DANS LE TABLEAU ARFA
            MN = MNARFA + L1ARFA * ( ABS(NAR) - 1 )

C           LE K-EME QTANGLE DE L'ARETE NAR EST STOCKE
            DO K=1,MXFAAR
               M = MN + 2 + K
               IF( MCN( M ) .EQ. 0 ) THEN
C                 STOCKAGE DU K-EME QTANGLE DE L'ARETE
                  MCN( M ) = NF
                  GOTO 30
               ENDIF
            ENDDO

C           ICI PROBLEME: >MXFAAR QTANGLES POUR CETTE ARETE NAR
C           AFFICHAGE DE L'ARETE APPARTENANT A >MXFAAR QTANGLES
            NBARERR = NBARERR + 1
            WRITE(KERR(MXLGER)(1:10),'(I1)' ) MXFAAR
            NBLGRC(NRERR) = 1
            WRITE(IMPRIM,*)

            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'gearfa: ARETE COMMUNE A PLUS de ' //
     %                   KERR(MXLGER)(1:1) // ' QTANGLES'
               WRITE(IMPRIM,*)'gearfa: face',NF,' arete',J,
     %              ' St:',MCN(MN), MCN(MN+1),
     %              ' des QTANGLES',(MCN(MN+2+K),K=1,MXFAAR),NF
               DO K=0,1
                  WRITE(IMPRIM,10025) MCN(MN+K),
     %                 (XYZSOM(L,MCN(MN+K)),L=1,3)
               ENDDO
               DO K=1,MXFAAR
                  NFF = ABS( MCN(MN+2+K) )
                  WRITE(IMPRIM,*)'QTANGLE',NFF,
     %              ' Sommets:',(NOSOEF(L,NFF),L=1,4)
               ENDDO
               WRITE(IMPRIM,*)'QTANGLE',NF,
     %              ' Sommets:',(NOSOEF(L,NF),L=1,4)

            ELSE

               KERR(1) = 'gearfa: COMMON EDGE IN MORE of  ' //
     %                    KERR(MXLGER)(1:1) // ' QTANGLES'
               WRITE(IMPRIM,*)'gearfa: face',NF,' edge',J,
     %                        ' Vertices:',MCN(MN), MCN(MN+1),
     %         ' in QTANGLES',(MCN(MN+2+K),K=1,MXFAAR),NF
               DO K=0,1
                  WRITE(IMPRIM,20025) MCN(MN+K),
     %                 (XYZSOM(L,MCN(MN+K)),L=1,3)
               ENDDO
               DO K=1,MXFAAR
                  NFF = ABS( MCN(MN+2+K) )
                  WRITE(IMPRIM,*)'QTANGLE',NFF,
     %              ' Vertices:',(NOSOEF(L,NFF),L=1,4)
               ENDDO
               WRITE(IMPRIM,*)'QTANGLE',NF,
     %              ' Vertices:',(NOSOEF(L,NF),L=1,4)

            ENDIF

C           PROBLEME: >MXFAAR QTANGLES POUR CETTE ARETE NAR
            M = MN + 2 + MXFAAR
C           VALEUR IMPOSEE NEGATIVE DU NO DE LA DERNIERE FACE STOCKEE
C           EN ECRASANT LA DERNIERE VALEUR
            MCN(M) = -NF

            IF( NBARERR .LE. 1 ) THEN
               CALL LEREUR
            ENDIF

 30      CONTINUE
 100  CONTINUE

10025 FORMAT(' SOMMET',I7,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
20025 FORMAT(' VERTEX',I7,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)

C     CALCUL DU NOMBRE D'ARETES APPARTENANT A 1, 2, 3, 4, 5, 6 QTANGLES
C     -----------------------------------------------------------------
      CALL NBARXFA( L1ARFA, L2ARFA, MCN(MNARFA),
ccc     %              NBSOEF, NOSOEF,
     %              NBARXF, NOARPB, IERR )

      RETURN
      END
