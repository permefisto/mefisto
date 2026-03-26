      SUBROUTINE TRISO2( NBISO,  VALISO, NTDISO,
     %                   NDIM  , NCAS,   NCAS0,  NCAS1, NTDL,
     %                   NTYP,   TEMPER, dptemp,
     %                   NUTYEL, NBELEM,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL,
     %                   NBCOOR, NBPOIT, XYZPOI,
     %                   MXSOMM, MXPOIS, SOLEL,  COPOE,
     %                   MXPILE, LAPILE, NOPOTR ,
     %                   COPOT , VALST , FBASE  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN 2D LES ISOTHERMES DES ELEMENTS FINIS D'UN TYPE
C -----    DONNE SOIT TRIANGLE SOIT QUADRANGLE

C ENTREES :
C ---------
C NBISO  : NOMBRE D'ISOVALEURS A TRACER
C VALISO : TABLEAU DES VALEURS DES ISO A TRACER
C NTDISO : TABLEAU DU NOMBRE DE TRACES DE SEGMENTS PAR ISOVALEURS
C          (0 EN ENTREE)
C NDIM   : ESPACE DE TRAVAIL 2
C NCAS   : NUMERO DU CAS A TRAITER
C NCAS0  : NUMERO DU PREMIER VECTEUR TEMPERATURE
C NCAS1  : NUMERO DU PREMIER VECTEUR TEMPERATURE
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NTYP   : =0 EMPLOI de TEMPER
C          =1 EMPLOI de dptemp
C TEMPER : TABLEAU DES NTDL * NCAS0:NCAS1 TEMPERATURES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NTDL) TEMPERATURES

C NUTYEL : NUMERO DU TYPE D'EF A TRAITER ICI
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT  FINI  DE CE TYPE
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS FINIS DE CE TYPE
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT  FINI  DE CE TYPE
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS FINIS DE CE TYPE

C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT (3 ou 6)
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE

C TABLEAUX AUXILIAIRES :
C ----------------------
C MXSOMM : NOMBRE MAXIMUM DE SOMMETS DES SOUS-TRIANGLES GENERES
C MXPOIS : NOMBRE MAXIMUM DE POINTS DES ISOVALEURS SUR LE TRIANGLE
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C MXPILE : NBRE MAXIMUM DE SOUS-TRIANGLES DANS LA PILE LAPILE
C LAPILE : PILE DES 3 SOMMETS DES SOUS-TRIANGLES
C NOPOTR : NUMERO DES POINTS SUR LES ISO D'UN TRIANGLE
C COPOT  : COORDONNEES DES POINTS ISO SUR LE  TRIANGLE
C VALST  : VALST(1 ET 2,.) 2 COORDONNEES DES SOMMETS DES SOUS TRIANGLES
C          VALST(3,.) VALEUR DE LA SOLUTION EN CE SOMMET
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1994
C23456---------------------------------------------------------------012
      PARAMETER   ( LIGCON=0, LIGTIR=1 )
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NUNDEL(NBELEM,NBNOEL),
     %                  NUPTEL(NBELEM,NBPOEL),
     %                  LAPILE(3,MXPILE),
     %                  NOPOTR(3,NBISO),
     %                  NTDISO(1:NBISO)
      REAL              XYZPOI(NBCOOR,NBPOIT)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp
      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1),
     %                  SOLEL(*),
     %                  COPOE(NBPOEL,NDIM)
      DOUBLE PRECISION  XD,YD,ZD,PROSCD,
     %                  FBASE(30),
     %                  XREF(3)
      REAL              VALISO(1:NBISO),
     %                  VALST(3,MXSOMM),
     %                  COPOT(2,MXPOIS)
      REAL              COSOTR(2,3), C(3), C1(3)
      INTEGER           NOS(6)
      CHARACTER*8       KNOM
      EQUIVALENCE      (XD,XREF(1)),(YD,XREF(2)),(ZD,XREF(3))
      DATA              ZD /0.D0/
C
C     CALCUL DU NOMBRE DE SUBDIVISIONS DE (0,1)
C     EN FONCTION DU NOMBRE DE NOEUDS DE L ELEMENT FINI
      IF( NBNOE .GT. 4 ) THEN
         NBSOUI = 2
      ELSE
         NBSOUI = 1
      ENDIF
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
      CALL ELINTE( 'THERMIQUE',NUTYEL,NDIMF,NOINTF,
     &              NBINVA,NUINVA,NUINTI,NBNDIN )
C     LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
      NOINTE = NUINTI(1)
C
CCC      IDEM POUR LE GRADIENT DE TEMPERATURE
CCC      CALL ELINT1( 'THER',NUTYEL,NDIMF,NOINFD,NBPOF,
CCC     &              NBINVA,NUINVA,NUINTI,NBNDIN,NUNOIN )
C
C     LA DIMENSION DE TRACE EST NDIM
      IF( NDIM .NE. NDIMF ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'INCOMPATIBILITE DIMENSION ESPACES POUR F ET X'
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NARET .LT. 3  .OR.  NARET .GT. 4 ) THEN
C        ELEMENT NON TRIANGULAIRE OU QUADRANGULAIRE
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR : ELEMENT FINI 2D'
         KERR(2) = 'NON TRIANGLE NON QUADRANGLE'
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TRACE DES ARETES
C     -------------------
      CALL XVTYPETRAIT( NTLAFR )
      CALL XVEPAISSEUR( 1 )
      CALL TRAREF( NDIM, NBPOIT, XYZPOI, NBELEM, NUPTEL )
C
C     LE TRACE DES ISOTHERMES
C     -----------------------
      CALL XVTYPETRAIT( LIGCON )
      CALL XVEPAISSEUR( 2 )
      NBCOUL = NDCOUL - N1COUL + 1
C
C     OUVERTURE DE LA PILE DU NO DES 3 SOMMETS DES SOUS-TRIANGLES
C     -----------------------------------------------------------
C     LHPIL0 POINTE SUR LE SOMMET DE LA PILE LAPILE
      LHPIL0 = 0
C     NOSOM0 NO DU DERNIER SOMMET GENERE DANS LES SOUS-TRIANGLES
      NOSOM0 = 0
C
C     GENERATION DES SOUS-TRIANGLES DU TRIANGLE (0,0),(1,0),(0,1)
C     =============================================================
C     LES COORDONNEES DU TRIANGLE DE REFERENCE
C     ----------------------------------------
      I = 1
      COSOTR(1,1) = 0.
      COSOTR(2,1) = 0.
      COSOTR(1,2) = 1.
      COSOTR(2,2) = 0.
      COSOTR(1,3) = 0.
      COSOTR(2,3) = 1.
C
C     LE NO DES 3 SOMMETS DES NBSOUI**2 SOUS-TRIANGLES DANS
C     LE TRIANGLE REFERENCE
C     -----------------------------------------------------
 20   CALL EMPITR( NBSOUI,NOSOM0,MXSOMM,MXPILE,LHPIL0,LAPILE )

C     LES COORDONNEES DES (NBSOUI+1) ** 2 SOMMETS DES SOUS-TRIANGLES
C     --------------------------------------------------------------
      CALL COPOTR( COSOTR, NBSOUI+1, 3, VALST(1,NOSOM0+1) )
      NOSOM0 = NOSOM0 + (NBSOUI+2) * (NBSOUI+1) / 2

C     SI L ELEMENT EST UN QUADRANGLE GENERATION DES SOUS-TRIANGLES
C     DU TRIANGLE COMPLEMENTAIRE DE SOMMETS (1,1), (0,1), (1,0)
C     ============================================================
      IF( NARET .EQ. 4  .AND.  I .EQ. 1 ) THEN
         I = 2
         COSOTR(1,1) = 1.
         COSOTR(2,1) = 1.
         COSOTR(1,2) = 0.
         COSOTR(2,2) = 1.
         COSOTR(1,3) = 1.
         COSOTR(2,3) = 0.
         GOTO 20
      ENDIF

C     =====================================
C     BOUCLE SUR LES EF COURANTS DE CE TYPE
C     =====================================
      DO 1000 K=1,NBELEM

C        FORMATION DES TABLEAUX COPOE(NBPOEL,NDIM) ET SOLEL(NBNOEL)
C        COORDONNEES DES POINTS ET VALEUR DE LA SOLUTION NOSOL AUX NOEUDS
C        ================================================================
C        EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER AUX NOEUDS DE L'EF
         DO L=1,NBNOEL
            NDL = NUNDEL(K,L)
            IF( NTYP .EQ. 0 ) THEN
               SOLEL( L ) = TEMPER( NDL, NCAS )
            ELSE
               SOLEL( L ) = dptemp( NCAS )%dptab( NDL )
            ENDIF
         ENDDO

C        EXTRACTION DES COORDONNEES DES POINTS DE L'EF
         DO 28 I=1,NBPOEL
            DO 27 J=1,NDIM
               COPOE(I,J) = XYZPOI( J, NUPTEL(K,I) )
 27         ENDDO
 28      ENDDO
C
C        ICI LE TABLEAU SOLEL DOUBLE PRECISION CONTIENT LA VALEUR
C        DE LA SOLUTION AUX NOEUDS DE SON INTERPOLATION DANS L EF
C        DE REFERENCE ET LE TABLEAU COPOE(NBPOEL,NDIM) LES COORDONNEES
C        DES POINTS DANS L'EF DE REFERENCE OU L INTERPOLATION EST EFFECTUEE
C
C        BOUCLE SUR LES SOMMETS DES SOUS-TRIANGLES DANS L'EF COURANT
C        VALEUR AUX SOMMETS DES SOUS-TRIANGLES DANS L'EF
C        ===========================================================
         NOSOMM = NOSOM0
         DO 78 I=1,NOSOMM
            XD = VALST(1,I)
            YD = VALST(2,I)
C           LA VALEUR DES NBN FONCTIONS DE BASE EN (XD,YD,0D0)
            CALL INTERP( NOINTE, XD,YD,ZD, NBN, FBASE )
C           LA VALEUR DE LA SOLUTION EN CE POINT
            VALST(3,I) = REAL( PROSCD( FBASE, SOLEL, NBN ) )
   78    ENDDO
C
C        LA HAUTEUR INITIALE DE LA PILE
         LHPILE = LHPIL0
C
C        **************************************
C        TANT QUE LA PILE3 N EST PAS VIDE FAIRE
C        **************************************
C
   80    IF( LHPILE .GT. 0 ) THEN
C
C           DANS LE SOUS-TRIANGLE : LES 3 SOMMETS SONT DEPILES
C           --------------------------------------------------
            FMIN =  1.E30
            FMAX = -FMIN
            DO 90 I=1,3
C              LE NUMERO DES 3 SOMMETS DU SOUS-TRIANGLE
               NOSI   = LAPILE(I,LHPILE)
               NOS(I) = NOSI
               F1     = VALST(3,NOSI)
C              VALEUR MIN ET MAX DANS LE SOUS-TRIANGLE
               IF( F1 .GT. FMAX ) FMAX = F1
               IF( F1 .LT. FMIN ) FMIN = F1
   90       ENDDO
            LHPILE = LHPILE - 1
C
C           SI LES 3 VALEURS AUX SOMMETS DU TRIANGLE SONT EXTERIEURES
C           A L INTERVALLE DES ISOS LE TRIANGLE EST ABANDONNE
C           =========================================================
            IF( FMIN .GT. VALISO(NBISO) .OR. FMAX .LT. VALISO(1))GOTO 80
C
C           NOP NOMBRE DE POINTS DES ISOVALEURS SUR CE TRIANGLE
C           ---------------------------------------------------
            NOP = 0
C
C           MISE A ZERO DU TABLEAU NOPOTR
C           -----------------------------
            CALL AZEROI( 3*NBISO, NOPOTR )
C
C           KMIN,KMAX NO DE L ISOVALEUR MIN ET MAX SUR LE TRIANGLE
C           ------------------------------------------------------
            KMIN = NBISO + 1
            KMAX = 0
C
C           BOUCLE SUR LES 3 ARETES DU SOUS-TRIANGLE
C           ========================================
            DO 200 I=1,3
C
C              NO DES EXTREMITES DE L ARETE
               N1 = NOS(I)
               IF( I .NE. 3 ) THEN
                  N2 = I + 1
               ELSE
                  N2 = 1
               ENDIF
               N2 = NOS(N2)
C
C               LA VALEUR DE LA SOLUTION AUX EXTREMITES
C
                F1 = VALST(3,N1)
                F2 = VALST(3,N2)

ccc                IF( F1 - F2 ) 120, 100, 110

                IF( F1 .EQ. F2 ) GOTO 100
                IF( F1 .GT. F2 ) THEN
                   GOTO 110
                ELSE
                   GOTO 120
                ENDIF
C
C               F1 = F2 LA SOLUTION EST ELLE STATIONNAIRE SUR L ELEMENT?
C               --------------------------------------------------------
  100           N3 = I + 2
                IF(N3 .GT. 3) N3 = N3 - 3
                N3 = NOS(N3)
                F3 = VALST(3,N3)
                IF( F1 .EQ. F3 ) GOTO 80
C
C               SI F1=F2=F3 PAS DE TRACE D ISOVALEUR.PASSAGE AU TRIANGLE
C                           SUIVANT  (=> 80)
C
C               SI F1=F2=/F3 LES 2 EXTREMITES FORMENT EVENTUELLEMENT UN
C                            SEGMENT D ISOVALEUR
C
                DO 105 KK=1,NBISO
                   IF( F1 .EQ. VALISO(KK) ) THEN
C                     F1=F2 EST LA KK-EME ISOVALEUR
C                     ADJONCTION DU SEGMENT N1-N2
C                     -----------------------------
C                     LES EXTREMITES DE L ARETE SONT AJOUTEES
                      CALL POKISO( KK,NOP,VALST(1,N1),COPOT,NOPOTR,IERR)
                      CALL POKISO( KK,NOP,VALST(1,N2),COPOT,NOPOTR,IERR)
                      IF( IERR .NE. 0 ) GOTO 300
                      KMIN = MIN0( KMIN, KK )
                      K2   = KK
                      GOTO 180
                   ENDIF
  105           ENDDO
                GOTO 200
C
C               F1>F2 PERMUTATION DE (N1,N2) ET (F1,F2)
C               ---------------------------------------
  110           KK = N1
                N1 = N2
                N2 = KK
                F3 = F1
                F1 = F2
                F2 = F3
C
C               F1<F2 . RECHERCHE DE LA PLUS PETITE ISOVALEUR SUR
C                       LE SOUS-TRIANGLE
C               -----   -----------------------------------------
  120           DO 130 KK=1,NBISO
                   K1 = KK
                   IF( VALISO(KK) .GE. F1 ) GOTO 140
  130           ENDDO
                GOTO 200
C
C               LE TRIANGLE INTERSECTE AU MOINS UNE ISOVALEUR
C               ---------------------------------------------
  140           KMIN = MIN0(KMIN,K1)
                K2   = K1
C
C               BOUCLE SUR LES ISOVALEURS K1 A NBISO
C               ------------------------------------
                DO 170 KK=K1,NBISO
                   IF( VALISO(KK) .GT. F2 ) GOTO 180
C
C                  L ISOVALEUR KK EST COMPRISE ENTRE F1 ET F2
                   K2 = KK
C
C                  ADJONCTION DU POINT D INTERSECTION ISO KK ET ARETE I
C                  ---------------------------------------------------
                   Z    = (VALISO(KK) - F1) / (F2 - F1)
                   C(1) = VALST(1,N1) + Z  * (VALST(1,N2)-VALST(1,N1))
                   C(2) = VALST(2,N1) + Z  * (VALST(2,N2)-VALST(2,N1))
                   CALL POKISO( KK, NOP, C, COPOT, NOPOTR, IERR )
                   IF( IERR .NE. 0 ) GOTO 300
  170           ENDDO
C
C               NO DE LA PLUS GRANDE ISOVALEUR TRAITEE SUR LE TRIANGLE
  180           KMAX = MAX0(KMAX,K2)
C
C               FIN DE LA BOUCLE SUR LES ARETES DU TRIANGLE
  200        ENDDO
C
C            TRACE DES SEGMENTS DE CHAQUE ISOVALEUR DU SOUS-TRIANGLE
C            =======================================================
             IF ( KMIN .LE. KMAX ) THEN
                DO 210 KK=KMIN,KMAX
                   IF(NOPOTR(1,KK) .LE. 1) GOTO 210
C                  LE TRACE DE L'ISOTHERME AVEC LA BONNE COULEUR
                   IF( NOIBLA .EQ. 0 ) THEN
C                     LA COULEUR DE L'ISOTHERME
                      NCOUL = NBISO + 1 - KK
                      NCOUL = MOD(NCOUL-1,NBCOUL) + N1COUL
                   ENDIF
C
C                  IL EXISTE UN SEGMENT A TRACER
C                  CALCUL DES COORDONNEES DANS LE SOUS-TRIANGLE COURANT
C                  A PARTIR DE LA TRANSFORMATION F: EF REFERENCE -> EF
C                  ----------------------------------------------------
                   N1 = NOPOTR(2,KK)
                   XREF(1) = COPOT(1,N1)
                   XREF(2) = COPOT(2,N1)
                   CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C )
C
                   N2 = NOPOTR(3,KK)
                   XREF(1) = COPOT(1,N2)
                   XREF(2) = COPOT(2,N2)
                   CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C1 )
C
C                  L'ARETE ISO-VALEUR SUR L'EF COURANT VA DE C A C1
                   CALL TRAIT2D( NCOUL, C(1), C(2),  C1(1), C1(2) )
C
C                  TRACE DU NO DE L ISOVALEUR
C                  --------------------------
                   NTDISO( KK ) = NTDISO( KK ) + 1
                   IF( MOD( NTDISO(KK), 16 ) .EQ. 3 ) THEN
                      WRITE( KNOM, '(I8)' ) KK
                      DO 208 J=7,1,-1
                         IF( KNOM(J:J) .EQ. ' ' ) THEN
                            I1 = J + 1
                            GOTO 209
                         ENDIF
  208                 ENDDO
                      I1 = 1
  209                 CALL TEXTE2D( NCOUL,
     %                             (C(1)+C1(1))*0.5,
     %                             (C(2)+C1(2))*0.5, KNOM(I1:8) )
                   ENDIF
  210           ENDDO
             ENDIF
C
C            REMONTEE POUR TRAITER LE SOUS-TRIANGLE SUIVANT
C            ==============================================
             GOTO 80
C
C*******************************************************************************
C            IL Y A AU MOINS 3 POINTS D INTERSECTION D UNE ISOVALEUR
C            AVEC LE TRIANGLE . IL EST SUBDIVISE A SON TOUR EN 4 SOUS-TRIANGLES
C            ==================================================================
C            ACTUELLEMENT CE CAS NE PEUT PAS SE PRODUIRE EN P1
C            OU LA RECHERCHE DU ZERO DE L INTERPOLATION EST FAITE
C            IL PEUT SEULEMENT Y AVOIR 2 POINTS D INTERSECTION PAR TRIANGLE
C            DANS UN AUTRE CONTEXTE, CETTE RECHERCHE DES ZEROS RESTE A PROGRAMME
C
  300        WRITE(IMPRIM,10300) LHPILE+1,(J,NOPOTR(1,J),J=1,NBISO)
10300 FORMAT('ERREUR ELEMENT',I7,' 3 POINTS D INTERSECTION ISO-ARETE'/
     &       5(' NOPOTR(1,',I4,')=',I5))
             CALL XVPAUSE
C
C            GENERATION DU MILIEU DES 3 ARETES DU SOUS-TRIANGLE
C            ==================================================
             DO 320 I=1,3
                NOSOMM = NOSOMM + 1
C
C               CONTROLE DE SATURATION
C               ----------------------
                IF( NOSOMM .GT. MXSOMM ) THEN
C                  TABLEAU VALST TROP PETIT.DIAGNOSTIC.ARRET
                   NBLGRC(NRERR) = 2
                   KERR(1) = 'TAILLE VALST INSUFFISANTE'
                   KERR(2) = 'AUGMENTER MXSOMM DANS SP TRISOT'
                   CALL LEREUR
                   RETURN
                ENDIF
C
C               LE NO DU MILIEU DE L ARETE ET SES COORDONNEES
C               ---------------------------------------------
                NOS(I+3) = NOSOMM
                J = I + 1
                IF( J .GT. 3 ) J = 1
                DO 315 KK=1,2
                   VALST(KK,NOSOMM) = ( VALST(KK,NOS(I))
     &                              +   VALST(KK,NOS(J)) ) * 0.5
  315           ENDDO
C
C               LA VALEUR DE LA SOLUTION EN CE MILIEU D ARETE
C               ---------------------------------------------
                XD = VALST(1,NOSOMM)
                YD = VALST(2,NOSOMM)
C               LA VALEUR DES FONCTIONS DE BASE EN CE POINT
                CALL INTERP( NOINTE,XD,YD,ZD,NBN,FBASE )
C               LA VALEUR DE LA SOLUTION EN CE POINT
                VALST(3,NOSOMM) = REAL( PROSCD( FBASE, SOLEL, NBN ) )
  320        ENDDO
C
C            LES 4 SOUS-TRIANGLES SONT EMPILES DANS LAPILE
C            ============================================
             CALL EMPIL3( LHPILE,MXPILE,LAPILE,NOS(1),NOS(4),NOS(6) )
             CALL EMPIL3( LHPILE,MXPILE,LAPILE,NOS(4),NOS(2),NOS(5) )
             CALL EMPIL3( LHPILE,MXPILE,LAPILE,NOS(5),NOS(6),NOS(4) )
             CALL EMPIL3( LHPILE,MXPILE,LAPILE,NOS(6),NOS(5),NOS(3) )
             GOTO 80
C*******************************************************************************
C
C            ******************************************************
C            FIN DU TRAITEMENT DES SOUS-TRIANGLES DE L ELEMENT FINI
C            ******************************************************
         ENDIF
 1000 ENDDO
C
C     RETOUR A LA VALEUR PAR DEFAUT
      CALL XVEPAISSEUR( 1 )
C
      RETURN
      END
