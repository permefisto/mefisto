      SUBROUTINE CODTEI ( TEINTE , PRCTEI , DEGTEI , CODEE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :   CODE UNE TEINTE SELON NOS NORMES, C.A.D. DANS LE SYSTEME DE
C  ====    MUNSELL, ICI C'EST LA TEINTE EN DEGRES QUI EST DEFINIE.
C
C  PARAMETRE D'ENTREE :
C  ===================
C  TEINTE : LA TEINTE A CODER, CHAINE DE CARACTERES, T. Q.
C           * SOIT : ELLE CONTIENT DEUX TEINTES ( CAS 1/)
C           * SOIT : ELLE NE CONTIENT QU'UNE TEINTE (CAS 2/)
C
C  PRCTEI : POURCENTAGE ( DE 0 % , A 100. % ),
C           * CAS 1/ C'EST LE POURCENTAGE DE LA TEINTE A
C                    DEFINIR PRISE ENTRE LES DEUX TEINTES
C                    STIPULEES PAR ''TEINTE'' ,
C                    0.%   SI L'ON DESIRE LA 1 IERE TEINTE
C                    100.% SI L'ON DESIRE LA 2 IEME.
C              RMQ : LES DEUX TEINTES (DANS ''TEINTE'')
C              ----- PEUVENT ETRES EGALES, ET SI ''PRCTEI''
C                    VARIE ENTRE 0% ET 100%, ON OBTIENT
C                    TOUTES LES TEINTES DU CONE DE COULEUR.
C
C            * CAS 2/ N'A AUCUN SENS, (EST IGNOREE)
C
C  PARAMETRES DE SORTIE :
C  =====================
C  CODEE  : = VRAI, SI ''TEINTE'' APPARTIENT A LA TABLE DES TEINTES.
C           = FAUX, SINON.
C  SI LA TEINTE A ETE CODEE :
C  DEGTEI            ; LA VALEUR DE CETTE TEINTE EN DEGRES (DE 0 A 360)
C                                                 ( 0 DEGRES --> ROUGE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : SVIGA        ANALYSE NUMERIQUE UPMC  PARIS         AVRIL 1986
C2345X7..............................................................012
      CHARACTER*(*)   TEINTE
      LOGICAL         CODEE
      CHARACTER*80    TEIN , TEINTS(2)
      REAL            DEG(2)
      LOGICAL         SEULE
C
      SEULE = .FALSE.
C
C     ENLEVONS LES BLANCS DEVANT.
C     ===========================
      TEIN = TEINTE
      IPLA = 1
      DO 1000 I = 1 , 80
         IF (TEIN(I:I) .NE. ' ') THEN
            IPLA = I
            GOTO 1001
         ENDIF
1000  CONTINUE
C
1001  TEIN = TEIN(IPLA : 80)
C
C     LA PREMIERE TEINTE.
C     ===================
      DO 1002 I = 1 , 80
         IF (TEIN(I:I) .EQ.' '.OR. TEIN(I:I) .EQ.'_') THEN
            IPLA = I - 1
            GOTO 1003
         ENDIF
1002  CONTINUE
      IPLA = 80
C
1003  TEINTS(1) = TEIN(1 : IPLA)
C
C     LA SECONDE TEINTE.
C     ==================
C
      IP = IPLA + 1
      IF ( IP .GE. 80 ) THEN
         SEULE = .TRUE.
      ELSE
C
C        ENLEVER LES BLANCS.
C        -------------------
         DO 1004 I = IP , 80
            IF ( TEIN(I:I) .NE. ' ' .AND.
     &           TEIN(I:I) .NE. '_' )    THEN
               IPLA = I
               GOTO 1005
            END IF
1004     CONTINUE
         SEULE = .TRUE.
C
C        LA SECONDE TEINTE.
C        ------------------
 1005    IF ( .NOT. SEULE ) THEN
                   DO 1006 I = IPLA , 80
                      IF (TEIN(I:I) .EQ.' '.OR.
     &                    TEIN(I:I) .EQ.'_') THEN
                                                   IPLA2 = I - 1
                                                   GOTO 1007
                      END IF
1006               CONTINUE
                   IPLA2 = 80
1007               CONTINUE
                   TEINTS(2) = TEIN(IPLA:IPLA2)
                END IF
         END IF
C
C=======================================================================
C     DETERMINATION DES ANGLES DES DEUX TEINTES.
C     ==========================================
      CODEE = .TRUE.
      I = 1
 11   CONTINUE
C     ------------------------------------------------------------------
      IF      ( TEINTS(I)(1:5) .EQ. 'ROUGE' )      THEN
                                                         DEG(I) = 0.
      ELSE IF ( TEINTS(I)(1:5) .EQ. 'JAUNE' )      THEN
                                                         DEG(I) = 60.
      ELSE IF ( TEINTS(I)(1:4) .EQ. 'VERT' )       THEN
                                                         DEG(I) = 120.
      ELSE IF ( TEINTS(I)(1:4) .EQ. 'CYAN' )       THEN
                                                         DEG(I) = 180.
      ELSE IF ( TEINTS(I)(1:4) .EQ. 'BLEU' )       THEN
                                                         DEG(I) = 240.
      ELSE IF ( TEINTS(I)(1:7) .EQ. 'MAGENTA' )    THEN
                                                         DEG(I) = 300.
      ELSE
             CODEE = .FALSE.
      END IF
C
C     LA SECONDE TEINTE.
C     ------------------
      IF ( I .LE. 1 .AND. CODEE .AND. .NOT. SEULE ) THEN
               I = I + 1
               GOTO 11
      END IF
C
C     DETERMINATION DE DEGTEI.
C     ========================
      IF ( CODEE ) THEN
         IF ( SEULE ) THEN
             DEGTEI = DEG(1)
         ELSE
C
C           NOIR OU BLANC DANS ''TEINTE'' ?
C           -------------------------------
C           TEINTE INDEFINIE.
C           -----------------
                 DEGTEI = -1.
C
C                LA 1 ERE TEINTE EST AVANT LA SECONDE
C                -------------------------------- (DANS LE SENS TRIGO) ?
                 IF (DEG(2) .LT. DEG(1)) DEG(2) = DEG(2) + 360.
C
C                DEUX FOIS LA MEME TEINTE ?
C                --------------------------
                 IF (TEINTS(1) .EQ. TEINTS(2)) DEG(2) = DEG(1) + 360.
C
C                LA TEINTE EN DEGRES.
C                --------------------
                 D      = DEG(2) - DEG(1)
                 DEGTEI = DEG(1) + 0.01 * PRCTEI * D
                 IF ( DEGTEI .GE. 360. ) DEGTEI = DEGTEI - 360.
         END IF
      END IF
      END
