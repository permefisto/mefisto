       CHARACTER*80 FUNCTION KINFO(MOTCLE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER UNE INFORMATION DEPENDANT MACHINE DE TYPE CHARACTER
C -----
C ENTREE :
C --------
C MOTCLE : NATURE DE L'INFORMATION A RETOURNER, EN CLAIR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : FREDERIC HECHT  PATRICK LAUG  ALAIN PERRONNET   NOVEMBRE 1994
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++072
      include"./incl/langue.inc"
      CHARACTER*(*)   MOTCLE
      CHARACTER*80    KTITRE, BUFFER
      CHARACTER*10    KKK
      COMMON /TRAVA1/ MTITRE(20), NDATE(2), NOMCRE(6)
      EQUIVALENCE    (MTITRE,KTITRE)
C
C     1. PARTIE PORTABLE
C     ------------------
      IF (MOTCLE .EQ. 'TITRE') THEN
         KINFO = KTITRE
C
C     2. PARTIE NON PORTABLE
C     ----------------------
      ELSE IF (MOTCLE .EQ. 'MACHINE') THEN
C       NOM DE LA MACHINE HOTE
        CALL NOMORDINATEURHOTE( BUFFER, NBCAR )
        KINFO = BUFFER(1:NBCAR) // '    '
C
      ELSE IF (MOTCLE .EQ. 'DATE') THEN
         CALL LADATE( IAA, IMM, IJJ )
C        IAA AN, IMM MOIS, IJJ JOUR
         WRITE(KKK,'(3I2)') mod(IAA,100), IMM, IJJ
         IF( LANGAG .EQ. 0 ) THEN
C           INVERSION JOUR AN
C           DATE DU JOUR (JJ/MM/AA) JJ JOUR, MM MOIS, AA ANNEE
            KINFO(1:2) = KKK(5:6)
            KINFO(3:3) = '/'
            KINFO(4:5) = KKK(3:4)
            KINFO(6:6) = '/'
            KINFO(7:8) = KKK(1:2)
         ELSE
C           DATE DU JOUR (AA/MM/JJ) AA ANNEE, MM MOIS, JJ JOUR
            KINFO(1:2) = KKK(1:2)
            KINFO(3:3) = '/'
            KINFO(4:5) = KKK(3:4)
            KINFO(6:6) = '/'
            KINFO(7:8) = KKK(5:6)
         ENDIF
C
      ELSE IF (MOTCLE .EQ. 'HEURE') THEN
C        HEURE (HHMMSSFFFFFF)
C        HH HEURE DE 00 A 23, MM MINUTE DE 00 A 59, SS SECONDE DE 00 A 59,
C        FFFFFF MICROSECONDE, DE 000000 A 999999
         CALL HEUREMINUTESECONDE( IHH, IMM, ISS, IMS )
         WRITE(KINFO,'(3(I2.2),I6.6)') IHH, IMM, ISS, IMS
C
      ELSE IF (MOTCLE .EQ. 'UTILISATEUR') THEN
C        NOM DE L'UTILISATEUR RECUPERE EN VARIABLE D'ENVIRONNEMENT
         CALL VALVARENV( 'LOGNAME'//CHAR(0), LEN(BUFFER),
     %                   BUFFER, LMS )
         IF( LMS .LT. 0 ) THEN
             WRITE (IINFO('I'),*) 'KINFO : VARIABLE LOGNAME NON DEFINIE'
             KINFO = ' '
         ELSE IF( LMS .EQ. 0 ) THEN
             WRITE (IINFO('I'),*) 'KINFO : VARIABLE LOGNAME VIDE'
             KINFO = ' '
         ELSE
            KINFO = BUFFER(1:LMS)
         ENDIF
         IF( INDEX(KINFO,'perronne ') .GT. 0  .OR.
     %       INDEX(KINFO,'ap ')       .GT. 0 ) THEN
            BUFFER='Alain Perronnet'
            KINFO = BUFFER
         ENDIF
C
      ELSE
C        ERREUR
         WRITE (IINFO('I'),*) 'KINFO : MOT-CLE INCONNU : ', MOTCLE
         KINFO = ' '
      END IF
      END
