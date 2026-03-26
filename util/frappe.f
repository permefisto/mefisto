      SUBROUTINE FRAPPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : OUVRIR LE FICHIER DE SAUVEGARDE DE LA FRAPPE
C ----- SAUVEGARDER LE NOMBRE D'EXECUTIONS EFFECTUEES
C
C       A EXECUTER APRES UN CALL MSOU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1990
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/nmproj.inc"
      COMMON /UNITES/   LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      CHARACTER*8       DATE
      CHARACTER*8       HEURE
      CHARACTER*28      KNOM
      CHARACTER*80      KINFO, KAUX, NOM
      CHARACTER*72      KLIG
C
C     OUVERTURE DU FICHIER DE SAUVEGARDE DE LA FRAPPE
      NFFRAP = 9
C
C     UNE EXECUTION DE PLUS
C     NUEXEC A ETE LU SUR LE SUPER FICHIER DANS MSOU
      NUEXEC = NUEXEC + 1
C
      KNOM   = 'frappe.'
      WRITE( KNOM(8:12) , '(I5)' ) NUEXEC
      DO 20 I=12,8,-1
         IF( KNOM(I:I) .EQ. ' ' ) GOTO 30
 20   CONTINUE
C
C     INSERTION DE 0 POUR AVOIR LES FICHIERS FRAPPE DANS L'ORDRE AVEC LS
C     JUSQU'A CONCURRENCE DE 1000 FICHIERS frappe.1   ...  frappe.999
 30   IF( I .EQ. 11 ) THEN
         KNOM = KNOM(1:7) // '00' // KNOM(12:12)
      ELSE IF( I .EQ. 10 ) THEN
         KNOM = KNOM(1:7) // '0' // KNOM(11:12)
      ELSE
         KNOM = KNOM(1:7) // KNOM(10:12)
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,10000) KNOM
      ELSE
         WRITE (IMPRIM,20000) KNOM
      ENDIF
10000 FORMAT('Les donnees frappees sont ecrites dans le fichier ',A)
20000 FORMAT('Typed data are written on file ',A)
C
      OPEN(FILE=KNOM,UNIT=NFFRAP,ACCESS='SEQUENTIAL',
     %     FORM='FORMATTED',STATUS='UNKNOWN')
C
C     ECRITURE DU NOMBRE D'EXECUTIONS SUR LA MS
      NPA    = 2
      NOFISF = 10
      WRITE( UNIT=NOFISF , REC=NPA ) NUEXEC,NMPROJ
C
C     UNE LIGNE DE =
      WRITE( NFFRAP , 10001 )
10001 FORMAT( '{ ',68('='),' }' )
C
C     ECRITURE DU NOM DE LA VERSION DE MEFISTO
      CALL VRSION( KNOM )
      IF( LANGAG .EQ. 0 ) THEN
         KLIG = '{ Logiciel MEFISTO     : ' // KNOM
      ELSE
         KLIG = '{ MEFISTO Software     : ' // KNOM
      ENDIF
      KLIG(71:72) = ' }'
      WRITE( NFFRAP, 10002 ) KLIG
10002 FORMAT(A72)
C
C     ECRITURE DU NOM DE L'UTILISATEUR
      NOM = KINFO( 'UTILISATEUR' )
      L   = NUDCNB( NOM )
      IF( LANGAG .EQ. 0 ) THEN
         KLIG = '{ Nom de l''Utilisateur : ' // NOM(1:L)
      ELSE
         KLIG = '{ USER''s Name          : ' // NOM(1:L)
      ENDIF
      KLIG(71:72) = ' }'
      WRITE( NFFRAP, 10002 ) KLIG

C     ECRITURE DU NOM DU PROJET
      L = NUDCNB( NMPROJ )
      IF( LANGAG .EQ. 0 ) THEN
         KLIG = '{ Nom du Projet        : ' // NMPROJ(1:L)
      ELSE
         KLIG = '{ Project Name         : ' // NMPROJ(1:L)
      ENDIF
      KLIG(71:72) = ' }'
      WRITE( NFFRAP, 10002 ) KLIG

C     ECRITURE DE LA DATE HEURE
      KAUX      = KINFO( 'DATE' )
      DATE(1:8) = KAUX(1:8)
      IF( DATE(4:4) .EQ. ' ' ) DATE(4:4) = '0'
      KLIG = '{ Date                 : ' // DATE(1:8)
      L    = NUDCNB( KLIG )
      KAUX       = KINFO( 'HEURE' )
      HEURE(1:8) = KAUX(1:8)
      KLIG       = KLIG(1:L) // '   ' //
     % HEURE(1:2) // 'h ' // HEURE(3:4) // 'm ' // HEURE(5:6) // 's'
      KLIG(71:72) = ' }'
      WRITE( NFFRAP, 10002 ) KLIG

C     UNE LIGNE DE =
      WRITE( NFFRAP, 10001 )

      RETURN
      END
