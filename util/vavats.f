      SUBROUTINE VAVATS( NOVATS , NCODEV , DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCUL DE LA VARIABLE TMS NOVATS DE NOM KVATS( NOVATS )
C -----
C
C ENTREES :
C ---------
C NOVATS : NUMERO DANS KVATS DU NOM DE LA VARIABLE TMS DE VALEUR
C          A CALCULER
C
C SORTIES :
C ---------
C NCODEV : 0 DBLVAL N'EST PAS INITIALISEE
C          1 DBLVAL EST INITIALISEE
C DBLVAL : VALEUR REELLE DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/msvaau.inc"
      COMMON /UNITES/LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NOAFTS,NUNIT(26)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      LOGICAL           LMCN(MOTMCN)
      DOUBLE PRECISION  DMCN(MOTMCN/2)
      EQUIVALENCE      (RMCN(1),MCN(1),LMCN(1),DMCN(1))
      DOUBLE PRECISION  DBLVAL
C
      IF( NOVATS .GT. 0 .AND. NOVATS .LE. NBVATS ) THEN
C
C        LE NOM DE LA VARIABLE TMS  NOVATS
         CALL VATSTD( KVATS(NOVATS) , NOTYPE , MNTMS , LDTMS , NOTYP )
         IF( MNTMS .GT. 0 ) THEN
C
C           VERIFICATION DE NOTYPE
            NOTYP = NOTYPE
            IF( NOTYP .EQ. 21 ) THEN
C              NUMERO DU TMS . VALEUR ENTIERE
               NOTYP = 4
            ELSE IF( NOTYP .LE. 0 .OR. NOTYP .GT. NBTYPV ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) =  'LU: SP VAVATS: TYPE INCORRECT'
               CALL LEREUR
               RETURN
            ENDIF
C
C           VERIFICATION DE L'ADRESSE
            IF( MNTMS .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(2)(1:12),'(I12)') MNTMS
               KERR(1) = 'LU: SP VAVATS:ADRESSE INCORRECTE MNTMS '//
     %                    KERR(2)(1:12)
               CALL LEREUR
               RETURN
            ENDIF
C
C           L'ADRESSE DANS MCN DU DEBUT DE LA VARIABLE
            MN = MNTMS + LDTMS
C
C           AFFICHAGE SELON LE TYPE . ICI MOT = ENTIER !
            GOTO(   10 , 9000 , 9000 ,  40 , 50 , 60 , 9000 , 9000 ,
     %            9000 ,   40 ,   40 , 9000, 9000 ) , NOTYP
C
C           LOGICAL
 10         IF( LMCN(MN) ) THEN
                DBLVAL = 1
            ELSE
                DBLVAL = 0
            ENDIF
            GOTO 500
C
C           ENTIER
 40         DBLVAL = MCN(MN)
            GOTO 500
C
C           REEL
 50         DBLVAL = RMCN(MN)
            GOTO 500
C
C           REEL2
 60         DBLVAL = DMCN( MN/2 + 1 )
C
C           AFFICHAGE DE LA LIGNE
 500        NCODEV = 1
            RETURN
         ENDIF
      ENDIF
C
C     SI PAR EXEMPLE "AFFICHER ligne>axe>xyzsommet(11)" EST ACTIF,
C     ALORS IL N'Y A PAS D'ERREUR!
C     EN EFFET, UN MENU PEUT ETRE ACTIF ET L'ENTREE D'UN ENTIER EST EN ATTENTE
C     AU LIEU DE CELA ARRIVE LA DONNEE DE AFFICHER UNE VARIABLE DE TYPE
C     DIFFERENT DE 1 4 5 6 10 11, PAR EXEMPLE 12 C A D LA VALEUR DE X ET Y ET Z
C     => ERREUR MAIS COMME L'AFFICHAGE A DEJA ETE FAIT IL SUFFIT DE
C     CONTINUER EN SUPPRIMANT L'ERREUR POUR NE PAS EFFRAYER L'UTILISATEUR
CCC 9000 IF( NOAFTS .NE. 0 ) GOTO 9900
CCCC
CCCC     ERREUR EFFECTIVE
CCC      NBLGRC(NRERR) = 1
CCC      KERR(1) ='LU: TYPE INCORRECT POUR LA VARIABLE TMS '//
CCC     %          KVATS(NOVATS)
CCC      CALL LEREUR
CCCC
CCC 9900 NCODEV = 0
C
C     DBLVAL VALEUR NON INITIALISEE
 9000 NCODEV = 0
      END
