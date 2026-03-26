       SUBROUTINE PFRONT( MNNTAB, NBEF, MNNSEF, MNXYZ, NUMPT, NSFR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :	   CHERCHE SI LE SOMMET NUMPT EST A LA FRONTIERE DE LA
C -----    TRIANGULATION DE LA SURFACE
C          NUMPT EST AU BORD SI TOUS LES AUTRES SOMMETS DES TRIANGLES
C          ADJACENTS SONT DANS 2 TRIANGLES
C          (ALGORITHME INCORRECT S'IL EXISTE UN QUADRANGLE)
C ENTREES:
C --------
C MNNTAB : TABLEAU NO DES EF A EXAMINER
C NBEF   : NOMBRE DES EF DE NO DANS MNNTAB
C MNNSEF : ADRESSE MCN DU TABLEAU NSEF DU MAILLAGE
C MNXYZ  : ADRESSE MCN DU TABLEAU XYZSOMMET	
C NUMPT  : NUMERO DU SOMMET A SUPPRIMER
C
C SORTIES:
C --------
C NSFR   : 0 SI PAS A LA FRONTIERE,
C          1 SI SOMMET A LA FRONTIERE DES TRIANGLES
C         -1 EN PRESENCE DE QUADRANGLE (ALGORITHME FAUX)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : KCEHA & NICO    LJLL UPMC                                2000
C MODIFS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
C
C     ON DECLARE UN TABLEAU TEMPORAIRE
      NBSOM=MCN(MNXYZ+WNBSOM)
      CALL TNMCDC( 'ENTIER', 1+NBSOM, MNTEMP )
      DO 10 I=0,NBSOM
         MCN(MNTEMP+I)=0
 10   CONTINUE
C
      MNN = MNNSEF + WUSOEF - 5
      DO 20 I=1,NBEF
C
C        NUMERO DE L'EF A EXAMINER
         NUMEF=MCN(MNNTAB+I-1)
C
C        NUMERO DES 4 SOMMETS
         MNN = MNN + 4
C
C        EF QUADRANGLE?
         NUS = MCN(MNN+4)
         IF( NUS .EQ. 0 ) THEN
            NSFR = -1
            GOTO 40
         ENDIF
C
         DO 15 NS=1,4
C           NUMERO DU SOMMET NS DE L'EF I
            NUS = MCN(MNN+NS)
            IF( NUS .NE. NUMPT ) THEN
C              SOMMET VU UNE FOIS DE PLUS
               MCN(MNTEMP+NUS) = MCN(MNTEMP+NUS) + 1
            ENDIF
 15      CONTINUE
C
 20   CONTINUE
C
C      SI ON EST PAS AU BORD, CHAQUE SOMMET DES EF AUTRES QUE NUMPT
C      EST APPARU DEUX ET SEULEMENT DEUX FOIS
       NSFR = 0
       DO 30 I=1,NBSOM
          NS=MCN(MNTEMP+I)
          IF( (NS.NE.0) .AND. (NS.NE.2) ) THEN
             NSFR = 1
             GOTO 40
          ENDIF
 30    CONTINUE
C
 40    CALL TNMCDS( 'ENTIER', 1+NBSOM, MNTEMP )
       RETURN
       END
