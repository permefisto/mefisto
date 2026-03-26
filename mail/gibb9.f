      SUBROUTINE GIBB9(NBNOE,NSND,NUM,ISDIR,NFLG,IDPTH,MAXDEG,
     %                 LISTVOI,LPVOIS,IALVL2,IANDEG,
     %                 IAIPFA,IALSTP,IALVLS,IASTKA,IASTKB,IASTKC,IASTKD,
     %                 MNRENU,LBD,LPRF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     RENUMEROTER LES NOEUDS SUIVANT LA STRUCTURE DE NIVEAU LVL2
C ----

C ENTREES:
C --------
C NBNOE  : NOMBRE DE NOEUDS
C NSND   : EXTREMITE DE PSEUDO-DIAMETRE DE + BAS DEGRE TROUVEE DANS GIBB3
C NUM    : NOUVEAU NUMERO DU PREMIER ELEMENT DE LVL2
C ISDIR  : SENS DE PARCOURS DE LA NUMEROTATION
C NFLG   : SENS DE PARCOURS DES NIVEAUX DE LVL2
C IDPTH  : HAUTEUR DE LA STRUCTURE DE NIVEAU LVL2
C MAXDEG : DEGRE MAXIMUM DES NOEUDS1
C LISTVOI: LISTE DES VOISINS DES NOEUDS
C LPVOIS : POINTEUR DERNIER NOEUD VOISIN DE CHAQUE NOEUD
C IALVL2 : ADRESSE MCN DE LA STRUCTURE DE NIVEAU TROUVEE DANS GIBB8
C IANDEG : ADRESSE MCN DE TABLEAU DES DEGRES DE CHAQUE NOEUD

C MODIFIEES:
C ----------
C IAIPFA : ADRESSE MCN DU TABLEAU DES DISTANCES MAXIMALES ENTRE UN NOEUD
C          ET SES VOISINS DANS LA NOUVELLE NUMEROTATION
C IALSTP : ADRESSE MCN DU TABLEAU DES ADRESSES DANS LVLS DU PREMIER NOEUD
C          DE CHAQUE NIVEAU
C IALVLS : ADRESSE MCN DU TABLEAU DES NOEUDS CLASSES PAR NIVEAU
C IASTKA : \
C IASTKB : \ ADRESSES MCN DE DIVERS TABLEAUX INTERMEDIAIRES
C IASTKC : /
C IASTKD : /

C SORTIES:
C --------
C MNRENU : ADRESSE MCN DU TABLEAU CONTENANT LA NOUVELLE NUMEROTATION
C LBD    : LA NOUVELLE LARGEUR DE BANDE
C LPRF   : LE NOUVEAU PROFIL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR   : BARGACH MOHAMED LAN189 PARIS             OCTOBRE  1980
C MODIFICATIONS : DEFAIX THIERRY                           DECEMBRE 1989
C MODIFICATIONS : Alain PERRONNET Saint Pierre du Perray    Janvier 2021
C23456---------------------------------------------------------------012
      INTEGER   LPVOIS(0:NBNOE), LISTVOI(*),
     %          NXA,NXB,NXC,NXD,NTEST,NSND,INCX
      include"./incl/pp.inc"
      COMMON    MCN(MOTMCN)

C     INITIALISATION ET CONSTRUCTION DE LVLS ET LSTP
C     ==============================================
      IANDEG1 = IANDEG - 1
      IALVLS1 = IALVLS - 1
      MNRENU1 = MNRENU - 1
      DO I=1,NBNOE
         MCN(IAIPFA-1+I)=0
      ENDDO

      NSTPT=1
      DO I=1,IDPTH
         MCN(IALSTP-1+I)=NSTPT
         DO J=1,NBNOE
            IF(MCN(IALVL2-1+J).EQ.I) THEN
               MCN(IALVLS1+NSTPT)=J
               NSTPT=NSTPT+1
            ENDIF
         ENDDO
      ENDDO
      MCN(IALSTP-1+IDPTH+1)=NSTPT

C     NUMEROTATION
C     ============
C     -LVLN=indice de parcours des niveaux
      LVLN=0
      IF(NFLG.LT.0) LVLN=IDPTH+1
C     -NXC indice d'ecriture dans STKC
      NXC=1
C     -On classe la racine dans STKC
      MCN(IASTKC-1+NXC)=NSND

C     BOUCLE SUR LES NIVEAUX
C     ----------------------
C     -On passe au niveau suivant
C     --INCX indice de parcours des noeuds deja classes dans STKC
 10   INCX=1
      NXD=0
      LVLN=LVLN+NFLG
      LST=MCN(IALSTP-1+LVLN)
      LND=MCN(IALSTP-1+LVLN+1)-1

C     BOUCLE SUR LES POINTS D'UN NIVEAU
C     .................................
C     -On passe au noeud suivant
 20   IPRO=MCN(IASTKC-1+INCX)
      MCN(MNRENU1+IPRO)=NUM
      NUM=NUM+ISDIR
      NXA=0
      NXB=0
ccc      MU1=MCN(MNLPVO-1+IPRO) + 1
ccc      MU2=MCN(MNLPVO-1+IPRO+1)
      MU1 = LPVOIS(IPRO-1) + 1
      MU2 = LPVOIS(IPRO  )

C     -On traite les voisins des noeuds deja classes
      DO J=MU1,MU2
ccc         NTEST=MCN(MNLIVO-1+J)  modif 2021/01/19 pour ALLOCATE
         NTEST=LISTVOI(J)
         INX=MCN(MNRENU1+NTEST)
         IF(INX.EQ.0) THEN
            MCN(MNRENU1+NTEST)=-1
            IF(MCN(IALVL2-1+NTEST).EQ.MCN(IALVL2-1+IPRO)) THEN
               NXA=NXA+1
               MCN(IASTKA-1+NXA)=NTEST
             ELSE
               NXB=NXB+1
               MCN(IASTKB-1+NXB)=NTEST
            ENDIF
         ENDIF
         IF(INX.GT.0) THEN
            NBW=(MCN(MNRENU1+IPRO)-INX)*ISDIR+1
            IF(ISDIR.GT.0) INX=MCN(MNRENU1+IPRO)
            IF(MCN(IAIPFA-1+INX).LT.NBW) MCN(IAIPFA-1+INX)=NBW
         ENDIF
      ENDDO

C     -On classe les noeuds qui sont dans le niveau courant
      IF(NXA.GT.0) THEN
         IF(NXA.EQ.1) THEN
            NXC=NXC+1
            MCN(IASTKC-1+NXC)=MCN(IASTKA-1+NXA)
          ELSE
            CALL GIBB5( MCN(IANDEG),NXC,NXA,MCN(IASTKA),
     %                  MCN(IASTKC))
         ENDIF
      ENDIF

      IF(NXB.GT.0) THEN
         IF(NXB.EQ.1) THEN
            NXD=NXD+1
            MCN(IASTKD-1+NXD)=MCN(IASTKB-1+NXB)
          ELSE
            CALL GIBB5( MCN(IANDEG),NXD,NXB,MCN(IASTKB),
     %                  MCN(IASTKD))
         ENDIF
      ENDIF
      INCX=INCX+1
      IF( NXC .GE. INCX ) GOTO 20

C     -On cherche parmi les noeuds du niveau courant non encore numerotes
C     -celui de plus bas degre et on le met dans NSND
      NMAX=MAXDEG+1
C     NSND VALEUR TEST DE NON MODIFICATION
      NSND=NBNOE+1
      DO 30 I=LST,LND
         NTEST=MCN(IALVLS1+I)
         IF( MCN(MNRENU1+NTEST) .NE. 0    ) GOTO 30
         IF( MCN(IANDEG1+NTEST) .GE. NMAX ) GOTO 30
         MCN(MNRENU1+NSND)=0
         MCN(MNRENU1+NTEST)=-1
         NMAX=MCN(IANDEG1+NTEST)
         NSND=NTEST
 30   ENDDO

      IF( NSND .NE. NBNOE+1 ) THEN
C     ---Le noeud NSND existe
         NXC=NXC+1
         MCN(IASTKC-1+NXC)=NSND
         GOTO 20
      ENDIF

      IF( NXD .NE. 0 ) THEN
C     ---Existence d'un niveau superieur
         DO I=1,NXD
            MCN(IASTKC-1+I)=MCN(IASTKD-1+I)
         ENDDO
         NXC=NXD
         GOTO 10
      ENDIF

C     FIN DE LA NUMEROTATION
C     ======================
      DO I=1,NBNOE
         IF(MCN(IAIPFA-1+I).GT.LBD) LBD=MCN(IAIPFA-1+I)
         LPRF=LPRF+MCN(IAIPFA-1+I)
      ENDDO

      RETURN
      END
