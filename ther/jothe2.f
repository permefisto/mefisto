      SUBROUTINE JOTHE2( NUOBSD, NUOBJE , NBJOIN , NBNOEU ,
     &                   JOIM ,  NBJOIM , JOIE ,   NBJOIE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DEFINIR POUR CHAQUE SOUS-DOMAINE LES JOINTS MAITRES ET
C ----- LES JOINTS ESCLAVES. LA STRATEGIE EST SDi EST MAITRE SUR
C       SDj SI LE NOMBRE DE NOEUDS DU COTE SDi EST PLUS GRAND QUE
C       DU COTE SDj, AVEC SDi ET SDj SOUS-DOMAINES VOISINS.
C
C ENTREES :
C ---------
C NOUBSD : LE NUMERO D'OBJET DU SOUS-DOMAINE TRAITE
C NUOBJE : LE TABLEAU DES NUMEROS D'OBJETS
C NBJOIN : NOMBRE DE JOINTS
C NBNOEU : LE NOMBRE DE NOEUDS SUR LE JOINT
C
C SORTIES :
C ---------
C JOIM   : LE TABLEAU DES JOINTS MAITRES ASSOCIES AU SOUS-DOMAINE
C NBJOIM : LE NOMBRE DE JOINTS MAITRES
C JOIE   : LE TABLEAU DES JOINTS ESCLAVES ASSOCIES AU SOUS-DOMAINE
C NBJOIE : LE NOMBRE DE JOINTS ESCLAVES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MARS 1993
C23456---------------------------------------------------------------012
      INTEGER  NUOBJE(NBJOIN,4),NBNOEU(NBJOIN,2),
     &         JOIM(NBJOIN),JOIE(NBJOIN)
      CHARACTER*10    NMTYOB,KNOMTY
      CHARACTER*24    KNOMOB,KNOMI,KNOMJ
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSIMTA / NOIMPR
C
      NBJOIM = 0
      NBJOIE = 0
      DO 100 NJ=1,NBJOIN
C        LES NUMEROS D'OBJET DE 2 SOUS-DOMAINES VOISINS
         NSDI = NUOBJE(NJ,2)
         NSDJ = NUOBJE(NJ,4)
         NBNI = NBNOEU(NJ,1)
         NBNJ = NBNOEU(NJ,2)
         NUTY = 5
         KNOMTY = NMTYOB(NUTY)
         CALL NMOBNU( KNOMTY , NSDI , KNOMI )
         CALL NMOBNU( KNOMTY , NSDJ , KNOMJ )
         KNOMOB = KNOMI
         I = INDEX( KNOMOB , '_' )
         KNOMI = KNOMOB(1:I-1)
         KNOMOB = KNOMJ
         I = INDEX( KNOMOB , '_' )
         KNOMJ = KNOMOB(1:I-1)
CW         WRITE (IMPRIM,4421) NJ,KNOMI,KNOMJ
         IF (NUOBSD.EQ.NSDI) THEN
            IF (NBNI.GT.NBNJ) THEN
C              JOINT MAITRE
               NBJOIM = NBJOIM + 1
               JOIM(NBJOIM) = NJ
            ELSE IF (NBNI.LT.NBNJ) THEN
C              JOINT ESCLAVE
               NBJOIE = NBJOIE + 1
               JOIE(NBJOIE) = NJ
            ELSE
C           EN CAS D'EGALITE : CHOIX D'APRES LE NUMERO DE SD
               IF (NSDI.LT.NSDJ) THEN
C                 JOINT MAITRE
                  NBJOIM = NBJOIM + 1
                  JOIM(NBJOIM) = NJ
               ELSE
C                 JOINT ESCLAVE
                  NBJOIE = NBJOIE + 1
                  JOIE(NBJOIE) = NJ
               END IF
            END IF
         ELSE IF (NUOBSD.EQ.NSDJ) THEN
            IF (NBNJ.GT.NBNI) THEN
C              JOINT MAITRE
               NBJOIM = NBJOIM + 1
               JOIM(NBJOIM) = NJ
            ELSE IF (NBNJ.LT.NBNI) THEN
C              JOINT ESCLAVE
               NBJOIE = NBJOIE + 1
               JOIE(NBJOIE) = NJ
            ELSE
C           EN CAS D'EGALITE : CHOIX D'APRES LE NUMERO DE SD
               IF (NSDJ.LT.NSDI) THEN
C                 JOINT MAITRE
                  NBJOIM = NBJOIM + 1
                  JOIM(NBJOIM) = NJ
               ELSE
C                 JOINT ESCLAVE
                  NBJOIE = NBJOIE + 1
                  JOIE(NBJOIE) = NJ
               END IF
            END IF
         END IF
 100  CONTINUE
C
C================================================== IMPRESSIONS  ==============
      IF (NOIMPR.GE.10)  THEN
      WRITE(IMPRIM,4422) NBJOIM,NBJOIE
      WRITE(IMPRIM,4423) (JOIM(K),K=1,NBJOIM)
      WRITE(IMPRIM,4424) (JOIE(K),K=1,NBJOIE)
CW 4421 FORMAT(' JOINT :',I6,/A80,/A80)
 4422 FORMAT( ' NOMBRE DE JOINTS MAITRES',I5,' ESCLAVES',I5)
 4423 FORMAT(' JOINTS MAITRES ',8I5)
 4424 FORMAT(' JOINTS ESCLAVES',8I5)
      ENDIF
C==============================================================================
C
      RETURN
      END
