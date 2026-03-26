      SUBROUTINE VEOPTE( NBTETRA, NOTETRA, NOTETR, PTXYZD, NBFANR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     VERIFIER L'ABSENCE DE TETRAEDRE  OPPOSE  NEGATIF (INCONNU)
C -----     VERIFIER L'ABSENCE DE TETRAEDRES OPPOSES DOUBLES
C           VERIFIER L'ABSENCE DE FACE COMMUNE A AU MOINS 3 TETRAEDRES
C           VERIFIER L'OPPOSITION DES TETRAEDRES 
C           DANS LA LISTE DES NBTETRA TETRAEDRES NOTETRA

C ENTREES:
C --------
C NBTETRA : NOMBRE DE TETRAEDRES DE LA LISTE NOTETRA
C NOTETRA : NUMERO NOTETR DES TETRAEDRES DE LA LISTE
C NOTETR  : LISTE DES TETRAEDRES
C           SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C           TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C           DE L'AUTRE COTE DE LA FACE
C           1: 123      2: 234      3: 341      4: 412
C PTXYZD  : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE

C SORTIE  :
C ---------
C NBFANR  : NOMBRE DE FACES DES TETRAEDRES DE TETRAEDRE OPPOSE A PROBLEME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2016
C2345X7..............................................................012
      PARAMETER        (MXTEVN=2048)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      CHARACTER*100     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*), VOLTET, VOLNTE, VOLNTOP, V
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)

      INTEGER  NOTETRA(NBTETRA), NOTETR(8,*), NOFANR(2,MXTEVN),
     %         NOTEVN(MXTEVN), NOSOTR(3), NOTEDB(10), NO3TEF(3)
      INTEGER  NOSOFATE(3,4)
      DATA     NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      TRACTE0= TRACTE
      NBFANR = 0
      NBTEVN = 0

C     RECHERCHE DE TETRAEDRES OPPOSES DOUBLES POUR UN TETRAEDRE
C     ---------------------------------------------------------
      NBTEDO = 0
      DO NT = 1, NBTETRA

C        NO NOTETR DU TETRAEDRE NT SUPPOSE >0
         NTE = NOTETRA( NT )
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

C           VERIFICATION: UN TETRAEDRE OPPOSE EST IL DOUBLE?
            DO 5 NF1=1,4

               NTOP = NOTETR( 4+NF1, NTE )
               IF( NTOP .LT. 0 ) THEN

C                 NO DE TETRAEDRE OPPOSE <0 => TETRAEDRE INCONNU
                  CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                          PTXYZD(1,NOTETR(2,NTE)),
     %                          PTXYZD(1,NOTETR(3,NTE)),
     %                          PTXYZD(1,NOTETR(4,NTE)),
     %                          ARMIN, ARMAX, SURFTR, VOLNTE, QUALTE )
                  PRINT*,'veopte: ATTENTION un TETRAEDRE OPPOSE=',NTOP,
     %                  ' du TETRAEDRE',NTE,':',(NOTETR(kk,NTE),kk=1,8),
     %                   ' V=',VOLNTE,' Q=',QUALTE
                  GOTO 5

               ENDIF

               IF( NTOP .GT. 0 ) THEN
                  DO NF2=NF1+1,4
                     NTOP2 = NOTETR( 4+NF2, NTE )
                     IF( NTOP .EQ. NTOP2 ) THEN

C                       UN TETRAEDRE OPPOSE DOUBLE POUR NTE
                        NBTEDO = NBTEDO + 1

                        CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                                PTXYZD(1,NOTETR(2,NTE)),
     %                                PTXYZD(1,NOTETR(3,NTE)),
     %                                PTXYZD(1,NOTETR(4,NTE)),
     %                         ARMIN, ARMAX, SURFTR, VOLNTE, QUALTE )

                        CALL QUATETD( PTXYZD(1,NOTETR(1,NTOP)),
     %                                PTXYZD(1,NOTETR(2,NTOP)),
     %                                PTXYZD(1,NOTETR(3,NTOP)),
     %                                PTXYZD(1,NOTETR(4,NTOP)),
     %                         ARMIN, ARMAX, SURFTR, VOLNTOP, QUALTOP )

                        PRINT*
                        NOTEDB(1) = NTE
                        PRINT*,'veopte: ATTENTION TETRAEDRE         ',
     %                          NTE,':',(NOTETR(kk,NTE),kk=1,8),
     %                         ' V=',VOLNTE,' Q=',QUALTE

                        NOTEDB(2) = NTOP
                        NBTEDB = 2
                        PRINT*,'veopte: AVEC DOUBLE TETRAEDRE OPPOSE',
     %                          NTOP,':',(NOTETR(kk,NTOP),kk=1,8),
     %                         ' V=',VOLNTOP,' Q=',QUALTOP

                        PRINT*
                        DO kkk=1,4
                           NTT = NOTETR( 4+kkk, NTE )
                           IF( NTT .GT. 0 ) THEN
                              NBTEDB = NBTEDB + 1
                              NOTEDB(NBTEDB) = NTT
                           ENDIF
                        ENDDO

                        PRINT*
                        DO kkk=1,4
                           NTT = NOTETR( 4+kkk, NTOP )
                           IF( NTT .GT. 0 ) THEN
                              NBTEDB = NBTEDB + 1
                              NOTEDB(NBTEDB) = NTT
                           ENDIF
                        ENDDO

C                       SUPPRESSION DES DOUBLONS
                        NBTD = 2
                        DO 7 KKK=3,NBTEDB
                           NTT = NOTEDB(KKK)
                           DO NN=1,NBTD
                              IF( NTT .EQ. NOTEDB(NN) ) GOTO 7
                           ENDDO
                           NBTD = NBTD + 1
                           NOTEDB( NBTD ) = NTT
 7                      ENDDO
                        NBTEDB = NBTD

                        PRINT*
                        DO kkk=1,NBTEDB
                           NTT = NOTEDB(kkk)
                           CALL QUATETD( PTXYZD(1,NOTETR(1,NTT)),
     %                                   PTXYZD(1,NOTETR(2,NTT)),
     %                                   PTXYZD(1,NOTETR(3,NTT)),
     %                                   PTXYZD(1,NOTETR(4,NTT)),
     %                                   ARMIN, ARMAX, SURFTR, V, Q )
                           PRINT*,'veopte:',NTT,':',
     %                           (NOTETR(kk,NTT),kk=1,8),' V=',V,' Q=',Q
                        ENDDO

C                       TRACE DES TETRAEDRES NTE NTOP et OPPOSES
ccc                        tracte = .true.
                 KTITRE='TRACE des TETRAEDRES NTE NTOP et leurs OPPOSES'
                        CALL TRFETO15( KTITRE, PTXYZD,
     %                                 NBTEDB, NOTEDB, NOTETR, 0, 0 )

                     ENDIF
                  ENDDO
               ENDIF
 5          ENDDO

         ENDIF

      ENDDO

      IF( NBTEDO .GT. 0 ) THEN
ccc         tracte = .true.
      KTITRE='veopte:           TETRAEDRES AVEC un TETRAEDRE OPPOSE A 2 
     %de SES FACES'
         WRITE(KTITRE(9:11),'(I3)') NBTEDO
         CALL SANSDBL( KTITRE, NBC )
         CALL TRAFNBTE( KTITRE(1:NBC), PTXYZD, NBTETRA, NOTETRA, NOTETR)
         print*
      ENDIF


C     VERIFICATION DES TETRAEDRES OPPOSES AUX FACES
C     ---------------------------------------------
      DO NT = 1, NBTETRA

C        NO NOTETR DU TETRAEDRE NT SUPPOSE >0
         NTE = NOTETRA( NT )
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

         DO 10 NF1=1,4

C           NO DU TETRAEDRE OPPOSE A LA FACE NF1 DE NTE
            NTEOP = NOTETR(4+NF1,NTE)
            IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN

               DO M=1,3
C                 NUMERO DU SOMMET M DE LA FACE NF1 DE NTE
                  NOSOTR(M) = NOTETR( NOSOFATE(M,NF1), NTE )
               ENDDO

C              LA FACE NF1 DE NTE NOSOTR EST ELLE UNE FACE NFTEOP DE NTEOP?
               CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP), NFTEOP )

               IF( NFTEOP .GT. 0 ) THEN

C                 OUI: LA FACE NF1 DE NTE EST LA FACE NFTEOP DE NTEOP
C                 LE TETRAEDRE OPPOSE A LA FACE NFTEOP DE NTEOP EST IL BIEN NTE?
                  NTE2 = NOTETR( 4+NFTEOP, NTEOP )
                  IF( NTE2 .EQ. NTE ) THEN

C                    OUI: LES 2 TETRAEDRES NTEOP NTE2 SONT BIEN OPPOSES
                     GOTO 10

                  ELSE

C                    NON: PROBLEME: UNE FACE APPARTIENT A AU MOINS 3 TETRAEDRES
                     PRINT*
                     PRINT*,'veopte: PROBLEME le TETRAEDRE OPPOSE',
     %                        NTEOP,' pour sa FACE',NFTEOP,
     %                       'N''EST PAS NTE=',NTE,' MAIS NTE2=',NTE2
                     PRINT*,'veopte: NOTETR(',NTE,  ')=',
     %                       (NOTETR(N,NTE),N=1,8)
                     PRINT*,'veopte: NOTETR(',NTEOP,')=',
     %                       (NOTETR(N,NTEOP),N=1,8)
                     PRINT*,'veopte: NOTETR(',NTE2,')=',
     %                       (NOTETR(N,NTE2),N=1,8)

C                    NTE2 A T IL BIEN LA FACE NOSOTR?
                     CALL NUFATRTE( NOSOTR, NOTETR(1,NTE2), NFTE2 )
                     IF( NFTE2 .GT. 0 ) THEN
                      PRINT*,'veopte: PB la FACE',NOSOTR,' EST la FACE',
     %                          NFTE2,' DE',NTE2,
     %                         ' => ELLE APPARTIENT aux 3 TETRAEDRES',
     %                          NTE,NTEOP,NTE2
                        print*
                     ENDIF

      KTITRE='veopte: LA FACE                             APPARTIENT A 3
     % TETRAEDRES'
                     WRITE(KTITRE(17:24),'(I8)') NOSOTR(1)
                     WRITE(KTITRE(26:33),'(I8)') NOSOTR(2)
                     WRITE(KTITRE(35:42),'(I8)') NOSOTR(3)
                     CALL SANSDBL( KTITRE, NBC )
                     NB3TEF = 3
                     NO3TEF( 1 ) = NTE
                     NO3TEF( 2 ) = NTEOP
                     NO3TEF( 3 ) = NTE2
                     tracte = .true.
                     CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, 0,
     %                                 NB3TEF, NO3TEF )

                     IF( NBFANR .LT. MXTEVN ) THEN
                        NBFANR = NBFANR + 1
C                       NO NOTETR DU TETRAEDRE DE FACE SANS OPPOSE CORRECT
                        NOFANR( 1, NBFANR ) = NTE
C                       NO DE LA FACE DANS LE TETRAEDRE
                        NOFANR( 2, NBFANR ) = NF1
                     ENDIF

                  ENDIF

               ELSE

C                 NTEOP N'A PAS LA FACE NF1 DE NTE
ccc                  tracte = .true.
                  PRINT*
                  PRINT *,'veopte: Pb LA FACE',NF1, ' de NTE=',NTE,
     %                    ' N''EST PAS UNE FACE de NTEOP=',NTEOP
                  PRINT *,'veopte: NOTETR(',NTE,  ')=',
     %                     (NOTETR(N,NTE),N=1,8)
                  PRINT *,'veopte: NOTETR(',NTEOP,')=',
     %                     (NOTETR(N,NTEOP),N=1,8)

                  IF( NBFANR .LT. MXTEVN ) THEN
                     NBFANR = NBFANR + 1
C                    NUMERO NOTETR DU TETRAEDRE DE FACE SANS OPPOSE
                     NOFANR( 1, NBFANR ) = NTE
C                    NUMERO DE LA FACE DANS LE TETRAEDRE
                     NOFANR( 2, NBFANR ) = NF1
                  ENDIF

      KTITRE='veopte: LA FACE                             N''APPARTIENT 
     %AU TETRAEDRE OPPOSE'
                  WRITE(KTITRE(17:24),'(I8)') NOSOTR(1)
                  WRITE(KTITRE(26:33),'(I8)') NOSOTR(2)
                  WRITE(KTITRE(35:42),'(I8)') NOSOTR(3)
                  CALL SANSDBL( KTITRE, NBC )
                  NB3TEF = 2
                  NO3TEF( 1 ) = NTE
                  NO3TEF( 2 ) = NTEOP
                  PRINT*,'veopte:',KTITRE(1:NBC)
                  tracte = .true.
                  CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, 0,
     %                              NB3TEF, NO3TEF )
               ENDIF

            ELSE IF( NTEOP .LT. 0 ) THEN

C              NTEOP=-1 FACE AVEC TETRAEDRE OPPOSE NON RETROUVE
ccc               tracte = .true.
               PRINT*
               PRINT*,'veopte: LA FACE',NF1, ' de NTE=',NTE,
     %                ' N''A PAS de TETRAEDRE OPPOSE NTEOP=',NTEOP
               PRINT*,'veopte: NOTETR(',NTE,  ')=',(NOTETR(N,NTE),N=1,8)

               IF( NBFANR .LT. MXTEVN ) THEN
                  NBFANR = NBFANR + 1
C                 NUMERO NOTETR DU TETRAEDRE DE FACE SANS OPPOSE
                  NOFANR( 1, NBFANR ) = NTE
C                 NUMERO DE LA FACE DANS LE TETRAEDRE
                  NOFANR( 2, NBFANR ) = NF1
               ENDIF

ccc         ELSE
cccC           NTEOP= 0 FACE FRONTALIERE

            ENDIF

 10      ENDDO
         ENDIF

      ENDDO

      IF( NBFANR .GT. 0 ) THEN

C        AFFICHAGE DES FACES NON RETROUVEES
C        ----------------------------------
         PRINT *,'veopte:',NBFANR,' FACES NON RETROUVEES dans les',
     %            NBTETRA,' TETRAEDRES'
         PRINT *

C        RECHERCHE DES TETRAEDRES DE VOLUME NEGATIF
C        ------------------------------------------
         NBTEVN = 0
         DO NT1=1,NBTETRA
C           NO NOTETR DU TETRAEDRE NT1
            NTE1 = NOTETRA( NT1 )
            IF( NTE1 .GT. 0 .AND. NOTETR(1,NTE1) .GT. 0 ) THEN
               VOLNTE = VOLTET( PTXYZD(1,NOTETR(1,NTE1)),
     %                          PTXYZD(1,NOTETR(2,NTE1)),
     %                          PTXYZD(1,NOTETR(3,NTE1)),
     %                          PTXYZD(1,NOTETR(4,NTE1)) )
               IF( VOLNTE .LT. 0D0 ) THEN
C                 STOCKAGE DU NO DE TETRAEDRE A VOLUME NEGATIF
                  NBTEVN = NBTEVN + 1
                  NOTEVN( NBTEVN ) = NTE1
               ENDIF
ccc               PRINT *,'veopte: ',NT1,' NOTETR(',NTE1,')=',
ccc     %                 (NOTETR(L,NTE1),L=1,8),' Volume=',VOLNTE
            ENDIF
         ENDDO

C        TRACE DES FACES NON RETROUVEES
         TRACTE = .TRUE.
         KTITRE='veopte:       FACES de TETRAEDRE OPPOSE NON RETROUVE'
         WRITE( KTITRE(8:12),'(I5)') NBFANR
         CALL SANSDBL( KTITRE, NBC )
         PRINT*,'veopte: ',KTITRE(1:NBC)
         CALL TRFETO12( KTITRE,  PTXYZD,
     %                  NBTETRA, NOTETRA, NOTETR, NBFANR, NOFANR )

         IF( NBTEVN .GT. 0 ) THEN
            KTITRE='veopte:       TETRAEDRES DE VOLUME <0'
            WRITE( KTITRE(8:12),'(I5)') NBTEVN
            CALL TRFETO13( KTITRE, PTXYZD, NBTEVN, NOTEVN, NOTETR )
         ENDIF

      ENDIF

      TRACTE = TRACTE0
      RETURN
      END

