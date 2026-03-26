      SUBROUTINE VERIFNS( TEXT, NBSOMM, NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %                    NBTECNS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    VERIFIER QUE TOUS LES TETRAEDRES CHAINES DU SOMMET NS
C -----    CONTIENNENT LE SOMMET NS
C
C ENTREES:
C --------
C TEXT   : POUR LOCALISER L'APPEL
C NBSOMM : NOMBRE DE SOMMETS A VERIFIER
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NOSOTE(>3)
C NOSOTE : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NOSOTE DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C
C SORTIES:
C --------
C NBTECNS: NOMBRE DE TETRAEDRES CHAINES CONTENANT LE SOMMET NS
C IERR   : 0 SI NS DANS TOUS LES TETRAEDRES CHAINES
C          1 SI AU MOINS UN TETRAEDRE CHAINE NE CONTIENT PAS NS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & St Pierre du PerrayAout 2012
C....................................................................012
      CHARACTER*(*)     TEXT
      INTEGER           NOSOTE(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*)
C
      IERR = 0
      DO 100 NS=MAX(1,NBSOMM-10),NBSOMM

         NBTECNS = 0
C
C        POSITION DANS NOTESO DU 1-ER TETRAEDRE DE SOMMET NS
         NDT = NO1TSO( NS )
C
C        TANT QU'IL EXISTE UN TETRAEDRE CHAINE DE SOMMET NS FAIRE
 10      IF( NDT .GT. 0 ) THEN
C
C           LE NUMERO DU TETRAEDRE DANS NOSOTE
            NT = NOTESO(1,NDT)
C
C           VERIFICATION DU CHAINAGE DES TETRAEDRES DU SOMMET NS
            DO J=1,4
C
               IF( NOSOTE(J,NT) .LE.  0 ) THEN
                  print *
                  print *,'verifns:',text,' probleme nosote(',j,',',
     %            nt,')=',nosote(j,nt),' parmi',(nosote(kk,nt),kk=1,4)
                  print *,'entrez un entier pour continuer'
                  read *,ierr
                  print*,'valeur lue=',ierr
               ENDIF
C
               IF( NOSOTE(J,NT) .EQ. NS ) GOTO 20
            ENDDO
c
C           LE SOMMET NS N'APPARTIENT PAS AU TETRAEDRE NT
            print *
            print *,'verifns:',text,' probleme ns=',ns,
     %      ' non sommet de nosote(',NT,')=(',(nosote(kk,nt),kk=1,4)
            print *,'les tetraedres chaines de sommet ns=',ns,' sont'
            ndt0 = no1tso( ns )
 17         if( ndt0 .gt. 0 ) then
c           le numero du tetraedre dans nosote
               ntt = noteso(1,ndt0)
               print *,'nosote(',ntt,')=',(nosote(kk,ntt),kk=1,4)
c              le tetraedre suivant
               ndt0 = noteso(2,ndt0)
               goto 17
            endif
            print *,'entrez un entier'
            read *,ierr
            print*,'valeur lue=',ierr
            IERR = 1
C
C           UN TETRAEDRE DE PLUS CONTIENT LE SOMMET NS
 20         NBTECNS = NBTECNS + 1
C
C           LE TETRAEDRE SUIVANT
            NDT = NOTESO(2,NDT)
            GOTO 10
         ENDIF
C
 100  CONTINUE
C
      RETURN
      END
