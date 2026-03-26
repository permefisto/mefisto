      program ajax
c                                        
c     =============================================================
c       transforme dans les fichiers .f les noms de X11T en xvue
c     =============================================================
c                   A. PERRONNET Septembre 1995
c
      character*80  nmbibMODIFS,nmfrt,newfrt,oldfrt,ligne
c         
c     le nom de la bibliotheque a traiter
      write(6,10000) nmbibMODIFS 
10000 format('nom de la bibliotheque Mefisto=?',a)
      read(5,'(a)') nmbibMODIFS
c     recherche du premier blanc
      i = index( nmbibMODIFS , ' ' ) 
      oldfrt = '/net/ef/efd2/nef/'
     %          // nmbibMODIFS(1:i-1) // '/'
      newfrt = '/net/ef/efd2/nef/'
     %          // nmbibMODIFS(1:i-1) // '/MODIFS/ ' 
      nmbibMODIFS = newfrt
      lbibl  = index( newfrt , ' ' ) - 1
c    
c     ouverture du fichier liste des fichiers fortran
      open( unit=20 , file='LISTFIC' , iostat=nerr )
      if( nerr .ne. 0 ) then
         write(6,10001) 'LISTFIC'
10001 format(' fichier non ouvrable',a ) 
         stop
      endif    
C
C     NOMBRE DE LIGNES INSTRUCTION FORTRAN ET C
      NBLIF = 0
      NBLIC = 0
c                         
C     BOUCLE SUR LES FICHIERS FORTRAN
C     *******************************
c     lecture du nom du fichier fortran
 10   read(20,'(a)', end=8000 , iostat=nerr ) nmfrt 
c
c     traitement si suffixe= .f
      i = index( nmfrt , '.f' )
      if( i .le. 0 ) then
         write(6,*) 'fichier non traite: ',oldfrt,nmfrt
         goto 10
      else
         newfrt(1:1) = nmfrt(i+2:i+2)
         if( newfrt(1:1) .GE. '0' .AND. newfrt(1:1) .LE. '9' ) then
            write(6,*) 'fichier non traite: ',oldfrt,nmfrt
            goto 10
         endif
      endif
      newfrt = nmfrt(1:i+4) 
      nmfrt  = newfrt
c      
c     ouverture de l'ancien fichier fortran
      open( unit=21 , file=nmfrt ,iostat=nerr )
      if( nerr .ne. 0 ) then
         write(6,10001) nmfrt
         goto 10
      endif 
      write(6,*) 'fichier ancien : ',nmfrt
c     
c     ouverture du nouveau fichier fortran
      newfrt = nmbibMODIFS(1:lbibl) // nmfrt
      open( unit=22 , file=newfrt ,iostat=nerr )
      if( nerr .ne. 0 ) then
         write(6,10001) newfrt
         goto 10
      endif
      write(6,*) 'fichier nouveau: ',newfrt
c         
c     BOUCLE SUR LES LIGNES FORTRAN 
C     =============================
 20   read(unit=21,fmt='(a)',err=9000,end=50) ligne
c                        
C ATTENTION A NE PAS ENGENDRER LA RECURSIVITE !
C      SI Chaine initiale INCLUSE dans Chaine Finale 
C      --------------------------------------------- 
CCC      call change(ligne, '' , 'X11TYPE',  'XVTYPE'  )
CCC      call change(ligne, '' , 'X11TEXTE', 'XVTEXTE' )
CCC      call change(ligne, '' , 'X11TRAIT', 'XVTRAIT' )
CCC      call change(ligne, '' , 'X11TEST',  'XVTEST'  )
CCC      call change(ligne, '' , 'X11T',   'xvue'  )
CCC      call change(ligne, '' , 'x11t',   'xvue'  )
CCC      call change(ligne, '' , 'X11',    'XV' )
CCC      call change(ligne, '' , 'x11',    'xv' )
CCC      call change(ligne, '' , 'xtmgraph', 'xvuelc' )
CCC      call change(ligne, '' , 'xtmg',   'xvue'  )
c
         call change(ligne, '' , 'ETGMA ',   'BTGEF '  )
         call change(ligne, '' , 'LBTGEF',   'NBTGEF '  )
C
C     recherche du dernier non blanc de la ligne pour 
C     limiter le stockage . suppression des lignes blanches
      do 45 i=80,1,-1
         if( ligne(i:i) .ne. ' ' ) goto 48
 45   continue
      goto 20
c
c     ecriture de la ligne 
 48   write(unit=22,fmt='(a)') ligne(1:i) 
      if( ligne(1:1) .eq. 'C' .or. 
     %    ligne(1:1) .eq. 'c' ) then
         NBLIC = NBLIC + 1
      else
         NBLIF = NBLIF + 1
      endif
      goto 20 
c
c     fin du fichier fortran
 50   close( unit=21 )
      close( unit=22 ) 
      goto 10
c
c     fin de lecture de la liste des fichiers
 8000 close( unit=20 )  
      print *,'NOMBRE DE LIGNES FORTRAN     =',NBLIF
      print *,'NOMBRE DE LIGNES COMMENTAIRE =',NBLIC
      print *,'NOMBRE TOTAL DE LIGNES       =',NBLIC+NBLIF
      stop 
c
 9000 write(6,*) 'le contenu de la ligne est incorrect ',ligne
      write(6,*) 'dans ',nmfrt
      goto 20
      end     
C
      subroutine CHANGE( ligne , excpt , nomanc , nomnou )
      character*(*) excpt,nomanc,nomnou
      character*80  ligne,ligne1 
C      
C     si exception 0 modification 
      if( excpt .ne. '' ) then
         k = index( ligne , excpt )
         if( k .gt. 0 ) return  
      endif
C             
 10   k = index( ligne , nomanc )
      if( k .gt. 0 ) then      
         l = len( nomanc )
         ligne1 = ligne(1:k-1) // nomnou // ligne(k+l:80)
         ligne  = ligne1 
         goto 10
      endif 
      end
C
      SUBROUTINE MAJUSC( CHAINE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER LES EVENTUELLES LETTRES MINUSCULES DE CHAINE
C ----- EN MAJUSCULES EN PASSANT PAR LE CODAGE DE LA TABLE ASCII
C
C ENTREE ET SORTIE :
C ------------------
C CHAINE : LA CHAINE DE CARACTERES A TRAITER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      CHARACTER*(*) CHAINE
      CHARACTER*1   LETMAJ
C
C     BOUCLE SUR LES CARACTERES
      DO 10 I = 1 , LEN( CHAINE )
         CHAINE(I:I) = LETMAJ( CHAINE(I:I) )
 10   CONTINUE
      END 
      CHARACTER*1 FUNCTION LETMAJ( CAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER L'EVENTUELLE LETTRE MINUSCULE CAR EN SA MAJUSCULE
C ----- EN PASSANT PAR LE CODAGE DE LA TABLE ASCII
C
C       CETTE FUNCTION FAIT APPEL AUX FONCTIONS CHAR ET ICHAR DE FORTRAN
C
C ENTREE :
C --------
C CAR    : LE CARACTERE A TRAITER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      CHARACTER*(*) CAR
      CHARACTER*1   CHAR
C
C     LA PLACE DU CARACTERE DANS LA TABLE ASCII
      N = ICHAR( CAR(1:1) )
C
C     CE CARACTERE EST IL UNE LETTRE MINUSCULE ?
      IF( N .GE. 97 .AND. N .LE. 122 ) THEN
C        OUI: LA MAJUSCULE CORRESPONDANTE
         LETMAJ = CHAR( N - 32 )
      ELSE
C        NON: PAS DE MODIFICATION
         LETMAJ = CAR(1:1)
      ENDIF
      END  
      CHARACTER*1 FUNCTION LETMIN( CAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER L'EVENTUELLE LETTRE MAJUSCULE CAR EN SA MINUSCULE
C ----- EN PASSANT PAR LE CODAGE DE LA TABLE ASCII
C
C       CETTE FUNCTION FAIT APPEL AUX FONCTIONS CHAR ET ICHAR DE FORTRAN
C
C ENTREE :
C --------
C CAR    : LE CARACTERE A TRAITER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      CHARACTER*(*) CAR
      CHARACTER*1   CHAR
C
C     LA PLACE DU CARACTERE DANS LA TABLE ASCII
      N = ICHAR( CAR(1:1) )
C
C     CE CARACTERE EST IL UNE LETTRE MAJUSCULE ?
      IF( N .GE. 65 .AND. N .LE. 90 ) THEN
C        OUI: LA MAJUSCULE CORRESPONDANTE
         LETMIN = CHAR( N + 32 )
      ELSE
C        NON: PAS DE MODIFICATION
         LETMIN = CAR(1:1)
      ENDIF
      END

