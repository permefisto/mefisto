      program modiftexte
c
c     transforme sousobjets en nsef et elements en npef
c     ==================================================== 
c                                   A. PERRONNET mars 1994
c
      character*80  nmbibl,nmftn,newftn,ligne
c         
c     le nom de la bibliotheque a traiter
      write(6,10000) nmbibl 
10000 format('nom de la bibliotheque mefisto=?',a)
      read(5,'(a)') nmbibl
c     recherche du premier blanc
      i = index( nmbibl , ' ' )  
      newftn = '//thym/users/peronet/mefisto/'
     %          // nmbibl(1:i-1) // '/f1/' 
      nmbibl = newftn
      lbibl  = 32 + i 
c    
c     ouverture du fichier liste des fichiers fortran
      open( unit=20 , file='../bibftn' , iostat=nerr )
      if( nerr .ne. 0 ) then
         write(6,10001) 'bibftn'
10001 format(' fichier non ouvrable',a ) 
         stop
      endif
c                  
c     lecture du nom du fichier fortran
 10   read(20,'(a)', end=8000 , iostat=nerr ) nmftn 
c
c     traitement si suffixe= .ftn ou .inc  
      i = index( nmftn , '.bak' )
      if( i .gt. 0 ) goto 10
c
      i = index( nmftn , '.ftn' )
      if( i .le. 0 ) then
         i = index( nmftn , '.inc' )
         if( i .le. 0 ) then
            write(6,*) 'fichier non traite: ',nmftn
            goto 10
         endif
      endif
      newftn = nmftn(1:i+4) 
      nmftn  = newftn
c      
c     ouverture de l'ancien fichier fortran
      open( unit=21 , file=nmftn ,iostat=nerr )
      if( nerr .ne. 0 ) then
         write(6,10001) nmftn
         goto 10
      endif 
      write(6,*) 'fichier ancien : ',nmftn
c     
c     ouverture du nouveau fichier fortran
      newftn = nmbibl(1:lbibl) // nmftn
      open( unit=22 , file=newftn ,iostat=nerr )
      if( nerr .ne. 0 ) then
         write(6,10001) newftn
         goto 10
      endif
      write(6,*) 'fichier nouveau: ',newftn
      isuite = 0
c
 20   read(unit=21,fmt='(a)',err=9000,end=50) ligne
c
c     MODIFICATIONS
c     =============                 
      call change ( ligne , '' , '//thym/users/peronet/mefisto' , '.' )
c
      call change ( ligne , '' , 'CONTINIT' , 'CONTRINIT' )
      call change ( ligne , '' , 'continit' , 'contrinit' )
c  
      call change(ligne,'','''CONTRAINTES''','''CONTRAINTE''')
      call change(ligne,'','_contraintes','_contrainte')  
c                                                           
      call change ( ligne , '' , 'ES SOUSOBJETS' , 'ES NO SOMMET' )
      call change ( ligne , '' , 'es sousobjets' , 'es no sommet' )
c                                                           
      call change ( ligne , '' , 'SOUSOBJETS' , 'NSEF' )
      call change ( ligne , '' , 'sousobjets' , 'nsef' ) 
c  
      call change ( ligne , '' , '''ELEMENTS''' , '''NPEF''')
      call change ( ligne , 'es elements' , 'elements' , 'npef' )
c  
      call change(ligne,'' , 'NONOEUDS' , 'NOUVNOEUD' )
      call change(ligne,'' , 'nonoeuds' , 'nouvnoeud' )
c  
      call change(ligne,'' , '''NOEUDS''' , '''XYZNOEUD''' )
      call change(ligne,'' , '_noeuds' , '_xyznoeud' )   
c  
      call change(ligne,'' , '''POINTS''' , '''XYZPOINT''' )     
      call change(ligne,'' , '_points.inc' , '_xyzpoint.inc' ) 
      call change(ligne,'' , '_points' , '_point' ) 
c  
      call change(ligne,'' , '''SOMMETS''' , '''XYZSOMMET''' )
      call change(ligne,'es sommets' , '_sommets' , '_xyzsommet' )
c  
      call change(ligne,'' , '''FONCTIONS''' , '''FONCTION''' )
      call change(ligne,'' , '_fonctions' , '_fonction' )  
c  
      call change(ligne,'' , '''LIGNES''' , '''LIGNE''' )
      call change(ligne,'' , '_lignes' , '_ligne' )      
c  
      call change(ligne,'' , '''OBJETS''' , '''OBJET''' )
      call change(ligne,'' , '_objets' , '_objet' )    
c  
      call change(ligne,'' , '''SURFACES''' , '''SURFACE''' )
      call change(ligne,'' , '_surfaces' , '_surface' )   
c  
      call change(ligne,'' , '''VOLUMES''' , '''VOLUME''' )
      call change(ligne,'' , '_volumes' , '_volume' ) 
c  
      call change(ligne,'' , '''ARETES''' , '''ARETE''' )
      call change(ligne,'' , '_aretes' , '_arete' )
c  
      call change(ligne,'' , '''FACES''' , '''FACE''' )
      call change(ligne,'' , '_faces' , '_face' )   
c
      call change(ligne,'','>CONTRAINTES','>CONTRAINTE')  
c
      call change(ligne,'','>ELEMENTS','>NPEF')
c                                                           
      call change(ligne,'' , '>NOEUDS' , '>XYZNOEUD' )
c  
      call change(ligne,'' , '>POINTS' , '>XYZPOINT' )     
c  
      call change(ligne,'' , '>SOMMETS' , '>XYZSOMMET' )
c  
      call change(ligne,'' , '>FONCTIONS' , '>FONCTION' )
c  
      call change(ligne,'' , '>LIGNES' , '>LIGNE' )
c  
      call change(ligne,'' , '>OBJETS' , '>OBJET' )
c  
      call change(ligne,'' , '>SURFACES' , '>SURFACE' )
c  
      call change(ligne,'' , '>VOLUMES' , '>VOLUME' )
c  
      call change(ligne,'' , '>ARETES' , '>ARETE' )
c  
      call change(ligne,'' , '>FACES' , '>FACE' )
c  
      call change(ligne, '' , 'BSSOB' , 'BEFOB' )    
      call change(ligne, '' , 'BSOSO' , 'BSOEF' )
      call change(ligne, '' , 'USOSS' , 'USOEF' )
      call change(ligne, '' , 'KNMOBS' , 'NMTYOB' )
c
c     recherche du dernier non blanc de la ligne pour 
c     limiter le stockage . suppression des lignes blanches
 33   do 45 i=80,1,-1
         if( ligne(i:i) .ne. ' ' ) goto 48
 45   continue
      goto 20
c
c     ecriture de la ligne 
 48   write(unit=22,fmt='(a)') ligne(1:i)
      goto 20 
c
c     fin du fichier fortran
 50   close( unit=21 )
      close( unit=22 ) 
      goto 10
c
c     fin de lecture de la liste des fichiers
 8000 close( unit=20 )  
      stop 
c
 9000 write(6,*) 'le contenu de la ligne est incorrect ',ligne
      goto 20
      end   
      subroutine change( ligne , excpt , nomanc , nomnou )
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
