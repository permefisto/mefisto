      program entierw

c     ===================================================== 
c     ajoute aux fichiers ~/incl/a_*.inc une instruction de
c     declaration des variables entieres debutant par W
c     ===================================================== 
c                            Alain PERRONNET  Novembre 1994

      character*80  nmftn,newftn
      character*80  ligne(100),iligne(100)
    
c     le nom de la bibliotheque a traiter
      newftn = './incl/bibinc'

c     ouverture du fichier liste des fichiers include
      open( unit=20 , file=newftn, iostat=nerr )
      if( nerr .ne. 0 ) then
         print 10001, newftn
10001 format(' Programme entierw: fichier non ouvrable ',a ) 
         stop
      endif

c     lecture du nom du fichier include
 10   read(20,'(a)', end=8000 , iostat=nerr ) nmftn  
      if( nmftn(1:1) .eq. ' ' ) goto 8000

c     traitement si suffixe= .inc  
      i = index( nmftn , '.bak' )
      if( i .gt. 0 ) goto 10

      i = index( nmftn , '~' )
      if( i .gt. 0 ) goto 10

      i = index( nmftn , '.inc' )
      if( i .le. 0 ) then
         print *, 'entierw: fichier non traite: ',nmftn
         goto 10
      endif
      newftn = nmftn(1:i+4) 
      nmftn  = newftn

c     ouverture de l'ancien fichier fortran
      open( unit=21 , file=nmftn ,iostat=nerr )
      if( nerr .ne. 0 ) then
         print 10001, nmftn
         goto 10
      endif 
      print *,'fichier traite : ',nmftn  

c     lecture du fichier include 
      nligne = 0
 20   nligne = nligne + 1  
      ligne(nligne) = ' '
      read(unit=21,fmt='(a)',iostat=nerr,end=30) ligne(nligne) 
      if( nerr .ne. 0 ) GOTO 9000
      goto 20  

c     fin de lecture du fichier
 30   nligne = nligne - 1  

c     ajout de la declaration INTEGER des variables 
      do nl=1,nligne  
         iligne(nl) = ligne(nl) 
         if( nl .eq. 1 ) then
            iligne(1)(1:16) = '      INTEGER   ' 
         endif
         jv = 1
 35      i = index( iligne(nl)(jv:80), '=' ) 
         if( i .gt. 0 ) then
            j = index( iligne(nl)(jv:80), ',' ) 
            if( j .le. 0 ) then   
               j = index( iligne(nl)(jv:80), ')' ) 
               j = j + 1
            endif 
            j = j - 1
            iligne(nl)(jv-1+i:jv-1+j) = ' ' 
            jv = jv + j + 1
            goto 35 
         endif
      enddo

c     reecriture du fichier include 
      rewind (unit=21)

c     ecriture des lignes 
      do l=1,nligne  
c        recherche du dernier non blanc de la ligne pour 
c        limiter le stockage . suppression des lignes blanches
         do i=80,1,-1
            if( iligne(l)(i:i) .ne. ' ' ) goto 40
         enddo
 40      write(unit=21,fmt='(a)') iligne(l)(1:i)  
      enddo  
    
      do l=1,nligne  
c        recherche du dernier non blanc de la ligne pour 
c        limiter le stockage . suppression des lignes blanches
         do i=80,1,-1
            if( ligne(l)(i:i) .ne. ' ' ) goto 50
         enddo
 50      write(unit=21,fmt='(a)') ligne(l)(1:i)  
      enddo

c     fin du fichier fortran
      close( unit=21 )
      goto 10

c     fin de lecture de la liste des fichiers
 8000 close( unit=20 ) 
      print *,'Fin correcte de l''execution de entierw' 
      stop 

 9000 print*, 'le contenu de la ligne est incorrect ',
     % ligne(nligne)
      goto 20
      end   
