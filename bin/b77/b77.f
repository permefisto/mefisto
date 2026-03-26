      logical function alpha(car)
      character*1 car
c ......................................................................
c dire si un caracter est alphanumerique(majuscule)   ou '_' 
c ......................................................................
c entree:
c     car 
c sortie:
c     nalpha = .true. si oui
c ......................................................................
c
      character*1 carac
      
      carac=car
      icod = ichar(carac)
c 
      if ( (icod.ge.65 .and. icod.le.90) .or.
     .      car .eq. '_' ) then
         alpha = .true.
      else
         alpha = .false.
      endif
      end
      logical function alphanu(car)
      character*1 car
c ......................................................................
c dire si un caracter est alphanumerique(majuscule)  ou 0 ... 9 ou '_' 
c ......................................................................
c entree:
c     car 
c sortie:
c     alphanu = .true. si oui
c ......................................................................
c

      icod = ichar(car)
c 
      if ( (icod.ge.65 .and. icod.le.90) .or.
     .     (icod.ge.48 .and. icod.le.57) .or.
     .      car .eq. '_' ) then
         alphanu = .true.
      else
         alphanu = .false.
      endif
      end
      program b77
c dans un fichier de nom donne remplace tous les noms par un autre noms
C la correspondance entre les noms se trouvera dans le fichier "XXX.rmp" 
c   utilisation:
c    b77  [[-in] input_files *[input_file]] [-out output_file] [-dico dico_file]
c         [-comment] [-indent] [-blancs] 
c                  [-idir prefix *[ prefix] ]   { recherche pour include }
c 
c
%include'std.incl'
      character*1440 ligne,linout(5)
      character*80 buffer,buf5*5,form*4,kl*1 
      integer lenbuf,l,error,llino(5)
      character*80 nom_remp, nomfic
      character*1 ouinon 
      logical alpha,new,open1,open2,open3
      external alpha ,gener_nom
      character*6 gener_nom
      integer maxin
      parameter(maxin=1024)
      integer nbin
      integer adrfile(maxin)

      logical idir,dico,in,out
%INCLUDE '/sys/ins/base.ins.ftn'
%INCLUDE '/sys/ins/pgm.ins.ftn'

        CHARACTER*128 arg_text
        INTEGER*4 argv_ptr, arg_ptr, argv
        INTEGER*2 argc, arg_len, ii
        POINTER /argv_ptr/argv(0:127) {lower bound of dimension}
        POINTER /arg_ptr/arg_len, arg_text



cx   10 format('>> nom du fichier ''.rmp'' : ',$)
cx   30 format('>> continuer ? ',$)
cx   40 format('>> supprimer les commentaires ? ',$)
cx   50 format('>> supprimer les indentations ? ',$)
cx   60 format('>> supprimer les blancs ? ',$)
c
      comment = .false.
      indent  = .false.
      blancs = .false.
      do i=1,10000
         flageti(i)=.false.
      enddo
      flag_implicit=.false.
      buffer=gener_nom(-1,' ') { initialisation de gener_nom }

*    arg_len followed by arg_text

        CALL pgm_$get_args(argc,argv_ptr)

       idir=.false. { pas de regles de recherches pour les includes}
                    { on prend celles par defaut }
       dico=.false.
       in=.true.
       out=.false.
       open1=.false.
       open2=.false.
       open3=.false.
       kfile = 0 
       nbprefix=0
       nbin=0
       DO  ii = 1, argc-1
        arg_ptr = argv(ii)
        if(arg_text(1:arg_len).eq.'-idir')then
          idir=.true.
          in=.false.
        elseif(arg_text(1:arg_len).eq.'-dico')then
          dico=.true.
          idir=.false.
          in=.false.
        elseif(arg_text(1:arg_len).eq.'-in')then
          in=.true.
          idir=.false.
        elseif(arg_text(1:arg_len).eq.'-out')then
          out=.true.
          idir=.false.
          in=.false.
        elseif(arg_text(1:arg_len).eq.'-comment')then
          comment=.true.
          idir=.false.
          in=.false.
        elseif(arg_text(1:arg_len).eq.'-indent')then
          indent=.true.
          idir=.false.
          in=.false.
        elseif(arg_text(1:arg_len).eq.'-blancs')then
          blancs=.true.
          idir=.false.
          in=.false.
        elseif(idir)then { apres l'option -idir, les noms suivants sont les prefix }
          if(nbprefix.lt.maxprefix)then
            nbprefix=nbprefix+1
            prefix(nbprefix)=arg_text(1:arg_len)
          else
            print*,'trop de regles de recherche pour les includes'
          endif
        elseif(dico)then
           dico=.false.   { car un seul dictionnaire }
           open(1,file=arg_text(1:arg_len),status='old',iostat=ioerr)
           if(ioerr.ne.0)then
             print*,'erreur: dico parametre'
             call error_$print(ioerr)
           endif
           open1=ioerr.eq.0
        elseif(in)then
           if(nbin.lt.maxin)then
             nbin=nbin+1
             adrfile(nbin)=ii
           else
             print*,'trop de parametres pour -in'
             stop
           endif
           open2=.true.
        elseif(out)then
           out=.false.
           open(3,file=arg_text(1:arg_len),iostat=ioerr)
           if(ioerr.ne.0)then
             print*,'erreur: out parametre'
             call error_$print(ioerr)
           endif                                   
           open3=ioerr.eq.0
        else    { si pas d'option explicites }
           kfile = kfile + 1
           if(kfile.eq.1) then
            open(2,file=arg_text(1:arg_len),status='old',iostat=ioerr)
            if(ioerr.ne.0)then
              print*,'1er parametre (in)'
              call error_$print(ioerr)
            endif
            open2=ioerr.eq.0
           elseif(kfile.eq.2) then
             open(3,file=arg_text(1:arg_len),iostat=ioerr)
             if(ioerr.ne.0)then
c              print*,'2ieme parametre (out)'
               call error_$print(ioerr)
             endif                                   
             open3=ioerr.eq.0
           endif
        endif
c        WRITE(*,*) 'argument ',ii, ' is ', arg_text(1:arg_len)
       enddo
       if(.not.open3)then
         open(3,file='-stdout')
         open3=.true.
       endif
       if(.not.idir)then
         nbprefix=9
         prefix(1)='include/'
         prefix(2)='/udd/saltel/include/' 
         prefix(3)='/udd/saltel/msh3/include/' 
         prefix(4)='/udd/f3d/include/'
         prefix(5)='/udd/f3d/'
         prefix(6)='//turbot/'
         prefix(7)='/udd/hecht/include/' 
         prefix(8)='ins/'
         prefix(9)='incl/' 
       endif
cx       do i=1,nbprefix
cx         print*,'-idir: ',prefix(i)
cx       enddo
cx       if(indent)then
cx          print*,'garde l''indentation'
cx       else
cx          print*,'suprime l''indentation'
cx       endif
cx       if(blancs)then
cx          print*,'garde les blancs'
cx       else
cx          print*,'suprime les blancs'
cx       endif
cx       if(comment)then
cx          print*,'garde les commentaires'
cx       else
cx          print*,'suprime les commentaires'
cx       endif
cx         do i=1,nbin
cx           arg_ptr = argv(adrfile(i))
cx           print*,i,'-in  "',arg_text(1:arg_len),'"'
cx         enddo
cx

cxc     sdt filein fileout dico
cx      open(2,file='^1',status='old',iostat=ioerr)
cx      if(ioerr.ne.0)then
cx         print*,'1er parametre'
cx         call error_$print(ioerr)
cx      endif
cx      open2=ioerr.eq.0
cx      open(3,file='^2',iostat=ioerr)
cx      if(ioerr.ne.0)then
cxc         print*,'2ieme parametre'
cx         call error_$print(ioerr)
cx      endif                                   
cx      open3=ioerr.eq.0
cx      if(open2.and..not.open3) open(3,file='-stdout')
cx      open(1,file='^3',status='old',iostat=ioerr)
cx      open1=ioerr.eq.0


      nbremp = 0 
      nulig = 0
      if (open1)then
        goto 6 {on lit le dico}
      else
        goto 5 {on saute le dico}
      endif
cx      print*,'========================================================='
cx      print*,'== remplacement d''un nom par un autre dans un fichier  =
cx     +='
cx      print*,'========================================================='
cx      print*,'                                A. Golgolab INRIA 1989'
cx      print*
cx      print* 
cx    1 continue
cx      nbremp = 0 
cx      nulig = 0
cx      print*,'>> nom du fichier  .rmp :'
cx      read(5,'(a)',err=1) nom_remp
cx      l = lenstr(nom_remp)
cx      if ( l .eq. 0 ) goto 5
cx      open(1,status='old',form='formatted',file=nom_remp(1:l)//'.remp'
cx     +     ,err=1) 
cxxx      open(1,status='old',form='formatted',file='noms.remp') 
c
C     --- Lecture du dico ---
    6 continue
      read(1,'(a)',end = 5 ) buffer 
      nulig = nulig + 1
      call correc(buffer)
      l = lenstr(buffer) 
      if ( l .eq. 0 ) goto 6
      ipos = 1
      dowhile ( buffer(ipos:ipos) .ne. ' ' .and. ipos.lt.l )
         ipos = ipos + 1
      enddo
      if ( buffer(ipos:ipos) .eq. ' ' ) then
C        --- premier mot ---
         nbremp = nbremp + 1
         mot_init(nbremp) = buffer(1:ipos-1)
         l_mot1(nbremp) = lenstr(mot_init(nbremp))
      else
C        --- premier mot non suivi du mot a remplacer --- 
         print*,'**** erreur a la ligne ',nulig 
         print*,'--> ',buffer
         print*,'     le remplacant n''est pas precise'
         close(1)
         stop
      endif
c
c     --- lecture du mot remplacant ---  
      ipos = ipos + 1
      dowhile ( buffer(ipos:ipos) .eq. ' ' .and. ipos.lt. l )
         ipos = ipos + 1
      enddo
c
      if ( ipos .eq. l .and. buffer(ipos:ipos).eq.' ' ) then
C        --- premier mot non suivi du mot a remplacer --- 
         print*,'>> erreur a la ligne ',nulig
         print*,'   le remplacant n''est pas precise'
         close(1)
         stop
      endif
c
      if ( (l - ipos + 1) .gt. 6 ) then
         print*,' **** warning a la ligne ',nulig
         print*,'--> ',buffer
         print*,'  longueur du remplacant superieur a 6'
         print*,'  seul les 6 premiers caracters seront consideres'
      endif
      
      new_mot(nbremp) = ' '
      new_mot(nbremp) = buffer(ipos:l)
      l_mot2(nbremp) = lenstr(new_mot(nbremp)) 
c
      if ( new_mot(nbremp) .eq. mot_init(nbremp) ) then
         print*,'**** warning a la ligne :',nulig
         print*,'--> ',buffer
         print*,' mot a remplacer par lui-meme'
         nbremp = nbremp - 1
      endif
      goto 6 
c
c
c     --- fin de la lecture du dictionnaire ---
    5 continue
      close(1)
      if ( nbremp .eq. 0 ) then
         if(.not.open2) print*,' *** dictionnaire vide'
cxxx         goto 1
      endif


cxc =======================================================================
cxc        options
cxc =======================================================================
cx      if(.not.open2) then {on blanchi par defaut}
cx41    print40   {commentaires}
cx      read(5,'(a)') ouinon 
cx      call mimaj(ouinon)
cx      if ( ouinon .eq. 'O'.or.ouinon .eq.' ' ) then
cx         comment = .false.
cx      else if ( ouinon .eq. 'N' ) then
cx         comment = .true.
cx      else
cx         call bip(1)  
cx         print*,'>> repondez par O ou N '
cx         goto 41
cx      endif 
cx
cx 51   print50   {indentations}
cx      read(5,'(a)') ouinon 
cx      call mimaj(ouinon)
cx      if ( ouinon .eq. 'O'.or. ouinon .eq.' ' ) then
cx         indent = .false.
cx      else if ( ouinon .eq. 'N' ) then
cx         indent = .true.
cx      else
cx         call bip(1)  
cx         print*,'>> repondez par O ou N '
cx         goto 51
cx      endif 
cx
cx 61   print60   {blancs}
cx      read(5,'(a)') ouinon 
cx      call mimaj(ouinon)
cx      if ( ouinon .eq. 'O'.or. ouinon .eq.' ' ) then
cx         blancs = .false.
cx      else if ( ouinon .eq. 'N' ) then
cx         blancs = .true.
cx      else
cx         call bip(1)  
cx         print*,'>> repondez par O ou N '
cx         goto 61
cx      endif 
cx      endif




c =======================================================================
c ========= lecture du fichier source ===================================
c =======================================================================
c     
      kfilin=0
2     if(nbin.eq.0) then 
        open(2,status='old',form='formatted',file='-stdin')
      else
        kfilin=kfilin+1
        if(kfilin.gt.nbin)goto 9999
        arg_ptr = argv(adrfile(kfilin))
        open(2,status='old',form='formatted',file=arg_text(1:arg_len))
c        print*,kfilin,' -out: "',arg_text(1:arg_len),'"'
      endif
      nblig=0

    4 continue 
c     --- lecture du fichier  ---
         read(2,'(a)',end=7) buffer
CCC BIDOUILLE POUR SE PLACER DANS LE CAS %include'    ' avec % en colonne 1
         ll = index( buffer,'      include"' )
         if ( ll .gt. 0 ) then
            ll  = ll+14
            lll = lenstr(buffer)-1
            nom_ include = ' '
            nom_include = buffer(ll:lll)
CCC            print *,nom_include
            lll = lll - ll + 2
            nom_include(lll:lll) = 1H'
            buffer = ' '
            buffer = '%include''' // nom_include(1:lll) 
CCC            print *,buffer
         endif 
CCC
         if ( buffer(1:1).eq.'%' ) then
c           --- expendre l'include ---
            call elimine_blancs(buffer)
            nom_ include = ' '
            nom_include = buffer(10:lenstr(buffer)-1)
            call supr_blancs_debut(nom_include) 
            call expand_include 
         else 
            if ( lenstr(buffer) .gt. 0 ) then
               nblig = nblig + 1 
cccc               call mimaj(buffer)
               text(nblig) = buffer 
               if ( buffer(1:1).ne.'C'.and.buffer(1:1).ne.'c'
     +         .and.buffer(1:1).ne.'*' ) then
                  ind = indexx(text(nblig),'{')  
               else
                  ind = 0
               endif
               if ( ind .ne. 0 )then
                  nblig=nblig+1
                  text(nblig)(1:1)='C'
                  text(nblig)(2:)=text(nblig-1)(ind:)
                  text(nblig-1)(ind:) = ' '
               endif
            endif
         endif
      goto 4 
c 
    7 continue
cc     --- met les cartes suites avant les commentaires
      do i=1,nblig
        if(text(i)(6:6).ne.' ') then
          j=i
9967      j=j-1
          if(j.gt.0) then
            if(text(j)(1:1).EQ.'C'.or.text(j)(1:1).eq.'c'
     +      .or.text(j)(1:1).eq.'*') goto 9967
          endif
          j=j+1
          buffer=text(i)
          text(i)=text(j)
          text(j)=buffer
        endif
      enddo
C ======================================================================
c      max des etiquettes utilisees 
c ======================================================================
      indice_max = 0 
      nb_do = 0
      do 11 i = 1 , nblig  
         buf5=text(i)(1:5)
         call supr_blancs_debut(buf5)
         l = lenstr(buf5) 
         if (l.eq. 0 ) goto 11
         write(kl,'(I1)') l
         form = '(i'//kl//')'
         read(buf5(1:l),form,err=11) int
         if(int.gt.0.and.int.lt.10000) flageti(int)=.true.
         indice_max = max(indice_max,int)
   11 continue
c
c =======================================================================
c ============== liste des noms de variable =============================
c =======================================================================  
      nu_max = 1
      call extract_vars(nberr) 
      do i = 1 , nbremp
         call diconu(mot_init(i),num,new,error)
         if ( error .gt. 1 ) nberr = nberr + 1  
         if ( new ) then
            typ_id(num) = external
            new_id(num) = '??????'
         else
            if ( typ_id(num) .eq. external ) then
               new_id(num) = new_mot(i)
            else
               print*,'**** warning : nom utilise est celui d''un s/p:',
     +          mot_init(i)
            endif
         endif
      enddo
         
c =======================================================================
C     standardisation 
C =======================================================================
      flag_implicit=.false.
      buffer=gener_nom(-1,' ') { initialisation de gener_nom }
      ligo = 0 
      nbinclude = 0 
      lig = 1  
      nb_do = 0
      nu_max = 1
c
      dowhile  ( lig .le. nblig )  
c
c        --- ON EST SUR LA LIGNE LIG ---
         if ( text(lig)(1:1).eq.'C'.or.text(lig)(1:1).eq.'c'
     +   .or.text(lig)(1:1).eq.'*')then 
            {c'est un commentaire}
            if (.not. comment ) goto 1000
c            ligo = ligo + 1
c            texto(ligo) = text(lig)
            write(3,'(80(a))')
     +      (text(lig)(i:i),i=1,min(80,lenstr(text(lig))))
            goto 1000
         endif
c
         if ( text(lig)(1:1).eq.'D'.or.text(lig)(1:1).eq.'d') then
c           --- C'EST UNE COMPILATION CONDITIONNELLE ---on suprime la ligne
            goto 1000
         endif 
C
         if ( text(lig)(1:10).eq.'C %include' ) then
c           --- INCLUDE NON TROUVE--- 
            ligo = ligo + 1
            texto(ligo) = text(lig)
            goto 1000
         endif
c        ---- LECTURE DE LA LIGNE COMPLETE ---
         llig = lenstr(text(lig))
         if ( llig .eq. 0 ) goto 1000
         llig = 72
c           --- vrai ligne a scanner ---
         ligne = text(lig)(1:llig) 
c
         dowhile ( lig.lt. nblig .and.
     +             text(lig+1)(1:5) .eq. '     ' .and. 
     +             text(lig+1)(6:6) .ne. ' ' ) 
c           --- lecture des cartes suites ---
            lig = lig + 1
            ligne(llig+1:llig+72-6) = text(lig)(7:)
            llig = llig + 72 - 6 {pour le chaine de plusieur cartes}
         enddo 
c
c 
cccc         print *,ligne(1:llig),';'
c          call elimine_blancs(ligne(7:llig))
c         llig=lenstr(ligne(1:llig))
         CALL EDIT_LINE(LIGNE(1:LLIG),LINOUT,LLINO,NBLIN,ERROR)
         nberr = nberr + error
         if (nblin .eq. 0 ) goto 1000
c
         call format_fortran(linout,llino,nblin)
1000  continue
      lig = lig + 1
      enddo
c
c 
      do i = 1 , ligo 
         l = lenstr(texto(i))
         if ( l .ne. 0 ) then
            write(3,'(a)') texto(i)(1:l)
         endif
      enddo
      close(2)
      if(nberr.ne.0) then
        print*,'>>> ',nberr,' erreur(s) dans la traduction de ',
     +  nomfic(1:lenstr(nomfic))//'.f'
        print*,'>>> fichier genere: ',nomfic(1:lenstr(nomfic))//'.ftn'
      endif
      if(nbin.ge.1)goto 2
9999  continue
      close(3) 
      end










      SUBROUTINE BIP(N)
C SIMULATION EN ATTENDANT QUE BIP.C FONCTIONNE 
      PRINT*,CHAR(7)
      END
      logical function bloc_do(ligne,label)
      character*(*) ligne
      integer label
c ======================================================================
c dire si c'est un bloc do et si c'est un do etiqutte ou non
c ......................................................................
c entree: ligne
c sortie: bloc_do .true. si bloc do
c         label = numero du label ou 0
c ......................................................................
c     A.Golgolab
c ......................................................................
c
      logical eq_ok,virg_ok,numerique
      character*4 form,kl*1
c
      l = len(ligne)
      eq_ok = .false.
      virg_ok = .false.
      bloc_do = .false.
c
      if ( ligne(1:min(l,2)) .ne. 'DO' ) return
      if ( numerique(ligne(3:3)) ) then
         il = 4
         dowhile ( numerique(ligne(il:il)) )
            il = il + 1
         enddo
         il = il - 3
         write(kl,'(i1)') il
         form = '(i'//kl//')'
         read(ligne(3:2+il),form) label
      else
         label = 0
      endif
c
      il = 4
      dowhile ( il.lt.l  .and. (.not.virg_ok .or. .not. eq_ok) ) 
         if ( ligne(il:il) .eq. '(' ) then
            call end_par(ligne(il:), ipos) 
            il = il + ipos
         else if ( ligne(il:il) .eq. '=' ) then
            eq_ok = .true. 
            il = il + 1
         else if ( ligne(il:il) .eq. ',' ) then
            virg_ok = .true.
            il = il + 1 
         else
            il = il + 1
         endif
      enddo
c
      if ( virg_ok .and. eq_ok ) then
         bloc_do = .true.
      endif
      end

            
      SUBROUTINE CORREC( MOT )
C ......................................................................
C BUT : CORRIGER(DANS UNE CERTAINE MESURE) LA SYNTAXE DES MOTS ENTRES
C       AU CLAVIER  ( ABREGE !!! )
C ......................................................................
C ENTREE : MOT  : CHAINE DE CARACTERES CONTENANT LE MOT
C SORTIE : MOT  : LE MOT CORRIGE
C                 CONVERSION EN MAJUSCULE
C ......................................................................
C AUTEURS : A.GOLGOLAB ET X.DENG  ENS-CACHAN  DEC 1987
C ......................................................................
C
      CHARACTER*(*) MOT
      INTEGER LONG,I,POSD,IC,OLDLON
      CHARACTER*1 C
C
      IF ( MOT .EQ. ' ' ) RETURN
      LONG = LEN(MOT)
C
    2 CONTINUE
      IF ( MOT(1:1) .EQ. ' ' ) THEN
         LONG = LONG -1
         MOT(1:LONG) = MOT(2:LONG)//' '
         GOTO 2
      ENDIF
C - CONVERTIR EN MAJUSCULE ET LA POSITION DU DERNIER CARACTERE -
      OLDLON = LONG
      DO 1 I= 1 , LONG
         C = MOT(I:I)
         IF ( C .NE. ' ' ) POSD = I
         IC = ICHAR(C)
         IF ( IC.GE.97 .AND. IC.LE.122 ) MOT(I:I) = CHAR(IC-32)
   1  CONTINUE
      END
      SUBROUTINE DICOFI(NAME,NUNAME,ERROR)
      CHARACTER*(*) NAME
      INTEGER NUNAME,ERROR 
C ......................................................................
C   retourner le numero du nom dans le dictionnaire s'il existe 
C ......................................................................
C entree:
C     NAME : chaine de caractere contenant le nom
C sortie:
C     NUNAME : numero associe au nom ou 0
C     ERROR : code de l'erreur (>0) ou warning (<0) ou 0
C ......................................................................
C     A.Golgolab INRIA mai 1989
C ......................................................................
C
%INCLUDE'dico.incl'
      INTEGER L,NPOSIT,IRANG,I,LENSTR
      EXTERNAL LENSTR
C
      L = LENSTR(NAME)
      IF ( L .EQ. 0 ) THEN
C         CALL BIP(1)
C         WRITE(IMP,*) '> warning DICOFI: Nom vide!' 
         NUNAME = 0
         ERROR = -3
         RETURN
      ENDIF
C
      CALL DICOSK(NAME,L,NPOSIT,IRANG)
C
      IF ( NPOSIT .EQ. 0 ) THEN
C        --- le mot existe dans le dico ---
         NUNAME = CODPNT(IRANG)
      ELSE
C        --- le mot n'existe pas dans le dico ---
         NUNAME = 0
      ENDIF
C
      ERROR = 0 
      END
      SUBROUTINE DICOIN
C ......................................................................
C but : initialisation du dictionnaire
C ......................................................................
C AUTEURS : A.GOLGOLAB & X.DENG  ENS-CACHAN  DEC 1987
C ......................................................................
C
%INCLUDE 'dico.incl'
      CHARACTER*20 NOMFIC,NOMFIR
      INTEGER INFO,INP,LONG
      EXTERNAL INFO 
      COMMON/COMFIC/NOMFIC,NOMFIR
C
C
      NBMOT = 0
      LTABPT = 0
      LONDIC = 1
      CODPNT(1) = 0
      CARACT(1) = '@' 
      CODPNT(2) = 0
      CARACT(2) = ' ' 
      END
      SUBROUTINE DICONA(NUNAME,NAME,ERROR)
      CHARACTER*(*) NAME
      INTEGER NUNAME,ERROR
C ......................................................................
C  retourne le nom a partir de son numero dans le dico  
C ......................................................................
C entree:
C     NUNAME : numero du nom dans le dico
C sortie:
C     NAME : le nom
C     ERROR : code de l'erreur (>0) ou warning (<0) ou 0
C ......................................................................
C     A.Golgolab INRIA mai 1989
C ......................................................................
C
%INCLUDE'dico.incl' 
      INTEGER IRANG,NB
C
      IF ( NUNAME .GT. LTABPT ) THEN
         CALL BIP(1)
         WRITE(IMP,*) '>> Erreur DICONA: numero de nom incorrect: ',
     +   nuname
         ERROR = 6
         RETURN
      ENDIF
C
      IRANG = PNTMOT(NUNAME)
C
      IF ( IRANG .LT. 0 .OR. IRANG .GE. LONDIC ) THEN
         CALL BIP(2)
         WRITE(IMP,*) '>>> ERREUR INTERNE DICONA: IRANG =',IRANG
         ERROR = 7
         RETURN
      ENDIF
C
      NAME = ' '
      IRANG = IRANG + 1
      NB = 0
C
    1 CONTINUE
      IF ( CARACT(IRANG) .NE. '@' ) THEN 
         NB = NB + 1 
         NAME(NB:NB) = CARACT(IRANG)
         IRANG = IRANG + 1
         GOTO 1 
      ENDIF 
      ERROR = 0
      END

      SUBROUTINE DICONU(NAME,NUNAME,NEW,ERROR)
      CHARACTER*(*) NAME
      INTEGER NUNAME,ERROR 
      LOGICAL NEW
C ......................................................................
C   retourner le numero du nom dans le dictionnaire s'il existe 
C   l'ajouter dans le dico s'il n'y est pas et retourner son numero
C ......................................................................
C entree:
C     NAME : chaine de caractere contenant le nom
C sortie:
C     NUNAME : numero associe au nom
C     NEW : true si c'est un nouveau nom
C     ERROR : code de l'erreur (>0) ou warning (<0) ou 0
C ......................................................................
C     A.Golgolab INRIA mai 1989
C ......................................................................
C
%INCLUDE'dico.incl'
      INTEGER L,NPOSIT,IRANG,I,LENSTR
      EXTERNAL LENSTR
C
      L = LENSTR(NAME) 
      error = 0
      IF ( L .EQ. 0 ) THEN
         CALL BIP(1)
         WRITE(IMP,*) '> warning DICOIN: Nom vide a mettre dans le ',
     +   'dico!'
         ERROR = -3
         RETURN
      ENDIF
C
      CALL DICOSK(NAME,L,NPOSIT,IRANG)
C
      IF ( NPOSIT .EQ. 0 ) THEN
C        --- le mot existe deja dans le dico ---
         NEW = .FALSE. 
         NUNAME = CODPNT(IRANG)
         ERROR = 0
         RETURN
      ENDIF
C
C === LE MOT EST NOUVEAU === 
C     ---  y a t il assez de place dans les tableaux ---
      IF ( LTABPT + 1 .GT. MAXMOT ) THEN
         CALL BIP(1)
         WRITE(IMP,*) '>> Erreur DICONU: plus de place dans le tableau'
     +   ,' PNTMOT pour l''ajout d''un nouveau nom dans le dico: ',
     +    NAME(1:L)
         ERROR = 4
         RETURN  
      ELSE IF ( LTABPT + 20 .GT. MAXMOT ) THEN
         CALL BIP(1)
         WRITE(IMP,*) '>> Warning DICONU: le tableau PNTMOT est bientot'
     +   ,' plein. Nombre de places libres: ', MAXMOT-LTABPT-1
         ERROR = -4
      ENDIF 
C
      IF ( LONDIC + L + 1 .GT. MXDICO ) THEN
         CALL BIP(1)
         WRITE(IMP,*) '>> Erreur DICONU: plus de place dans le tableau'
     +   ,' CARACT pour l''ajout d''un nouveau nom dans le dico: ',
     +    NAME(1:L)
         ERROR = 5
         RETURN 
      ELSE IF  ( LONDIC + L +  200 .GT. MXDICO ) THEN 
         CALL BIP(1)
         WRITE(IMP,*) '>> Warning DICONU: le tableau CARACT est bientot'
     +   ,' plein. Nombre de mots libre inferieur a :',MXDICO-LONDIC-L
         ERROR = -5
      ENDIF
C
C --- On peut ajouter le nom --- 
      NEW = .TRUE.
      NBMOT = NBMOT + 1
      LTABPT = LTABPT + 1
      PNTMOT(LTABPT) = LONDIC 
      CODPNT(LONDIC) = NBMOT 
C
      DO I = LONDIC+1 , LONDIC+L
         CODPNT(I) = 0
         CARACT(I) = NAME(I-LONDIC:I-LONDIC)
      ENDDO 
C 
      IF ( NBMOT .GT. 1 ) CODPNT(IRANG) = LONDIC + NPOSIT
      LONDIC = LONDIC + L + 1 
C
      CODPNT(LONDIC) = 0
      CARACT(LONDIC) = '@'
      NUNAME = NBMOT
      END



      SUBROUTINE DICOSK(MOT,LONMOT,NPOSIT,IRANG)
C ......................................................................
C BUT : RECHERCHE DE LA POSITION D'UN MOT DANS LE DICTIONNAIRE ET SI IL
C       N'EXISTE PAS , DONNER LA POSITION DU POINTEUR DE CE MOT DANS LE
C       TABLEAU 
C ......................................................................
C ENTREE : 
C     MOT    : MOT RECHERCHE
C     LONMOT : LONGUEUR DE CE MOT
C SORTIE :
C     NPOSIT : = 0 --> MOT TROUVE
C              > 0 --> POSITION DU CARACTERE DANS LE MOT QUI DOIT 
C                      POINTER SUR 'MOT' COMPTEE A PARTIR DU DEBUT 
C                      DE CE MOT
C
C     IRANG  : POINTEUR SUR LE MOT SI NPOSIT = 0 ( SUR SON CODE )
C              RANG DANS LE SUPER-TABLEAU DU CARACTERE DESIGNE PAR 
C              NPOSIT ( DU MOT QUI DOIT POINTER ) 
C ......................................................................
C AUTEURS : A.GOLGOLAB & X.DENG  ENS-CAHAN  DEC 1987
C ......................................................................
C 
%INCLUDE 'dico.incl' 
C * L = Longueur commune trouvee *
C * P = Pointeur dans le super-tableau * 
      INTEGER LONMOT,NPOSIT,IRANG,L,P
      CHARACTER*(*) MOT 
C * C = Caractere test *
      CHARACTER*1 C,C2*2 
C
      L = 0
      P = 2 
C - Recherche des caracteres communs -
      DOWHILE ( L .LT. LONMOT )
         C = MOT(L+1:L+1)
         DOWHILE( CARACT(P).NE.C .AND. CARACT(P).NE.'@' .AND.
     &            CODPNT(P).NE.0 )
            P = CODPNT(P)
         ENDDO
         IF ( CARACT(P) .EQ. C )  THEN 
C           - Un caractere commun trouve -
            L = L + 1
            P = P + 1 
            GOTO 1
         ENDIF
         IF ( CARACT(P) .EQ. '@' ) THEN
C           - Fin du mot courant , on recule d'un step -
            L = L-1
            P = P-1
         ENDIF
         IF ( CODPNT(P) .EQ. 0 ) THEN
C           - Mot n'existe pas -
            NPOSIT = L+1
            IRANG = P 
            RETURN
         ELSE IF ( CARACT(P+1) .EQ. '@' ) THEN
C           - Fin du mot , le pointeur suivant existe -
            P = CODPNT(P)
         ENDIF 
    1    CONTINUE
      ENDDO
C - Test du dernier caractere -
      C2 = MOT(LONMOT:LONMOT) // '@'
      P = P - 1
      DOWHILE( CARACT(P)//CARACT(P+1) .NE. C2 )
         IF ( CODPNT(P) .NE. 0 ) THEN
            P = CODPNT(P)
         ELSE
C           - Mot n'existe pas -
            NPOSIT = L
            IRANG = P
            RETURN
         ENDIF
      ENDDO
C - Mot trouve -
      NPOSIT = 0
      IRANG = P - LONMOT
      END
      subroutine edit_line(ligne,linout,llino,nblin,error)
      character*(*) ligne,linout(*)
      integer llino(*)
c ======================================================================
c  edition d'une ligne fortran
c ======================================================================
c
%include'std.incl'
      external gener_nom
      character*6 gener_nom
c 
      character*1440 line
      character*35 nom_var(500),nom
      character*5 buf5,instr*8,kl*1,form*4
      integer typ_var(500),adr_var(500),error,typ_do
      logical new 
   11 format(a35,': ',i4,' --> ',a)

c 
      nblin = 0
      error = 0
      nberr = 0  
      erreur = 0
      call elimine_blancs(ligne(7:))
      
c
      if ( nb_do .ne. 0 ) then
c        --- PEUT ETRE LA FIN D'UN DO --- 
         buf5 = ligne(1:5)
         instr = ligne(7:)
c
         call supr_blancs_debut(buf5)
         l = lenstr(buf5) 
         if (l.ne. 0 ) then
c           --- IL EXISTE UNE ETIQUETTE ---
            write(kl,'(I1)') l
            form = '(i'//kl//')'
            read(buf5(1:l),form,err=999) int 
            if ( int .ne. 0 ) then
c              --- IL Y A UNE ETIQUETTE ---  
               if ( pile_do(nb_do) .eq. int ) then
c                 --- C'EST BIEN LA FIN D'UN BLOC DO ETIQUETTE ---
                  nb_do = nb_do - 1
                  if ( instr(1:5) .eq. 'ENDDO' .or.
     +                 instr .eq. 'CONTINUE' ) then 
                     write(linout(1)(1:5),'(I5)') int 
                     linout(1)(6:14) = ' CONTINUE'
                     llino(1) = 14
                     nblin = 1 
                  else 
                     print*,'**** Erreur edit_line: terminer la fin',
     +               ' d''un bloc DO etiquette par un CONTINUE ou un',
     +               ' ENDDO'
                     print*,'ligne: ',ligne 
                     error = 1
                     nblin = 0
                  endif 
                  return
               else if ( instr(1:5) .eq. 'ENDDO' ) then
c                 --- ENDDO MAIS DONT L'ETIQUETTE EST DIFFERENTE ---
                  write(linout(1)(1:5),'(I5)') int 
                  linout(1)(6:14) = ' CONTINUE'
                  llino(1) = 14
c
                  write(linout(2)(1:5),'(I5)') pile_do(nb_do) 
                  linout(2)(6:14) = ' CONTINUE'
                  llino(2) = 14
c
                  nblin = 2 
                  nb_do = nb_do + 1
                  return
               endif
            else
               print*,'**** erreur interne ligne: ',ligne
               nberr = nberr + 1
               nblin = 0
               return
            endif
c
         else if ( instr(1:5) .eq. 'ENDDO' ) then
c           --- IL FAUT DIFFERENCIER DO DE DOWHILE --- 
            if ( nb_do .le. 0 ) then
               print*,'*** Erreur ENDDO sans DO: ',ligne  
               nblin = 0
               return
            endif
c
            if ( pile_do(nb_do) .gt. 0 ) then
c              --- BLOC DO ---
               write(linout(1)(1:5),'(I5)') pile_do(nb_do)
               linout(1)(6:14) = ' CONTINUE'
               llino(1) = 14
               nblin = 1
            else
c              --- FIN D'UN BLOC DOWHILE ---
               linout(1)(1:10) = '      GOTO'
               write(linout(1)(11:15),'(I5)') -pile_do(nb_do) 
               call elimine_blancs(linout(1)(7:15))
               llino(1) = lenstr(linout(1)(1:15))
c
               linout(2)(1:11) = '      ENDIF'
               llino(2) = 11
               nblin = 2
            endif 
            nb_do = nb_do - 1
            return
         endif
      endif
c
c     --- toute la ligne est dans 'ligne' ---
      le = lenstr(ligne(7:))
      linout(1)(1:6) = ligne(1:6)
      linout(1)(7:7) = ' '
      line(1:le) = ligne(7:)
      call scan_vars_lig(line(1:le),nbvars,nom_var,typ_var,
     +                   adr_var,typ_do,error)
      nberr = error + nberr 
c     --- edition de la nouvelle ligne ---
      l = 0
      nuvar = 1
      adr_var(nbvars+1) = 0 
      nblin = 1
      llino(1) = 6 
      if ( typ_do .eq. do ) then
         linout(1)(7:8) = 'DO'
         write(linout(1)(9:13),'(I5)') pile_do(nb_do)
         call elimine_blancs(linout(1)(7:13))
         llino(1) = lenstr(linout(1)(1:13))
         l = 2
      else if ( typ_do .eq. dowhile ) then 
         write(linout(1)(1:5),'(I5)') -pile_do(nb_do)
         linout(1)(7:8) = 'IF'
         llino(1) = 8
         l = 7
      endif
c
      dowhile ( l .lt. le ) 
c
         l = l + 1 
         if ( l .eq. adr_var(nuvar) ) then   
            call diconu(nom_var(nuvar),num,new,error) 
            if (new) then
               print*,'*** erreur interne ligne: ',ligne
               print*,'    identificateur : ',nom_var(nuvar) 
               erreur = erreur + 1
               nblin = 0
               return
            endif
c
            if ( new_id(num)(1:1) .ne. '!' .and. 
     +           new_id(num)(1:1) .ne. '?'      ) then 
c              --- A  REMPLACER ---
               li = lenstr(new_id(num))
               if(blancs)then
                 llino(1)=llino(1)+1
                 linout(1)(llino(1):llino(1)) = ' '
               endif
               do i = 1 , li
                  llino(1) = llino(1) + 1
                  linout(1)(llino(1):llino(1)) = new_id(num)(i:i)
               enddo
               if(blancs)then
                 llino(1)=llino(1)+1
                 linout(1)(llino(1):llino(1)) = ' '
               endif
               l = l + lenstr(nom_var(nuvar)) - 1
               nuvar = nuvar + 1
            else if ( new_id(num)(1:1) .eq. '?' ) then
               print*,'*** erreur interne ligne: ',ligne
               print*,'    on a trouve : ',nom_var(nuvar)
               error = error + 1
               nblin = 0
               return
            else
c              ---  laisser tel quel ----
               li = lenstr(nom_var(nuvar))
               if(blancs)then
                 llino(1)=llino(1)+1
                 linout(1)(llino(1):llino(1)) = ' '
               endif
               do i = 1 , li
                  llino(1) = llino(1) + 1
                  linout(1)(llino(1):llino(1)) = nom_var(nuvar)(i:i)
               enddo
               if(blancs)then
                 llino(1)=llino(1)+1
                 linout(1)(llino(1):llino(1)) = ' '
               endif
               l = l + lenstr(nom_var(nuvar)) - 1
               nuvar = nuvar + 1
            endif
c
         else
            llino(1) = llino(1) + 1
            linout(1)(llino(1):llino(1)) = ligne(6+l:6+l)
         endif
      enddo
c
      if ( typ_do .eq. dowhile ) then
         linout(nblin)(llino(nblin)+1:llino(nblin)+4) = 'THEN'
         llino(nblin) = llino(nblin) + 4
      endif
c     
      if(line(1:4).eq.'END')then
        do i=1,10000
          flageti(i)=.false.
        enddo
        flag_implicit=.false.
        nom=gener_nom(-1,' ') { initialisation de gener_nom }
        indice_max = 0 
        nb_do = 0
        nu_max = 1
        ligo = 0 
      endif
      return


999   continue
       print*,'erreur edit_line'
       print*,'form= "',form,'" buf5= "',buf5(1:l),'"'
      end



      subroutine elimine_blancs(ligne)
      character*(*) ligne
c     ellimine les blancs et converti en majuscule
c
      l = lenstr(ligne)
      ipos = 0
      nb = 0 
      ideb = 1
1     continue
      do i = ideb,l
         if (ligne(i:i) .eq. ' ' ) then
            nb = nb + 1 
         else if ( ligne(i:i) .eq. '''' ) then
            call find_string(ligne(i:),idecal)
            do k = i , i+idecal-1
               ipos = ipos + 1
               ligne(ipos:ipos) = ligne(k:k)
            enddo
            ideb = i + idecal
            goto 1
         else 
            ipos = ipos + 1
            ic=ichar(ligne(i:i))
            if(ic.ge.97.and.ic.le.122)ic=ic-32 { vrais pour l' ASCII }
            ligne(ipos:ipos) = char(ic)
         endif
      enddo
c
      if(nb.gt.0) ligne(l-nb+1:) = ' '
      end
      subroutine end_par(ligne,ipos)
      character*(*) ligne
      integer ipos
c ......................................................................
c position de la parentese fermante
c ......................................................................
c 
      character*1 car
c
c      print*,'end_par "',ligne,'"'
      nbp = 1
      ipos = 1
      dowhile ( nbp .ne. 0. ) 
         ipos = ipos + 1
         if(ipos.gt.len(ligne))then
           ipos=len(ligne)
           return
         endif
         car = ligne(ipos:ipos) 
         if ( car .eq. '''' ) then
            call  find_string(ligne(ipos:),newpos) 
            ipos = ipos + newpos
            car = ligne(ipos:ipos)
         endif
c
         if ( car .eq. '(' ) nbp = nbp + 1
         if ( car .eq. ')' ) nbp = nbp - 1
c
      enddo  
      end
      subroutine expand_include
c ......................................................................
C expansion de l'include RECURSIVEMENT 
C l'include est expandu tel quel
C ......................................................................
%include'std.incl' 
      character*80 buffer
      integer nuprefix
c
      nbinclude = nbinclude + 1
      nunit = 9 + nbinclude
      nuprefix=0
c     
99    continue
      if(nuprefix.eq.0)then
        open(nunit,status='old',form='formatted'
     +     ,file=nom_include,err=9)
      else
        buffer=prefix(nuprefix)(1:lenstr(prefix(nuprefix)))//nom_include
        open(nunit,status='old',form='formatted'
     +     ,file=buffer(1:lenstr(buffer)),err=9)
      endif
      goto 4
9     continue { pas trouve l'include }
      if(nuprefix.gt.nbprefix)goto 1
      nuprefix=nuprefix+1
      goto 99
      
c     on essai avec un prefix
c
    4 continue 
c     --- lecture du fichier  ---
         read(nunit,'(a)',end=2) buffer
         if ( buffer(1:1).eq.'%' ) then
c           --- expendre l'include --- 
            call elimine_blancs(buffer)
            nom_ include = ' '
            nom_include = buffer(10:lenstr(buffer)-1)
            call supr_blancs_debut(nom_include)
            call expand_include 
         else 
            if ( lenstr(buffer) .gt. 0 ) then
               nblig = nblig + 1 
              call mimaj(buffer)
               text(nblig) = buffer
               if (buffer(1:1).ne.'C'.and.buffer(1:1).ne.'c'
     +         .and.buffer(1:1).ne.'*') then
                  ind = indexx(text(nblig),'{')  
               else
                  ind = 0
               endif
               if ( ind .ne. 0 ) text(nblig)(ind:) = ' ' 
            endif
         endif
      goto 4 
c
    2 continue  
      close(nunit)
      nbinclude = nbinclude - 1
      return
c
c
1     continue
      call bip(1)
      print*,'*** erreur include non retrouve : "',nom_include,'"'
      nbinclude = nbinclude - 1 
      ligo = ligo + 1
      texto(ligo) = ' '
      texto(ligo) = 'C %include'''//nom_include(1:lenstr(nom_include))
     +             //''''
      print*,'PREFIXES ESSAYES:'
      do i=1,nbprefix
        print*,'     ',prefix(i)
      enddo
      stop 'erreur: include non retrouve'
      end
      subroutine extract_vars(nberr)
c ======================================================================
c  extraction des noms de variable et mise dans le dico
c ======================================================================
c
%include'std.incl'
      external gener_nom,char_typ
      character*6 gener_nom,char_typ*11
c                          

      character*1440 ligne
      character*35 nom_var(500),nom
      character*5 buf5
      integer typ_var(500),adr_var(500),error,typ_do
      logical new 
   11 format(a35,': ',i4,' --> ',a)

c 
      error = 0
      nulig = 0
      nb_id = 0 
      nbint = 0
      nbrea = 0
      nbaut = 0
      nbcar = 0
      call dicoin
c
      dowhile ( nulig .lt. nblig )
         nulig = nulig + 1
         llig = lenstr(text(nulig))-6 
c
         if ( text(nulig)(1:1).ne.'C' .and. 
     +        text(nulig)(1:1) .ne. 'c' .and.
     +        text(nulig)(1:1) .ne. '*' .and.
     +        text(nulig)(1:1) .ne. 'D' .and.
     +        text(nulig)(1:1) .ne. 'd' .and.
     +        llig .ge. 1                      ) then
c           --- vrai ligne a scanner ---
            ligne(1:llig+1) = text(nulig)(7:) 
c
            dowhile ( nulig.lt. nblig .and.
     +                text(nulig+1)(1:5) .eq. '     ' .and. 
     +                text(nulig+1)(6:6) .ne. ' ' ) 
c              --- lecture des cartes suites ---
               nulig = nulig + 1
               l = lenstr(text(nulig))
               if ( l .ge. 7 ) then
                  ligne(llig+1:llig+l) = text(nulig)(7:)
                  llig = llig + l - 6
               endif
            enddo 
         else
            goto 1000
         endif
c
c        --- toute la ligne est dans 'ligne' ---
         call elimine_blancs(ligne(1:llig)) 
d         print*
d         print*,'12345678901234567890123456789012345678901234567890'//
d    +   '12345678901234567890123456789012345678901234567890'
d         print*,ligne(1:llig)
         call scan_vars_lig(ligne(1:llig),nbvars,nom_var,typ_var,adr_var
     +                      ,typ_do,error)
         nberr = error + nberr

d         do i = 1 , nbvars
d            print*,'--> ', typ_var(i),' -- ',adr_var(i),' -> ',
d    +                     nom_var(i)
d         enddo
c        --- mise dans le dico ---
         do i = 1 , nbvars 
            call diconu(nom_var(i),num,new,error)
            nuv = 0
            nb_id = max(nb_id,num)
            if ( error .gt. 0 ) then
               print*,'erreur dans le dico. mot: ',nom_var(i) 
               nberr = nberr + 1
            endif 
            ity = typ_var(i)
c
            if ( ity .eq. external ) then
               new_id(num) = '!!!!!!' 
               typ_id(num) = external
               nuv = 0  
c
            else if (ity .eq. integer .or. ity.eq.ivar) then
               if (new) then
                  new_id(num) = 'I00000' 
                  nbint = nbint + 1 
                  nuv = nbint 
                  typ_id(num) = ity
               endif
c
            else if ( ity .eq. real .or. ity.eq.rvar) then
               if (new) then
                  new_id(num) = 'R00000'
                  nbrea = nbrea + 1 
                  nuv = nbrea 
                  typ_id(num) = ity
               endif
c
            else if ( ity .eq. character ) then
               if (new) then
                  new_id(num) = 'K00000'
                  nbcar = nbcar + 1 
                  nuv = nbcar 
                  typ_id(num) = ity
               else
                  typ_id(num) = ity
               endif
c
            else if ( ity .eq. logical .or.ity .eq. double  ) then
               if (new) then
                  new_id(num) = 'A00000' 
                  nbaut = nbaut + 1
                  nuv = nbaut
                  typ_id(num) = ity
               endif 
c
            else if (ity .eq. integer_tab ) then
               if ( new ) then
                  new_id(num) = 'I00000'
                  nbint = nbint + 1
                  nuv = nbint
                  typ_id(num) = integer_tab
               else
                  typ_id(num) = integer_tab
               endif
c
            else if ( ity .eq. real_tab ) then
               if ( new ) then
                  new_id(num) = 'R00000'
                  nbrea = nbrea + 1
                  nuv = nbrea
                  typ_id(num) = real_tab
               else
                  typ_id(num) = real_tab
               endif
c
            else if ( ity .eq. logical_tab ) then
               if ( new ) then
                  new_id(num) = 'A00000'
                  nbaut = nbaut + 1
                  nuv = nbaut
                  typ_id(num) = logical_tab
               else
                  typ_id(num) = logical_tab
               endif
c
            else if ( ity .eq. double_tab ) then
               if ( new ) then
                  new_id(num) = 'A00000'
                  nbaut = nbaut + 1
                  nuv = nbaut
                  typ_id(num) = double_tab
               else
                  typ_id(num) = double_tab
               endif
c  
            else if ( ity .eq. char_tab ) then
               if ( new ) then
                  new_id(num) = 'A00000'
                  nbaut = nbaut + 1
                  nuv = nbaut
                  typ_id(num) = char_tab
               else
                  typ_id(num) = char_tab
               endif
c
c
            else if ( ity .eq. itab ) then
               if ( new ) then
                  new_id(num) = 'I00000'
                  nbint = nbint + 1
                  nuv = nbint
                  typ_id(num) = itab
               else 
                  typ_id(num) = typ_id(num) + 100
               endif
c
            else if ( ity .eq. rtab ) then
               if ( new ) then
                  new_id(num) = 'R00000'
                  nbrea = nbrea + 1
                  nuv = nbrea
                  typ_id(num) = rtab
               else 
                  typ_id(num) = typ_id(num) + 100
               endif
c
c 
c
c
            else if ( typ_var(i) .eq. ext_tab ) then
               nuv = 0
               if ( new ) then
                  new_id(num) = '!!!!!!'
                  typ_id(num) = external
               else
                  ityp = typ_var(i) 
                  idt = typ_id(num)
                  if ( idt.eq.integer .or. idt.eq.real .or.
     +                 idt.eq.logical .or. idt.eq.character .or.
     +                 idt.eq.double  .or. idt.eq.ivar .or.
     +                 idt.eq.rvar ) then
c                    --- c'est une fonction ---
                     new_id(num) = '!!!!!!'
                     typ_id(num) = external
                  endif
               endif
c
c
c
            else
             call bip(1)
             print*,'**** extract_vars:Erreur interne: ',ligne(1:llig) 
             nuv = 0
            endif   
c 
            if ( nuv .ne. 0 ) then
                new_id(num)(1:6)=gener_nom(typ_id(num),nom_var(i)(1:1))
ccxx               write(buf5,'(i5.5)') nuv
ccxx               call supr_blancs_debut(buf5)
ccxx               l = lenstr(buf5)
ccxx               new_id(num)(6-l+1:6) = buf5(1:l) 
c               print*,'ancien nom:"',nom_var(i)(1:12)
c     +               ,'" -->:"',new_id(num)(1:6),'"'
c     +               ,' type="',char_typ(typ_id(num)),'"'
            endif
         enddo 
1000  continue
      enddo
c
c
cx      do i = 1 , nb_id  
cx         call dicona(i,nom,error) 
cx         print11,nom,typ_id(i),new_id(i)
cx      enddo 
c
      end




      subroutine find_string(string,ipos) 
      character*(*) string
      integer ipos
c ........................................................
c     trouve la fin d'une chaine de caracteres
c ........................................................
c entree:
c     string commencant par un '
c sortie:
c     ipos: position de la fin  de la chaine de 
c ........................................................
c A.Golgolab
c ........................................................
c 
      character*1 car
c
      ll = len(string)
      nbc = 1
      ipos =  1
      dowhile(nbc.ne.0 .and. ipos.lt.ll) 
        ipos = ipos + 1
        if ( string(ipos:ipos) .eq. '''' ) then
            nbc = mod(nbc+1,2) 
            if ( ipos .lt. ll ) then
               if ( string(ipos+1:ipos+1) .eq. '''' ) then
                  nbc = mod(nbc+1,2)
                  ipos = ipos +1
               endif 
            endif
         endif
      enddo 
      end

      subroutine format_fortran(line,llin,nblin)
c
      character*(*) line(*)
      integer llin(*)
%include'std.incl'
      integer count
      save count
      data count=0
c 
      do i = 1 , nblin
         count=count+1
         l = llin(i)
         ii = max(72 - l,0)
         if(indent)then
           ii=0
         else
           ii = nint( (sin(float(ii*count))**2)*ii)
         endif
         write(3,'(72(a))')  line(i)(1:min(l,6))
     +                      ,(' ',k=1,ii)
     +                       ,(line(i)(k:k),k=7,min(l,72))
         l = l -72 
         lp = 72
         dowhile(l.gt.0)
            write(3,'(a,a)') '     .',line(i)(lp+1:lp+min(l,66))
            l = l - 66
            lp = lp + 66
         enddo 
      enddo
      end
      character*6 function gener_nom(typ,char1)
      integer typ
      character*1 char1
c ======================================================================
c     gener un nom de 6 character en fonction de typ :
c          char1 est le 1er caracter de l'ancien nom
c ======================================================================
%include'std.incl'

      integer base,nuv
      parameter(base=5)
      character*1 digit(0:base-1),ctyp(9)
      save nuv
      data digit/'O','0','1','I','Q'/
      data ctyp/'I','O','O','O','O','O','O','I','O'/

      if(typ.lt.0)then
        nuv=-1 { initialisation }
        return
      endif
      nuv = nuv + 1
      if(nuv.ge.5**base)then
        print*,'GENER_NOM: ERREUR TROP DE VARIABLES',nuv
        print*,'LE MAXIMUM EST:',5**BASE
        print*,'MODIFIER LA BASE DE NUMEROTATION!'
        stop 'BASE TROP PETITE'
      endif
      n=nuv
      if(flag_implicit)then
        gener_nom(1:1)=char1
      else                             
        if(typ.lt.100)then
          gener_nom(1:1)=ctyp(typ)
        else
          gener_nom(1:1)=ctyp(typ-100)
        endif
      endif

      do i=1,5
         k=mod(n,base)
         gener_nom(i+1:i+1)=digit(k)
         n=n/base
      enddo
c      print*,' gener_nom="',gener_nom,'"'
c     +      ,' nuv=',nuv,' typ=',typ
      end

      character*11 function char_typ(typ)
      integer typ
c ======================================================================
c     converti le typ en caracter
c ======================================================================
                                 
      character*11 ctyp1(9)
      character*11 ctyp2(9)
      data ctyp1/'integer'
     +          ,'real'
     +          ,'logical'
     +          ,'double'
     +          ,'character'
     +          ,' '
     +          ,' '
     +          ,'ivar'
     +          ,'rvar'/
      data ctyp2/'integer_tab'
     +          ,'real_tab'
     +          ,'logical_tab'
     +          ,'double_tab'
     +          ,'char_tab'
     +          ,'external'
     +          ,'ext_tab'
     +          ,'itab'
     +          ,'rtab'/

      if(typ.lt.100)then
        char_typ=ctyp1(typ)
      else
        char_typ=ctyp2(typ-100)
      endif
      end
      integer function indexx(lig,txt)
      character*(*) lig,txt
c ......................................................................
c index mais pas dans les chaines de caracteres
c ......................................................................
c 
      indexx = index(lig,txt)
      if ( indexx .eq. 0 ) return  {pas d'ambiguite}
c
      ind = index(lig,'''')
      if ( ind .eq. 0 ) return   {pas d'ambiguite}
c
c     --- il existe une chaine de carcteres ---
      il = 1
      llig = len(lig) 
      ltxt = len(txt)
      indexx = 0
c
      dowhile ( il .le. llig-ltxt )
         if ( lig(il:il) .eq. '''' ) then
            call find_string(lig(il:),ipos)
            il = il + ipos
         else
            if ( lig(il:il+ltxt-1) .eq. txt ) then
               indexx = il
               return  
            else
               il = il + 1
            endif
         endif
      enddo
      end




      subroutine isa_mot_cle(mot,oui)
      character*(*) mot
      logical oui
c ......................................................................
c     dit si mot est un mot cle du data
c ......................................................................
      integer maxn
      parameter (maxn=22)

      character*12 liste(maxn),nom
      data liste/ 'UNIT=       ', 'IOSTAT=     ', 'ERR=        ',
     +            'FILE=       ', 'STATUS=     ', 'ACCESS=     ',
     +            'FORM=       ', 'RECL=       ', 'BLANK=      ',
     +            'NUMBER=     ', 'NAMED=      ', 'NAME=       ',
     +            'SEQUENTIAL= ', 'DIRECT=     ', 'UNFORMATTED=',
     +            'FORMATTED=  ', 'NEXTREC=    ', 'END=        ',
     +            'EXIST=      ', 'OPENED=     ',
     +            'FMT=        ', 'REC=        '
     +        / 
C  
C         

      do i = 1 , maxn
         if ( mot .eq. liste(i)(1:len(mot)) ) then
            oui = .true.
            return
         endif
      enddo
c
      oui = .false.
      end
      INTEGER FUNCTION LENSTR(STR)
C ======================================================================
C BUT : RENVOYER LA LONGUEUR D'UNE CHAINE DE CARACTERES SANS COMPTER
C       LES BLANCS DE LA FIN
C ======================================================================
C ENTREE : STR : LA CHAINE DE CARACTERES
C SORTIE :
C     LONSTR : LONGUEUR SANS COMPTER LES BLANCS
C ======================================================================
C PROGRAMMEUR : ARDESHIR GOLGOLAB  INRIA AVRIL 1988
C ======================================================================
C
      CHARACTER*(*) STR
C
      LONG = LEN(STR)
C
      DO 1 I = LONG , 1 , -1
         IF ( STR(I:I) .NE. ' ' ) THEN
            LENSTR = I
            RETURN
         ENDIF
    1 CONTINUE
      LENSTR = 0
      END
      subroutine mimaj( mot )
c ......................................................................
c but : corriger(dans une certaine mesure) la syntaxe des mots entres
c       au clavier  ( abrege !!! )
c ......................................................................
c entree : mot  : chaine de caracteres contenant le mot
c sortie : mot  : le mot corrige
c                 conversion en majuscule
c ......................................................................
c auteurs : a.golgolab et x.deng  ens-cachan  dec 1987
c ......................................................................
c
      character*(*) mot
      integer long,i,ic
c
      if ( mot .eq. ' ' ) return
      long = len(mot)
c
c - convertir en majuscule et la position du dernier caractere -
      do 1 i= 1 , long
         ic = ichar(mot(i:i))
         if ( ic.ge.97 .and. ic.le.122 ) mot(i:i) = char(ic-32)
   1  continue
      end
      logical function numerique(car)
      implicit integer (i-n)
      character*1 car
c ......................................................................
c dire si un caracter est numerique  0 ... 9 
c ......................................................................
c entree:
c     car 
c sortie:
c     nalpha = .true. si oui
c ......................................................................
c
      icod = ichar(car)
c 
      if ( icod.ge.48 .and. icod.le.57 ) then
         numerique = .true.
      else
         numerique = .false.
      endif
      end
      subroutine remplace(ligne,nbtrsf)
      character*(*) ligne
c ======================================================================
c remplace dans une ligne une liste de mot par une autre
c ====================================================================== 
c entree: ligne
c sortie : nbtrsf = nombre de transformations effectuees
c ======================================================================
c
%include'std.incl'
      logical alphanu
      external alphanu 
      character*80 buffer
c
      do num = 1 , nbremp
c        --- num est le numero du mot a tester --- 
         ind = index( ligne , mot_init(num)(1:l_mot1(num)) )
c
         dowhile ( ind .ne. 0 .and. 
     +            .not. alphanu(ligne(ind+l_mot1(num):
     +                              ind+l_mot1(num)) ) )
            ligne(ind:ind+l_mot1(num)-1) = ' '  
            if ( l_mot2(num) .gt. l_mot1(num) ) then 
c              --- le remplacant est plus long ---
               buffer = ' '
               buffer = ligne(ind+l_mot1(num):)
               ligne(ind:ind+l_mot2(num)-1) = new_mot(num) 
               ligne(ind+l_mot2(num):)  = buffer
               ind = index(ligne,mot_init(num)(1:l_mot1(num))) 
               nbtrsf = nbtrsf + 1
            else
               ligne(ind:ind+l_mot2(num)-1) = new_mot(num) 
               buffer = ' '
               buffer = ligne(ind+l_mot1(num):)
               ligne(ind+l_mot2(num):)  = ' '
               ligne(ind+l_mot2(num):)  = buffer
               ind = index(ligne,mot_init(num)(1:l_mot1(num))) 
               nbtrsf = nbtrsf + 1 
            endif
         enddo 
      enddo 
      end
      SUBROUTINE scan_vars_lig(ligne,nbvars,nom_var,typ_var,adr_var,
     +                         typ_do,error)
      integer nbvars,llig,error,typ_var(*),adr_var(*),typ_do
      character*(*) ligne,nom_var(*)
c ......................................................................
c scanner une ligne fortran et sortir toutes les variables
c ......................................................................
c entree:
c     ligne : ligne fortran (cartes suites supprimes)
c sortie:
c     nbvars: nombres de variables dans la ligne (identificateurs)
c     nom_var : nom de ces variables
c     typ_var : type des ces variables
c ......................................................................
c A.Golgolab INRIA juillet 1989
c ......................................................................
c
%include'std.incl'
c
      character*50 buffer,buf*1440,sep*1,tok*35,car*1,ch4*4,ch5*5
      logical alpha,numerique,alphanu,dctab,oui,bloc_do
      external alpha,numerique,alphanu
c 
      llig = lenstr(ligne)
      nbvars = 0
      error = 0 
      idecal = 0  
      typ_do = pas_un_do
c =====================================================================
c cas speciaux 
c ===================================================================== 
c      print*,'scan_vars_lig: ',ligne(1:lenstr(ligne))
      buf = ligne(1:llig)  
      l = lenstr(buf)
      lbuf=l
      call mimaj(buf(1:lbuf))
      if ( l.eq. 0 ) return
      call supr_blancs_debut(buf)
      ind = indexx(ligne(1:llig),'=')
         if (
c     +        buf(1:4) .eq. 'SAVE' .or.
c     +        buf(1:4) .eq. 'ELSE' .or.
     +       (buf(1:3) .eq. 'END'.and.l.eq.3)  .or.
     +       (buf(1:5) .eq. 'ENDIF'.and.l.eq.5)  .or.
     +       (buf(1:5) .eq. 'ENDDO'.and.l.eq.5)  .or.
     +       (buf(1:8) .eq. 'IMPLICIT'.and.ind.eq.0)  .or. 
c     +        buf(1:5).eq. 'STOP'''  .or.
c     +        buf(1:4).eq. 'STOP'  .or.
c     +        buf(1:6).eq. 'PAUSE'''  .or.
c     +        buf(1:5) .eq. 'PAUSE'  .or.
     +       (buf(1:8) .eq. 'CONTINUE'.and.l.eq.8) ) then
c            if(buf(1:8).eq.'IMPLICIT')then
c              flag_implicit=.true.  
c            endif
            return
         endif  
      buf = ligne(1:llig)  
c
c ======================================================================
c type d'istruction
c ======================================================================
c
1000  continue
      ind = indexx(ligne(1:llig),'=')
c
         if (buf(1:min(llig,15)) .eq. 'INTEGERFUNCTION'  ) then
           ibloc = 11 
           mode = call
           ld = 15
         else if(buf(1:min(llig,15)) .eq. 'LOGICALFUNCTION') then
           ibloc = 11 
           mode = call
           ld = 15
         else if(buf(1:min(llig,12)) .eq. 'REALFUNCTION') then
           ld = 12
           ibloc = 11 
           mode = call
         else if(buf(1:min(llig,14)) .eq. 'REAL*8FUNCTION')then
           ld = 14
           ibloc = 11 
           mode = call
         else if(buf(1:min(llig,23)).eq.'DOUBLEPRECISIONFUNCTION')then
           ld = 23
           ibloc = 11 
           mode = call
         else if(buf(1:min(llig,10)) .eq. 'ENTRY'.and.ind.eq.0 ) then
           ld=5
           ibloc = 11 
           mode = call
         elseif(buf(1:min(llig,10)).eq.'SUBROUTINE'.and.ind.eq.0)then
            ld = 10
           ibloc = 11 
           mode = call
         elseif(buf(1:min(llig,8)).eq.'FUNCTION'.and.ind.eq.0)then
            ld = 10
           ibloc = 11 
           mode = call
c
      else if ( buf(1:min(llig,9)) .eq. 'CHARACTER' .and. 
     +          indexx(buf(1:llig),'FUNCTION') .ne. 0 .and.
     +          ind .eq. 0 ) then
         ibloc = 11 
         mode = call
         ld = indexx(buf(1:llig),'FUNCTION') + 7
      else if ( buf(1:min(llig,7)) .eq. 'INTEGER' .and. ind.eq.0) then
         mode = integer 
         ibloc = 1 
         ld = 7

      else if ( buf(1:min(llig,7)).eq.'LOGICAL' .and. ind.eq.0 ) then
         mode = logical
         ibloc = 1
         ld = 7
      else if ( buf(1:min(llig,7)).eq.'COMPLEX' .and. ind.eq.0 ) then 
         mode = real
         ibloc = 1
         ld = 7
      else if ( buf(1:min(llig,6)) .eq.'REAL*8' .and. ind.eq.0) then
         mode = double
         ibloc = 1  
         ld = 6
      else if ( buf(1:min(llig,15)) .eq. 'DOUBLEPRECISION' 
     +           .and. ind.eq.0) then 
         mode = double
         ibloc = 1
         ld = 15
      else if ( buf(1:min(llig,4)) .eq.  'REAL' .and. ind.eq.0) then
         mode = real
         ibloc = 1 
         ld = 4
      else if ( buf(1:min(llig,9)).eq.'CHARACTER'.and. ind.eq.0) then
         mode = character
         ibloc = 1
         ld = 9
         dowhile ( .not. alpha(buf(ld+1:ld+1)) )
            ld = ld + 1
         enddo
      else if ( buf(1:min(llig,9)).eq.'DIMENSION' .and. ind.eq.0) then
         ibloc = 1
         mode = declaration 
         ld = 9
      else if ( buf(1:min(llig,9)).eq.  'PARAMETER' ) then
         mode = declaration
         ibloc = 7  
         ld = 9
      else if (buf(1:min(llig,11)).eq.'EQUIVALENCE'.and.ind.eq.0)then
         mode = declaration
         ibloc = 7  
         ld = 11
      else if ( buf(1:min(llig,6)) .eq.  'COMMON' .and. ind.eq.0) then
         mode = declaration
         ibloc = 8  
         ld = 6
      else if ( buf(1:min(llig,7)) .eq.  'POINTER' .and. ind.eq.0) then
         mode = declaration
         ibloc = 2  
         ld = 7
      else if ( buf(1:min(llig,9)).eq.'BLOCKDATA'.and.ind.eq.0)then
         mode = declaration
         ibloc = 8  
         ld = 9
      else if ( buf(1:min(llig,8)) .eq.'EXTERNAL' .and. ind.eq.0) then
         mode = external
         ibloc = 1
         ld = 8
      else if ( buf(1:min(llig,9)).eq.'INTRINSIC' .and. ind.eq.0)then
         mode = external
         ibloc = 1
         ld = 9
      else if ( buf(1:min(llig,4)) .eq.  'DATA'   .and. ind.eq.0) then
         mode = declaration
         ibloc = 2  
         ld = 4
      else if ( buf(1:min(llig,6)) .eq.  'FORMAT' ) then
         mode = format 
         ibloc = 0
      else if ( buf(1:min(llig,3)) .eq.  'IF(' ) then
         mode = if
         ibloc = 3 
         ld = 3
      else if ( buf(1:min(llig,7)) .eq.'ELSEIF(' .and. ind.eq.0) then
         mode = elseif
         ibloc = 3 
         ld = 7
      else if ( buf(1:min(llig,4)) .eq.'ELSE' .and. ind.eq.0) then
         mode = elseif
         ibloc = 3 
         ld = 7
         return
      else if ( buf(1:min(llig,4)) .eq.  'CALL' .and. ind.eq.0) then
         mode = call
         ld = 4 
         ibloc = 3
      else if ( buf(1:min(llig,8)) .eq.'DOWHILE(' .and. ind.eq.0) then
         mode = dowhile
         ibloc = 3
         ld = 8
         nb_do = nb_do + 1
         nu_max = nu_max + 1
         dowhile (flageti(nu_max))
           nu_max = nu_max + 1
         enddo
         pile_do(nb_do) = - nu_max 
         typ_do = dowhile
c      else if ( buf(1:min(llig,4)) .eq.  'GOTO' .and. ind.eq.0) then
c         mode = goto
c         ld=4
c         ibloc = 0
      else if ( buf(1:min(llig,4)) .eq.  'SAVE' .and. ind.eq.0) then
         ibloc = 3 
         mode = save
         ld = 4
      else if ( buf(1:min(llig,5)) .eq.  'PAUSE' .and. ind.eq.0) then
         ibloc = 3 
         mode = save
         ld = 5
      else if ( buf(1:min(llig,4)) .eq.  'STOP' .and. ind.eq.0) then
         ibloc = 3 
         mode = save
         ld = 4
      else if ( buf(1:min(llig,4)).eq.'IMPLICIT')then
         ld=llig
         flag_implicit=.true.  
         return
      else if ( buf(1:min(llig,5)) .eq. 'PRINT' .or. 
     +          buf(1:min(llig,5)) .eq. 'WRITE'  .or. 
     +          (buf(1:min(llig,6)) .eq. 'RETURN'.and.ind.eq.0)  .or.
     +          (buf(1:min(llig,6)) .eq. 'ASSIGN'.and.ind.eq.0) .or.
     +          (buf(1:min(llig,4)) .eq. 'GOTO'.and.ind.eq.0)  .or.
     +          buf(1:min(llig,4)) .eq. 'READ'      ) then 
         if ( buf(1:min(llig,5)) .eq. 'PRINT' .or.
     +        buf(1:min(llig,5)) .eq. 'WRITE'      ) then
            ld = 5
         elseif(buf(1:min(llig,4)).eq.'READ')then
            ld = 4
         elseif(buf(1:min(llig,6)).eq.'RETURN')then
            ld=6
         elseif(buf(1:min(llig,4)).eq.'GOTO')then
            ld=4
         elseif(buf(1:min(llig,6)).eq.'ASSIGN')then
            ld=6
         endif
c        
         if ( ld .ge. llig ) goto 88888
         if ( buf(ld+1:ld+1) .eq. '*' ) then
            ld = ld + 2  
            ibloc = 3
         else 
            if ( buf(1:1) .eq. 'P' ) then
             if(buf(ld+1:ld+1).eq.'''')then
               ibloc = 3
             else  if(buf(ld+1:ld+1).eq.'*')then
               ld = ld + 1
               if ( ld .gt. llig ) goto 88888
               ibloc = 3
             else  if(.not.numerique(buf(ld+1:ld+1)))then
               ibloc = 3
             else  
               dowhile ( numerique(buf(ld:ld) ) )
                  ld = ld + 1 
                  if ( ld .gt. llig ) goto 88888
               enddo 
               ibloc = 3
             endif
            elseif(buf(1:min(llig,6)).eq.'RETURN')then
              ibloc=3
            elseif(buf(1:min(llig,6)).eq.'ASSIGN')then
              do while(buf(ld:ld) .ne.'T')  { on saute le label du ASSIGN }
                ld=ld+1
                if ( ld .gt. llig ) goto 88888
              enddo 
              ld=ld+1
              ibloc=3
            elseif(buf(1:min(llig,4)).eq.'GOTO')then
               if ( buf(ld+1:ld+1) .eq. '(' ) then
                  ibloc = 13  
                  new_pos=llig-ld
                  llig_save =0
               else
                do while(numerique(buf(ld:ld)))
                  ld=ld+1
                  if ( ld .gt. llig ) goto 88888
                enddo 
                ibloc=3
               endif
            else   {READ ou WRITE}
               if ( buf(ld+1:ld+1) .eq. '(' ) then
                  ibloc = 13  
                  call end_par(buf(ld+1:), new_pos )
                  llig_save = llig - new_pos -ld  
                  llig = ld + new_pos
               else if ( buf(ld+1:ld+1) .eq. '*' ) then
                  ld = ld + 2
                  ibloc = 3
               else 
c                  print*,'>>> erreur : ',buf(1:llig)
c                  error = 1
c                  goto 88888
                 ibloc=3
               endif 
            endif
         endif
         mode = save 
c
      else if ( buf(1:min(llig,6)) .eq. 'CLOSE(')then
         ibloc = 12
         mode =  i_o
         ld=6
      elseif  (buf(1:min(llig,7))  .eq. 'REWIND(')then
         ibloc = 12
         mode =  i_o
         ld=6
      elseif  (buf(1:min(llig,6))  .eq. 'REWIND'
     +          .and.buf(min(llig,7):min(llig,7)).ne.'=')then
         ibloc = 6
         mode =  aff
         ld=6
      elseif  (buf(1:min(llig,10))  .eq. 'BACKSPACE(')then
         ibloc = 12
         mode =  i_o
         ld=10
      elseif  (buf(1:min(llig,9))  .eq. 'BACKSPACE' 
     +          .and.buf(min(llig,10):min(llig,10)).ne.'=')then
         ibloc = 6
         mode =  aff
         ld=9
      elseif  (buf(1:min(llig,8))  .eq. 'ENDFILE(')then
         ibloc = 12
         mode =  i_o
         ld=8
      elseif  (buf(1:min(llig,7))  .eq. 'ENDFILE'
     +          .and.buf(min(llig,8):min(llig,8)).ne.'=')then
         ibloc = 6
         mode =  aff
         ld=7
      elseif  (buf(1:min(llig,5))  .eq. 'OPEN(')then
         ibloc = 12
         mode =  i_o
         ld=5
      elseif  (buf(1:min(llig,8))  .eq. 'INQUIRE('  ) then
         ibloc = 12
         mode =  i_o
         ld=8
      elseif ( buf(1:min(llig,6)) .eq. 'RETURN' ) then
            ld = 6
            ibloc = 6
            mode = aff 
c
      else if ( buf(1:min(llig,7)) .eq. 'PROGRAM' .and. ind.eq.0 .and.
     +          index(buf(min(llig,7):min(llig,7)),'(') .eq. 0 ) then
         ibloc = 0 
         mode = 0
c
      else if ( bloc_do(buf(1:llig),label) ) then
         mode = save
         ld = 2
         ibloc = 3 
         if ( label .eq. 0 ) then
            typ_do = do
            nb_do = nb_do + 1
            nu_max = nu_max + 1
            dowhile (flageti(nu_max))
              nu_max = nu_max + 1
            enddo
            pile_do(nb_do) = nu_max
         else
c            nb_do = nb_do + 1
c            pile_do(nb_do) = label
         endif
c
      else if ( ind .ne. 0 ) then
         mode = aff 
         ld = 0
         ibloc = 6 
c 
      else if ( buf(1:min(llig,6)).eq. 'RETURN' ) then
         ibloc = 0
      else
c         ibloc = 0 
         ibloc = 3
c         error = 1
         ipos=1
         llig=lenstr(buf)
c         print*,'>>>>> type inconnu:"',buf(1:llig),'"'
115      continue
c         print*,'"',buf(ipos:llig),'"'
         if(ipos.gt.llig)goto 88888
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
         ipos = ipos+idpos -1
         if ( ltok .gt. 0 ) then
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal +ipos - ltok
            call set_mode(tok(1:ltok),
     +                       buf(min(ipos,llig):llig), mode2)
            typ_var(nbvars) = mode
            goto 115 
         else
            goto 88888
         endif 
      endif  
c
c
1001  continue
c ======================================================================
c    traitement
c ======================================================================
c
      if ( ibloc .eq. 0 ) then 
c        --------------------------------------------------------------
c        RIEN
c        --------------------------------------------------------------
         goto 88888
c
      else if ( ibloc .eq. 1 ) then
c        --------------------------------------------------------------
c        --- DECLARATIONS --- 
c        --------------------------------------------------------------
         ipos = ld+1 
    1    continue 
         if ( ipos .gt. llig ) goto 88888
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
c
         if ( ltok .gt. 0 ) then
            ipos = ipos+idpos   {***??? -1 ???***}
            if (buf(min(llig,ipos-1):min(llig,ipos-1)).eq.'(' ) then
c              --- c'est une declaration de tableau ---
               dctab = .true.
            else
               dctab = .false.
            endif
c
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal + ipos - ltok-1
c  
            if ( mode .eq. declaration ) then
c              --- type implicite ---  
               if ( ichar(tok(1:1)) .ge. 73 .and.
     +              ichar(tok(1:1)) .le. 78    ) then
                  mode2 = ivar
               else
                  mode2 = rvar
               endif 
            else
               mode2 = mode
            endif 
c
            if ( dctab ) then
               typ_var(nbvars) = mode2 + 100
            else
               typ_var(nbvars) = mode2
            endif 
            goto 1 
c
c
         else
            goto 88888
         endif
c 
      else if ( ibloc .eq. 8 ) then
c        --------------------------------------------------------------
c        --- COMMON ---  
c        --------------------------------------------------------------
         ipos = ld + 1
         if ( buf(ipos:ipos) .eq. '/' ) then 
            ipos = ipos + 1
            ltok=0
            dowhile ( buf(ipos:ipos) .ne. '/' ) 
               ltok=ltok+1
               tok(ltok:ltok)=buf(ipos:ipos)
               ipos = ipos + 1
            enddo
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal +ipos - ltok
c  
            call set_mode(tok(1:ltok),
     +                       buf(min(ipos,llig):llig), mode2)
            typ_var(nbvars) = external
            mode = save
            ipos = ipos + 1
         endif
c        --- debut des varaibles ---
    2    continue 
         if ( ipos .gt. llig ) goto 88888
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
c
         if ( ltok .gt. 0 ) then
            ipos = ipos+idpos 
            if (buf(min(llig,ipos-1):min(llig,ipos-1)).eq.'(' ) then
c              --- c'est une declaration de tableau ---
               dctab = .true.
            else
               dctab = .false.
            endif
c
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal + ipos - ltok -1
c  
c           --- type implicite ---  
            if ( ichar(tok(1:1)) .ge. 73 .and.
     +           ichar(tok(1:1)) .le. 78    ) then
               mode = ivar
            else
               mode = rvar
            endif
c
            if ( dctab ) then
               typ_var(nbvars) = mode + 100
            else
               typ_var(nbvars) = mode
            endif 
            goto 2 
c
         else
            goto 88888
         endif
c
      else if ( ibloc .eq. 7 ) then
c        --------------------------------------------------------------
c        --- PARAMETER ---  
c        --------------------------------------------------------------
         ipos = ld + 1
         if ( buf(ipos:ipos) .eq. '(' ) then 
            ipos = ipos + 1
cx            dowhile ( .not. alpha(buf(ipos:ipos)) ) 
cx               ipos = ipos + 1
cx            enddo
         endif
c        --- debut des varaibles ---
    7    continue 
         if ( ipos .gt. llig ) goto 88888
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
         ipos = ipos+idpos -1
c
         if ( ltok .gt. 0 ) then
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal +ipos - ltok
c
cx            dowhile ( .not. alpha(buf(ipos:ipos)) 
cx     +                .and. ipos.lt.llig) 
cx               ipos = ipos + 1
cx            enddo 
c  
c           --- type implicite ---  
            if ( ichar(tok(1:1)) .ge. 73 .and.
     +           ichar(tok(1:1)) .le. 78    ) then
               mode = ivar
            else
               mode = rvar
            endif
            typ_var(nbvars) = mode
            goto 7 
         else
            goto 88888
         endif
c
      else if ( ibloc  .eq. 2 ) then
c        --------------------------------------------------------------
c        --- DATA --- 
c        --------------------------------------------------------------
         ipos = ld+1
122      continue
         if ( ipos .gt. llig ) goto 88888
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
         ipos = ipos+idpos-1 
         if ( ltok .eq. 0 ) then
            goto 88888
         else
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal+ipos-ltok
c  
            if ( ichar(tok(1:1)) .ge. 73 .and.
     +           ichar(tok(1:1)) .le. 78    ) then
               mode = ivar
            else
               mode = rvar
            endif
            typ_var(nbvars) = mode + 100
            goto 122
         endif
c
c
      else if ( ibloc .eq. 3 ) then
c        --------------------------------------------------------------
c        --- IF( --- 
c        --------------------------------------------------------------
         lif = lenstr(buf(1:llig)) 
         if ( mode .eq. if ) then
            call end_par(buf(3:llig),iposit) 
         else
            iposit = 1
         endif 
c
c         print*,'ibloc=3   "',buf(1:llig),'" iposit=',iposit
         if (( mode.eq.if .and. 
     +         buf(2+iposit:min(iposit+6,llig)).eq.')THEN') .or.
     +       mode .eq. dowhile .or.
     +       mode .eq. save .or.
     +       mode .eq. call .or.
     +       mode .eq. elseif ) then

c           --------------------------------------------------------------
c           --- BLOC IF  ou DOWHILE ou SAVE--- 
c           --------------------------------------------------------------
            ipos = ld + 1 
            if ( lif .gt. 4 ) then
               if ( buf(lif-4:lif) .eq. ')THEN'  )  llig = lif - 4
            endif
   13       continue 
c            print*,'   "',buf(ipos:llig),'"',' ipos=',ipos
c            print*,'13: decal=',idecal,' ipos=',ipos,' llig=',llig
            if ( ipos .gt. llig ) goto 88888
            call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
            ipos = ipos+idpos-1 
c
            if ( ltok .gt. 0 ) then
c               print*,'token="',tok(1:ltok),'"'
               nbvars = nbvars + 1
               nom_var(nbvars) = ' '
               nom_var(nbvars) = tok(1:ltok)
               adr_var(nbvars) = idecal+ ipos - ltok
c  
               call set_mode(tok(1:ltok),
     +                       buf(min(ipos,llig):llig), mode2)
               if ( mode .eq. call ) then
                   typ_var(nbvars) = external
                   mode = save
               else
                  typ_var(nbvars) = mode2 
               endif

c
               goto 13 
c
            else
               goto 88888
            endif
c
         else
c           --------------------------------------------------------------
C           --- IF LOGIQUE --- 
c           --------------------------------------------------------------
            llig = lenstr(buf(1:llig))
            llig_save = llig
            if(buf(1:4).eq.'GOTO')then
              call end_par(buf(5:),ipos2)
              ipos2 = ipos2 + 4
              llig_new = llig - ipos2 + 4
              ipos = 5 
            else
              call end_par(buf(3:),ipos2)
              ipos2 = ipos2 + 2
              llig_new = llig - ipos2 + 2
              ipos = 3 
            endif
C           --- TRAITEMENT DU TEST --- 
            llig = ipos2 
            idecal = llig 
c            print*,'ipos=',ipos,' llig=',llig,' idecal=',idecal
   14       continue 
c            print*,'   "',buf(ipos:llig),'"',' ipos=',ipos
            if ( ipos .gt. llig ) then
               buf   = buf  (ipos2+1:) 
               ligne = ligne(ipos2+1:) 
               llig = llig_new
c               print*,'ligne="',ligne(1:llig),'"'
               goto 1000 {reanalyse si instruction avec mot cles( print)}
            endif

            call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
            ipos = ipos+idpos-1 
c
            if ( ltok .gt. 0 ) then
               nbvars = nbvars + 1
               nom_var(nbvars) = ' '
               nom_var(nbvars) = tok(1:ltok)
               adr_var(nbvars) = ipos - ltok
c              --- type implicite ---  
               call set_mode(tok(1:ltok),
     +                       buf(min(ipos,llig):llig), mode)
               typ_var(nbvars) = mode
               goto 14 
c
            else
               goto 14
            endif
         endif
c 
      else if ( ibloc .eq. 11 ) then
c        ---------------------------------------------------------------
c        SUBROUTINES OU FONCTIONS
c        ---------------------------------------------------------------
         ipos = ld+1   
   15    continue 
         if ( ipos .gt. llig ) goto 88888
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
         ipos = ipos+idpos -1
c
         if ( ltok .gt. 0 ) then
c
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal +ipos - ltok
c  
            call set_mode(tok(1:ltok),
     +                       buf(min(ipos,llig):llig), mode2)
            if ( mode .eq. call ) then
                typ_var(nbvars) = external
                mode = save
            else
                typ_var(nbvars) = mode2 
            endif
c
            goto 15 
c
         else
            goto 88888
         endif 
c
c
      else if ( ibloc .eq. 6 ) then
c        ---------------------------------------------------------------
c        AFFECTATION
c        ---------------------------------------------------------------
         ipos = ld + 1  
c
   17    continue 
         if ( ipos .gt. llig ) return
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
         ipos = ipos+idpos  - 1
c
         if ( ltok .gt. 0 ) then
c
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal +ipos - ltok
c  
            call set_mode(tok(1:ltok),buf(min(llig,ipos):llig),mode)
            typ_var(nbvars) = mode
c
            goto 17 
c
         else
            goto 88888
         endif 
c
c
      else if ( ibloc .eq. 12 .or. ibloc .eq. 13 ) then
c        ---------------------------------------------------------------
c        OPEN ou CLOSE ou 1ere partie de READ ou WRITE
c        ---------------------------------------------------------------
         ipos = 1  
c
   18    continue 
         if ( ipos .gt. llig ) then
            if ( ibloc .eq. 12 ) goto 88888
            if ( llig_save .eq. 0 ) goto 88888
            ligne(1:llig_save) = ligne(llig+1:llig+llig_save)
            buf  (1:llig_save) = buf  (llig+1:llig+llig_save)
            idecal = llig 
            llig = llig_save  
            ibloc = 3
            mode = save    {magouille!!!} 
            ld = 0 
            goto 1001
         endif
c
c         print*,'decal=',idecal,' tok:buf:"',buf(ipos:llig),'"'
         call token(buf(ipos:llig),llig-ipos+1,tok,ltok,idpos)
         ipos = ipos+idpos- 1
c
         if ( ltok .gt. 0 ) then
c            print*,'ipos=',ipos,' ltok=',ltok,' llig=',llig
c     +            ,' token="',tok(1:ltok),'"'
c            print*,' token="',tok(1:ltok),'" param="'
c     +                       ,buf(ipos-ltok:ipos),'"'
            call isa_mot_cle(buf(ipos-ltok:min(ipos,llig)),oui)
            if ( oui ) then 
               ipos = ipos + 1
               goto 18
            endif
            nbvars = nbvars + 1
            nom_var(nbvars) = ' '
            nom_var(nbvars) = tok(1:ltok)
            adr_var(nbvars) = idecal +ipos - ltok
c  
            call set_mode(tok(1:ltok),
     +                       buf(min(ipos,llig):llig), mode)
            typ_var(nbvars) = mode
c
            goto 18 
c
         else
c            print*,'token=""'
            if ( ibloc .eq. 12 ) goto 88888
            if ( llig_save .eq. 0 ) goto 88888
            ligne(1:llig_save) = ligne(llig+1:llig+llig_save)
            buf  (1:llig_save) = buf  (llig+1:llig+llig_save)
c            print*,' decal=',idecal,' llig=',llig
            idecal = idecal+llig 
            llig = llig_save
            ibloc = 3
            mode = save    {magouille!!!} 
            ld = 0
c            print*,' decal=',idecal,' buf:"',buf(1:llig),'"'
            goto 1001
         endif 
c
c
      else if ( ibloc .eq. 14)then  { goto }
c
      else 
C 888888**********
         goto 88888
      endif
88888 continue
c      print*,'buf="',buf(1:lbuf),'"'
c      print*,'     12345678901234567890123456789012345678901234567890'//
c     +   '12345678901234567890123456789012345678901234567890'
c      do i=1,nbvars
c        print*,i,' typ_var=',typ_var(i),' adr_var=',adr_var(i)
c     +        ,' nom_var="',nom_var(i)(1:6)
c      enddo
      end





      subroutine set_mode(tok,ligne,mode)
      character*(*) tok,ligne
      integer mode
c ......................................................................
c but donner le type a priorie (entier reel character fonction)
c ......................................................................
c entree:
c     tok: l'identificateur
c     ligne : suite dans l'instruction
c sortie:
c     mode
c ......................................................................
c A.Golgolab
c ......................................................................
c
%include'params.incl' 
c 
      ltok = len(tok)
      llig = len(ligne)
c
      if ( ligne(1:1) .eq. '(' ) then
         call end_par(ligne(1:llig),ipos)
         ind = indexx(ligne(1:ipos),':')
         if ( ind .eq. 0 ) then 
            mode = ext_tab
         else 
c           --- verif plus fine ---
            nbp = 0
            il = 2
            dowhile ( il .lt. ipos )
               if ( ligne(il:il) .eq. '(' ) nbp = nbp + 1
               if ( ligne(il:il) .eq. ')' ) nbp = nbp - 1
               if ( nbp .eq. 0 .and. ligne(il:il) .eq. ':' ) then
                  mode = character 
                  return
               endif 
               il = il + 1
            enddo
            mode = ext_tab
         endif
      else
         if ( ichar(tok(1:1)) .ge. 73 .and.
     +        ichar(tok(1:1)) .le. 78    ) then
            mode = ivar
         else
            mode = rvar
         endif 
      endif 
      end


      subroutine supr_blancs_debut(mot)
      character*(*) mot
c
      if ( mot .eq. ' ' ) return
      long = len(mot)
c
    2 continue
      if ( mot(1:1) .eq. ' ' ) then
         mot(1:long) = mot(2:long)//' '
         long = long -1
         goto 2
      endif  
      end
      subroutine token(buf,llig,tok,ltok,newpos)
      character*(*) buf,tok
      integer ltok
c ======================================================================
c extrait une suite de caracteres alphanu pouvant correspondre a
c     un identificateur
c ======================================================================
c entree:
c     ligne : ligne courante sera lue a partir du debut
c     llig  : longueur de la ligne
c sortie:
c     tok   : identificateur lu
c     ltok  : longueur de l'identificateur
c     newpos : nouvelle position dans la ligne (caractere special)
c ======================================================================
c     A.Golgolab
c ======================================================================
c 
      logical alphanu,alpha,numerique
      external alphanu,alpha,numerique
      character*1440 ligne

      l = 1  
      ligne=buf
c 
    1 continue 
      if (l.gt.llig) then
         newpos = llig
         ltok = 0
         goto 999
      endif
      if ( ligne(l:l) .eq. '''') then
c        --- on est dans une chaine de caractere --- 
         call find_string(ligne(l:),ipos)
         l = l + ipos 
      endif
c
      if ( l .ge. llig .and. ligne(llig:llig) .eq. '''') then
         newpos=l
         ltok = 0
         goto 999
      endif
      dowhile ( .not.alpha(ligne(l:l)) .and. l.lt.llig .and.
     +          ligne(l:l) .ne. '''' ) 
c        on saute par dessus les non alpha
      if(ligne(l:l).eq.'.') then
          if(   ligne(l:l+3).eq.'.EQ.'
     +      .or.ligne(l:l+3).eq.'.GT.'
     +      .or.ligne(l:l+3).eq.'.LT.'
     +      .or.ligne(l:l+3).eq.'.OR.'
     +      .or.ligne(l:l+3).eq.'.NE.'
     +      .or.ligne(l:l+3).eq.'.GE.'
     +      .or.ligne(l:l+3).eq.'.LE.' )then
            l=l+4 
            goto 1
        elseif(ligne(l:l+4).eq.'.NOT.'
     +  .or.ligne(l:l+4).eq.'.AND.'
     +  .or.ligne(l:l+4).eq.'.EQV.'
     +  )then
          l=l+5 
          goto 1
        elseif(ligne(l:l+5).eq.'.NEQV.' 
     +  .or.ligne(l:l+5).eq.'.TRUE.'
     +  )then
           l=l+6
           goto 1
        elseif(ligne(l:l+6).eq.'.FALSE.' )then
          l=l+7
          goto 1
        endif
      endif
      l = l + 1
      enddo 
c
      if(    ligne(l:l).eq.'E'.or.ligne(l:l).eq.'D'
     +   .or.ligne(l:l).eq.'e'.or.ligne(l:l).eq.'d')then
        if(l.gt.2)then
          if(numerique(ligne(l-1:l-1))
     +    .or.(ligne(l-1:l-1).eq.'.'.and.numerique(ligne(l-2:l-2)))
     +    )then
            if(l.lt.llig)then
              if(ligne(l+1:l+1).eq.'+'.or.ligne(l+1:l+1).eq.'-'.or.
     +        numerique(ligne(l+1:l+1)) )then
                l=l+1
                goto 1
              endif
            endif
          endif
        endif
      endif
      if ( ligne(l:l) .eq. '''' ) goto 1
c
      if ( l .ge. llig .and. .not. alpha(ligne(l:l)) ) then
         newpos = llig+1
         ltok = 0
         goto 999
      endif
c     --- a la position l un caractere alphabetique 
      if ( l .gt. 1 ) then
         if ( ligne(l-1:l-1) .eq. '''' ) then
c           --- on est dans une chaine de caractere --- 
            l = l - 1
            goto 1
         endif
      endif
c
      tok(1:1) = ligne(l:l)
      ltok = 1
      l = l + 1 
      if ( l.gt.llig)  goto 2
c
      dowhile (  alphanu(ligne(l:l)) )
         ltok = ltok + 1
         tok(ltok:ltok) = ligne(l:l)
         l = l + 1
         if ( l.gt.llig)  goto 2
      enddo   
c
    2 continue   
      newpos=l
999   continue
c      if(ltok.gt.0)then
c         print*,'token "',tok(1:ltok),'"'
c      else
c         print*,'token ""'
c      endif   
      end
