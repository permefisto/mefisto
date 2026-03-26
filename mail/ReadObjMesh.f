c ----------------------------------------------------------------------
c DECLARATIONS OF VARIABLES and ARRAYS of a xyznpef MEFISTO MESH
c ----------------------------------------------------------------------
c     mxnode:  maximum number of nodes of the mesh
c     mxMater: maximum number of materials (for 1d  1 material= 1 line)
c     mxBC  :  maximum number of boundary conditions
c          (1d 1 boundary condition = 1 point)
c     mxtyfe:  maximum number of finite element types (here =1)
c     mxnofe:  maximum number of nodes         of a finite element
c     mxvpfe:  maximum number of vertex-point  of a finite element
c     mxelfe:  maximum number of edge-line     of a finite element
c     mxfsfe:  maximum number of face-surface  of a finite element
c     mxvvfe:  maximum number of volume-volume of a finite element
c     mxfe  :  maximum number of finite elements
c      parameter (mxnode=10 000, mxMater=16, mxBC=32, mxtyfe=3,
c     %           mxnofe=20, mxvpfe=8, mxelfe=12, mxfsfe=6, mxvvfe=1,
c     %           mxfe=20 000 )
c
c      real         xyznode(3,mxnode)
c
c     nuLineMater: Mefisto number of each line-material
c     nmLineMater: Mefisto name of each line-material
c      integer      nuLineMater(mxMater)
c      character*24 nmLineMater(mxMater)
c
c     nuPointBC: Mefisto number of each boundary-point
c     nmPointBC: Mefisto name   of each boundary-point
c      integer      nuPointBC(mxBC)
c      character*24 nmPointBC(mxBC)
c
c     nbfe  :  number of finite elements of the mesh
c     nbnofe:  number of nodes        of one finite element
c     nbvpfe:  number of vertex-point of one finite element
c     nbelfe:  number of edge-line    of one finite element
c     nbfsfe:  number of face-surface of one finite element
c     nbvvfe:  number of volume       of one finite element
c      integer  nbfe  ,
c     %         nbnofe,
c     %         nbvpfe,
c     %         nbelfe,
c     %         nbfsfe,
c     %         nbvvfe
c     nunofe:  number of each node         of each finite element
c     nuvpfe:  number of each vertex-point of each finite element
c     nuelfe:  number of each edge-line    of each finite element
c     nufsfe:  number of each face-surface of each finite element
c     nuvvfe:  number of the volume        of each finite element
c      integer  nunofe(mxnofe,mxfe),
c     %         nuvpfe(mxvpfe,mxfe),
c     %         nuelfe(mxelfe,mxfe),
c     %         nufsfe(mxfsfe,mxfe),
c     %         nuvvfe(mxvvfe,mxfe)
      subroutine ReadObjMesh( NameMesh, mxnode,  nbnode,  xyznode,
     %                        mxMater, nbMater,nuLineMater, nmLineMater,
     %                        mxBC,     nbBC,    nuPointBC,   nmPointBC,
     %                        mxfe,     nbfe,
     %                        mxnofe,   nbnofe,  nunofe,
     %                        mxvpfe,   nbvpfe,  nuvpfe,
     %                        mxelfe,   nbelfe,  nuelfe,
     %                        mxfsfe,   nbfsfe,  nufsfe,
     %                        mxvvfe,   nbvvfe,  nuvvfe,
     %                        ierr )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Goal : Read from a xyznpef Mefisto file all characteristics of a mesh
C ------
c input :
c -------
c NameMesh: name of the file of the mesh
c mxnode:  maximum number of nodes of the mesh
c mxMater: maximum number of materials (for 1d  1 material= 1 line)
c mxBC  :  maximum number of boundary conditions
c          (1d: 1 boundary condition = 1 point)
c mxnofe:  maximum number of nodes        of one finite element
c mxvpfe:  maximum number of vertex-point of one finite element
c mxelfe:  maximum number of edge-line    of one finite element
c mxfsfe:  maximum number of face-surface of one finite element
c mxvvfe:  maximum number of volume       of one finite element
c
c output:
c -------
c nbnode:   number of nodes of the mesh
c xyznode:  3 coordinates of the nbnode nodes
c
c nbMater: number of materials
c nuLineMater: Mefisto Line of each Line-Material
c nmLineMater: Mefisto name of each Line-Material
c
c nbBC: number of boundary conditions
c nuPointBC: Mefisto Point number of each boundary-point
c nmPointBC: Mefisto Point name   of each boundary-point
c
c nbtyfe:  number of finite element types (here =1)
c
c mxfe  :  maximum number of finite elements of the mesh
c nbfe  :  number of finite elements of the mesh
c nbnofe:  number of nodes           of a finite element
c nbvpfe:  number of vertex-point    of a finite element
c nbelfe:  number of edge-line       of a finite element
c nbfsfe:  number of face-surface    of a finite element
c nbvvfe:  number of volume          of a finite element
c
c nunofe:  number of each node         of each finite element
c nuvpfe:  number of each vertex-point of each finite element
c nuelfe:  number of each edge-line    of each finite element
c nufsfe:  number of each face-surface of each finite element
c nuvvfe:  number of the volume        of each finite element
c
c ierr  :  error code number 0 without error, else not zero
C23456---------------------------------------------------------------012
      character*24 NameMesh
      integer      ierr
      real         xyznode(3,mxnode)
      integer      nuLineMater(mxMater)
      character*24 nmLineMater(mxMater)
      integer      nuPointBC(mxBC)
      character*24 nmPointBC(mxBC)
      integer      nbfe  ,
     %             nbnofe,
     %             nbvpfe,
     %             nbelfe,
     %             nbfsfe,
     %             nbvvfe
      integer      nunofe(mxnofe,mxfe),
     %             nuvpfe(mxvpfe,mxfe),
     %             nuelfe(mxelfe,mxfe),
     %             nufsfe(mxfsfe,mxfe),
     %             nuvvfe(mxvvfe,mxfe)
c
      nfile = 10
      open( unit=nfile, file=NameMesh, status='old',
     %      access='sequential', form='formatted',
ccc     %      err=9001,
     %      iostat=ierr)
      if( ierr .ne. 0 ) then
         print *,'The file ',NameMesh,' can not open'
         goto 9001
      endif
c
      read(unit=nfile,fmt=*) NameMesh
      read(unit=nfile,fmt=*) nbcoor
      read(unit=nfile,fmt=*) nbdim
      read(unit=nfile,fmt=*) nbnode
      read(unit=nfile,fmt=*) nbtg
c
c     display these values
      print *, 'NameMesh=',NameMesh
      print *, 'nbcoor  =',nbcoor
      print *, 'nbdim   =',nbdim
      print *, 'nbnode  =',nbnode
      print *, 'nbtg    =',nbtg
c
c     the abscissa of the nbnode
      print *,'Coordinate x of',nbnode,' nodes:'
      do i=1, nbnode
         read (unit=nfile,fmt=*) (xyznode(k,i),k=1,3)
         print *, 'xyz',i,'=',(xyznode(k,i),k=1,3)
      enddo
c
c     the number of materials
      read (unit=nfile,fmt=*) nbMater
      print *,'Number of material lines=',nbMater
      do 10 k=1, nbMater
         read (unit=nfile,fmt=*) nuLineMater(k)
         read (unit=nfile,fmt=*) nmLineMater(k)
         print *,'Material',k,': line ',nuLineMater(k),
     %           ' of name ', nmLineMater(k)
 10   continue
c
c     the number of boundary condition support
      read (unit=nfile,fmt=*) nbBC
      print *, 'Number of Boundary Conditions =',nbBC
c
c     the number of points with a boundary condition
      read (unit=nfile,fmt=*) nbbcpo
      print *, 'Number of BC on Points=',nbbcpo
      do k=1, nbbcpo
         read (unit=nfile,fmt=*) nuPointBC(k)
         read (unit=nfile,fmt=*) nmPointBC(k)
         print *,'BC',k,': Point ',nuPointBC(k),
     %           ' of name ', nmPointBC(k)
      enddo
c
c     the number of finite element type
      read (unit=nfile,fmt=*) nbtyfe
      print *, 'Number of finite element type=',nbtyfe
c
c     loop on the different types of finite element
      do nutyfe = 1, nbtyfe
         read (unit=nfile,fmt=*) nbfe
         print *, 'Number of finite elements   =',nbfe
c
         read (unit=nfile,fmt=*) nbnofe
         print *, 'Node number     of this fe type=',nbnofe
c
         read (unit=nfile,fmt=*) nbvpfe
         print *, 'Vertices number of this fe type=',nbvpfe
         read (unit=nfile,fmt=*) nbelfe
         print *, 'Edges number    of this fe type=',nbelfe
         read (unit=nfile,fmt=*) nbfsfe
         print *, 'Faces number    of this fe type=',nbfsfe
         read (unit=nfile,fmt=*) nbvvfe
         print *, 'Volumes number  of this fe type=',nbvvfe
c
c        loop on the finite elements of this type
10050    FORMAT( 10I8 )
         do 50 nef = 1, nbfe
            print *
            print *,'Finite Element ',nef
            read (unit=nfile,fmt=10050)
     %           (nunofe(j,nef),j=1,nbnofe),
     %           (nuvpfe(j,nef),j=1,nbvpfe),
     %           (nuelfe(j,nef),j=1,nbelfe),
     %           (nufsfe(j,nef),j=1,nbfsfe),
     %           (nuvvfe(j,nef),j=1,nbvvfe)
            print *, 'fe nodes     ',
     %           (nunofe(j,nef),j=1,nbnofe)
            print *, 'vertex-point ',
     %           (nuvpfe(j,nef),j=1,nbvpfe)
            print *, 'edge-line    ',
     %           (nuelfe(j,nef),j=1,nbelfe)
            print *, 'face-surface ',
     %           (nufsfe(j,nef),j=1,nbfsfe)
            print *, 'volume       ',
     %           (nuvvfe(j,nef),j=1,nbvvfe)
 50      continue
      enddo
c
c     the total mesh number of finite elements
      read (unit=nfile,fmt=*) nbtofe
      print *, 'Total mesh number of finite elements=',nbtofe
c
 9001 return
      end
