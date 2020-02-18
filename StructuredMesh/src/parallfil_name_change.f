!***********************************************************************

!     FILTER OPERATIONS
      
!***********************************************************************
!c-----------------------------------------------------------------------
!      subroutine filter_init
!      include 'SIZE'
!      include 'INPUT'           ! IF3D
!      include 'GEOM'
!      include 'MASS'
!      include 'LEVELSET'        ! IBM_filter_R
!      include 'GHOSTBUSTER'
!      include 'FILTERSUPPORT'
!      integer gsh_ext
!      integer nlayers, ntot_ext
!      integer iel,ixyz,nxyz     !--> just for testing
!      real field_in(lx1,ly1,lz1,lelt)
!      real elemsize
!      
!!     calculate how many layers i need
!      elemsize= abs(XM1(NX1,1,1,1)-XM1(1,1,1,1))
!      nlayers = int(IBM_filter_R/elemsize)+1
!         
!      if (nlayers.gt.llay) then
!         if(NID.eq.0) write(*,*)'Filter radius=',IBM_filter_R
!         if(NID.eq.0) write(*,*)'Elem size=',elemsize
!         if(NID.eq.0) write(*,*)'Layer needed=',nlayers
!         if(NID.eq.0) write(*,*)'Max layer allowed=',llay
!         if(NID.eq.0) write(*,*)'Increase llay in GHOSTBUSTER!'
!         call exitt()
!      endif
!!     define the communication between processes for knowledge about the ghost cells
!      call handle_ext(gsh_ext,nlayers)
!!     expand the coordinates
!      call expand_data(x_ext,gsh_ext,XM1)
!      call expand_data(y_ext,gsh_ext,YM1)
!      if (IF3D) call expand_data(z_ext,gsh_ext,ZM1)
!      call expand_data(bm1_ext,gsh_ext,BM1)
!      
!!     save in FILTERSUPPORT the handle used for this specific extension
!      hndl_ext=gsh_ext
!
!      return
!      end
c-----------------------------------------------------------------------
!     define the filter
      subroutine blob(field_fil,x,y,z)
      implicit none
      include 'SIZE'
      include 'INPUT'           ! IF3D
      include 'TSTEP'           ! IF3D
      include 'ADJOINT'         ! IFADJ
      include 'LEVELSET'        ! IBM_filter_R, KERNEL_ !1-top-hat,2-bspline
      include 'GHOSTBUSTER'
      include 'FILTERSUPPORT'
!     argument list
      real x,y,z
      real field_fil
!     local variable
      integer ix,iy,iz,iel
      real distance,xdiff,ydiff,zdiff
      real fkernel, norm
!-----------------------------------------------------------------------
      field_fil = 0.0
      norm = 0.0
      do iel=1,nelt_ext 
         do iz=1,NZ1
            do iy=1,NY1
               do ix=1,NX1
                  xdiff = (x-x_ext(ix,iy,iz,iel))**2
                  ydiff = (y-y_ext(ix,iy,iz,iel))**2
                  if (IF3D) then
                     zdiff    = (z-z_ext(ix,iy,iz,iel))**2
                     distance = xdiff+ydiff+zdiff
                  else
                     distance = xdiff+ydiff
                  endif
                  if (IBM_KERNEL.eq.1) then
                     distance=abs(sqrt(distance)/(IBM_filter_R))
                     if (distance.le.1) then
                        fkernel=1.0
                     else
                        fkernel=0.0
                     endif
                  elseif (IBM_KERNEL.eq.2) then
                     distance=abs(1.5*(sqrt(distance)/(IBM_filter_R)))
                     if (distance.le.0.5) then
                        fkernel=1.0 - (distance**2/0.75)
                     elseif(distance.gt.0.5.and.distance.le.1.5)then
                        fkernel=(0.5*(1.5-distance)**2/0.75)
                     else
                        fkernel=0.0
                     endif
                  endif
                  norm = norm + fkernel*bm1_ext(ix,iy,iz,iel)
                  if (.not.IFADJ) then !direct filter operation
                     field_fil = field_fil + fkernel
     $                    *field_ext(ix,iy,iz,iel)
     $                    *bm1_ext(ix,iy,iz,iel)
                  else !adjoint filter operation
                     field_fil = field_fil + fkernel
     $                    *field_ext(ix,iy,iz,iel)
     $                    *bm1_ext(ix,iy,iz,iel)
     $                    *C1weights_ext(ix,iy,iz,iel)
                  endif
               enddo
            enddo
         enddo
      enddo
      if(.not.IFADJ) field_fil = field_fil/norm
      if(IF_FILINI) field_fil = 1.0/norm !if it's before the simulation starts, initializes the weigths
      return
      end
!***********************************************************************
      subroutine filter_weightave(field_in)
      implicit none
      include 'SIZE'
      include 'GEOM'            ! XM1,YM1
      include 'INPUT'           ! IF3D
      include 'SOLN'
      include 'GHOSTBUSTER'
      include 'FILTERSUPPORT'
!     argument list
      real field_in(nx1,ny1,nz1,nelt)
!     local variables
      integer ix, iy, iz, iel
      real field_fil
      real x(3)
!-----------------------------------------------------------------------
!     expand the data to filter
      call expand_data(field_ext,hndl_ext,field_in)

!     loop over ALL the elements
      do iel=1,NELT            
         do iz=1,NZ1
            do iy=1,NY1
               do ix=1,NX1
                  x(1) = XM1(ix,iy,iz,iel)
                  x(2) = YM1(ix,iy,iz,iel)
                  x(3) = 0.0
                  if (IF3D) x(3) = ZM1(ix,iy,iz,iel) 
!     kernel of the filter
                  call blob(field_fil,x(1),x(2),x(3))
!     save the filtered variable
                  field_in(ix,iy,iz,iel) = field_fil
               enddo
            enddo
         enddo
      enddo

      return
      end
!***********************************************************************
      subroutine filter_openclose(var)
      implicit none
      include 'SIZE'
      include 'PARALLEL'        ! GLLNID,GLLEL
      include 'GEOM'            ! XM1,YM1
      include 'MASS'            ! BM1
      include 'INPUT'           ! IF3D
      include 'LEVELSET'        ! IBM_MSKLS
      include 'SOLN'
!     argument list
      real var(1)
!     local variables
      integer ix, iy, iz, iel, nt, iii
      real filter(LX1,LY1,LZ1,LELT)
      real mskls_(LX1,LY1,LZ1,LELT)
      real x(3), tst(3), mass1, mass3,beta_
!     functions
      real glsc3, glsc2
!-----------------------------------------------------------------------
      nt = NX1*NY1*NZ1*NELT
      beta_ = IBM_FLTR_beta
      
      call copy(mskls_,var,nt)
      
!     erode:
!     f1(rho)
      do iii=1,nt
         mskls_(iii,1,1,1)=1/(mskls_(iii,1,1,1) + beta_)
      enddo
!     averaging
      call filter_weightave(mskls_)
      call copy(IBM_MSKF1,mskls_,nt)
!     g1(W*f1(rho))
      do iii=1,nt
         mskls_(iii,1,1,1)=(1./mskls_(iii,1,1,1)) - beta_
      enddo
      call copy(IBM_MSKG1,mskls_,nt)

!     dilate:
!     f2(rho)
      do iii=1,nt
         mskls_(iii,1,1,1)=1./((1.-mskls_(iii,1,1,1)) + beta_)
      enddo
!     averaging
      call filter_weightave(mskls_)
      call copy(IBM_MSKF2,mskls_,nt)
!     g2(W*f2(rho))
      do iii=1,nt
         mskls_(iii,1,1,1)=1.-(1./mskls_(iii,1,1,1) - beta_)
      enddo
      call copy(IBM_MSKLS,mskls_,nt)
      
      return
      end
!***********************************************************************
!     adjoint filter + localization
      subroutine adjfilter_openclose(var)
      implicit none
      include 'SIZE'
      include 'PARALLEL'        ! GLLNID,GLLEL
      include 'GEOM'            ! XM1,YM1
      include 'MASS'            ! BM1,BINVM1
      include 'INPUT'           ! IF3D
      include 'LEVELSET'        ! MSKLS
!      include 'TOPOLOOPD'
      
!     argument list
      real var(1)
!     local variables
      integer ix, iy, iz, iel, nt, iii
      real filter(LX1,LY1,LZ1,LELT)
      real mskls_(LX1,LY1,LZ1,LELT)
      real x(3), tst(3), mass1, mass3, beta_
!     functions
      real glsc4
!-----------------------------------------------------------------------
      nt = NX1*NY1*NZ1*NELT
      beta_ = IBM_FLTR_beta
      
      call copy(mskls_,var,nt)

!     dilate-1
      do iii=1,nt 
!     dJ/drhodilate * drhodilate/drhog2
         mskls_(iii,1,1,1)= mskls_(iii,1,1,1)/
     $        (IBM_MSKF2(iii,1,1,1))**2
      enddo
!     dJ/drhog1
      call filter_weightave(mskls_)
!     dJ/drho
      do iii=1,nt
         mskls_(iii,1,1,1)=
     $        mskls_(iii,1,1,1)/((1-IBM_MSKG1(iii,1,1,1))+beta_)**2
         
!     erode-1
!     dJ/drhoerode * drhoerode/drhog1
         mskls_(iii,1,1,1)= -mskls_(iii,1,1,1)/(IBM_MSKF1(iii,1,1,1))**2

      enddo
!     dJ/drhog2
      call filter_weightave(mskls_)
!     dJ/drhoerode
      do iii=1,nt
         mskls_(iii,1,1,1)=
     $        -mskls_(iii,1,1,1)/((IBM_MSKNF(iii,1,1,1))+beta_)**2
      enddo   
!     restrict to Omega_0 (localization 2)
      call col3(var,mskls_,IBM_LOCALIZATION,nt)
      
      return
      end
!***********************************************************************
      subroutine define_localization
      implicit none
      include 'SIZE'
      include 'GEOM'            ! XM1,YM1
      include 'INPUT'           ! IF3D
      include 'PARALLEL'           ! IF3D
      include 'LEVELSET'        ! IBM_LOCALIZATION,OPT_XMAX,OPT_YMIN,OPT_YMAX
                                ! IC_XMAX,C1_weights
      include 'GHOSTBUSTER'
      include 'FILTERSUPPORT'
!     argument list
      real field_in(nx1,ny1,nz1,nelt)
!     local variables
      integer ix, iy, iz, iel, nt
      real field_fil
      real x(3)
!-----------------------------------------------------------------------
      nt = NX1*NY1*NZ1*NELT
      call rzero(IBM_C1weights,nt)
!     loop over ALL the elements
      do iel=1,NELT            
         do iz=1,NZ1
            do iy=1,NY1
               do ix=1,NX1
                  x(1) = XM1(ix,iy,iz,iel)
                  x(2) = YM1(ix,iy,iz,iel)
                  x(3) = 0.0
                  if (IF3D) x(3) = ZM1(ix,iy,iz,iel) 
!     kernel of the filter
                  call blob(field_fil,x(1),x(2),x(3))
!     save the filtered variable
                  IBM_C1weights(ix,iy,iz,iel) = field_fil
               enddo
            enddo
         enddo
      enddo
      call expand_data(C1weights_ext,hndl_ext,IBM_C1weights)
      
      !define where the optimization box starts
      if(NID.eq.GLLNID(1)) then
         IBM_IC_XMAX=XM1(NX1,1,1,LGLEL(1))
      else
         IBM_IC_XMAX=0.0
      endif
      call gop(IBM_IC_XMAX,IBM_IC_XMAX,'+  ',1)
      
!     define the optimization domain
      call rone(IBM_LOCALIZATION,nt)
      do iel=1,NELT 
         do iz=1,NZ1
            do iy=1,NY1
               do ix=1,NX1
                  if (XM1(ix,iy,iz,iel).gt.IBM_OPT_XMAX
     $                 .or.XM1(ix,iy,iz,iel).lt.IBM_IC_XMAX) then
                     IBM_LOCALIZATION(ix,iy,iz,iel) = 0.0
                  endif
                  if(YM1(ix,iy,iz,iel).gt.IBM_OPT_YMAX
     $                 .or.YM1(ix,iy,iz,iel).lt.IBM_OPT_YMIN) then
                     IBM_LOCALIZATION(ix,iy,iz,iel) = 0.0
                  endif
                  if(ZM1(ix,iy,iz,iel).gt.0.25
     $                 .or.ZM1(ix,iy,iz,iel).lt.-0.25) then
                     IBM_LOCALIZATION(ix,iy,iz,iel) = 0.0
                  endif
               enddo
            enddo
         enddo
      enddo
      
      return
      end
!!***********************************************************************
!
!!     EXPAND THE PARTITION ON EACH PROCESS TO THE SURROUNDING GHOST CELLS
!
!c-----------------------------------------------------------------------
!      subroutine expand_data(field_out,
!     $     gsh_ext,field_in)
!      implicit none
!      include 'SIZE'
!      include 'GHOSTBUSTER'
!!     argument list
!      real field_out(1)
!      real field_in(1)
!      integer gsh_ext
!!     local variables
!      integer ntot_ext
!      ntot_ext = nx1*ny1*nz1*nelt_ext
!      call rzero(field_out,ntot_ext)
!      call copy (field_out,field_in,nx1*ny1*nz1*nelt)
!      call gs_op(gsh_ext,field_out,1,1,0) ! sum
!
!      return
!      end
!c-----------------------------------------------------------------------
!      subroutine handle_ext(gsh_ext,
!     $     nlay)
!      implicit none
!      include 'SIZE'
!      include 'PARALLEL'        ! lglel
!      include 'GHOSTBUSTER'
!!     argument list
!      integer nlay,gsh_ext
!!     local variables
!      integer ilay
!
!c     Initialize ieg_ext, lglel_ext, nelt_prev and nelt_ext
!      nelt_ext = nelt
!      call izero(lglel_ext,nelt_ext*(1+(3**ndim-1)*nlay))
!      call icopy(lglel_ext,lglel,nelt)
!      call init_nei
!
!c     Extend arrays by one layer
!      if (nlay .gt. llay) then
!         write(6,*) 'Increase llay'
!         call exitt()
!      endif
!
!c     Who you gonna call?
!      do ilay=1,nlay
!         call ghostbuster ! Extend partition with one layer of surrounding ghost elements
!      enddo
!
!c     Compute final gs handle
!      call gs_setup_ext(gsh_ext,nx1)
!
!      return
!      end
!c-----------------------------------------------------------------------
!      subroutine init_nei
!      implicit none
!      include 'SIZE'
!      include 'INPUT'           ! if3d
!      include 'PARALLEL'        ! lglel,nelgt
!      include 'GHOSTBUSTER'
!      
!c     local variables
!      integer ix,iy,iz,dx,dy,dz,glob_nei
!      integer iel,ieg,ic
!      integer nelx,nely,nelz
!
!      nelx = 22!nint(nelgt**(1./ndim))
!      nely = 22!nelx
!      nelz = 1
!      if (if3d) nelz = 11!nelx
!      ! CLIO'S VERSION
!
!      do iel=1,nelt
!         ic = 0
!         if (if3d) then
!            do dz = -1,1
!            do dy = -1,1
!            do dx = -1,1
!               ic = ic+1
!               ieg=lglel(iel)
!               if( (mod(ieg,nelx).eq.1.and.dx.eq.-1) .or.
!     $              (mod(ieg,nelx).eq.0.and.dx.eq.1) .or.
!     $              (ieg.ge.1.and.ieg.le.nelx*nely.and.dy.eq.-1) .or.
!     $              (ieg.ge.(nelx*nely*(nelz-1))+1.and.
!     $               ieg.le.nelx*nely*nelz.and.dy.eq.1) .or.
!     $              (mod(ieg,(nelx*nely)).ge.1
!     $              .and.mod(ieg,(nelx*nely)).le.nelx.and.dz.eq.-1).or.
!     $              (mod(ieg,(nelx*nely)).ge.(nelx*(nely-1))+1
!     $              .and.mod(ieg,(nelx*nely)).le.(nelx*nely)-1
!     $              .and.dz.eq.1)
!     $              .or.(mod(ieg,(nelx*nely)).eq.0.and.dz.eq.1) ) then
!                  ieg_nei(ic,iel)=0
!               else
!                  ieg_nei(ic,iel)=ieg+(nelx*nely*dz)+(nelx*dy)+dx
!               endif
!            enddo
!            enddo
!            enddo
!         else
!            do dy = -1,1
!            do dx = -1,1
!               ic = ic+1
!               ieg=lglel(iel)
!               if( (mod(ieg,nelx).eq.1.and.dx.eq.-1) .or.
!     $              (mod(ieg,nelx).eq.0.and.dx.eq.1) .or.
!     $              (ieg.ge.1.and.ieg.le.nelx.and.dy.eq.-1) .or.
!     $              (ieg.ge.(nelx*(nelx-1))+1.and.ieg.le.nelx**2
!     $              .and.dy.eq.1) ) then
!                  ieg_nei(ic,iel)=0
!               else
!                  ieg_nei(ic,iel)=ieg+(nelx*dy)+dx
!               endif
!            enddo
!            enddo
!         endif
!      enddo
!c$$$
!c$$$!     NICOLAS' VERSION
!c$$$      do iel=1,nelt
!c$$$         ic = 0
!c$$$         ieg = lglel(iel)
!c$$$         if (if3d) then
!c$$$            do dz = -1,1
!c$$$            do dy = -1,1
!c$$$            do dx = -1,1
!c$$$               ic = ic+1
!c$$$               ieg_nei(ic,iel) = glob_nei(ieg,dx,dy,dz,nelx,nely,nelz)
!c$$$            enddo
!c$$$            enddo
!c$$$            enddo
!c$$$         else
!c$$$            dz = 0
!c$$$            do dy = -1,1
!c$$$            do dx = -1,1
!c$$$               ic = ic+1
!c$$$               ieg_nei(ic,iel) = glob_nei(ieg,dx,dy,dz,nelx,nely,nelz)
!c$$$            enddo
!c$$$            enddo
!c$$$         endif
!c$$$      enddo
!c$$$
!c$$$      return
!c$$$      end
!c$$$c-----------------------------------------------------------------------
!c$$$      function glob_nei(ieg,dx,dy,dz,nx,ny,nz)
!c$$$c     Assuming a regular box mesh with nx*ny*nz elements, the routine
!c$$$c     computes the global number of an element shifted by (dx,dy,dz)
!c$$$c     from the element with global number ieg
!c$$$      implicit none
!c$$$      include 'SIZE_DEF'
!c$$$      include 'SIZE'
!c$$$      include 'INPUT_DEF'
!c$$$      include 'INPUT'
!c$$$      integer ieg,dx,dy,dz,nx,ny,nz,glob_nei
!c$$$      integer itmp,igx,igy,igz
!c$$$
!c$$$      igz = 1
!c$$$      if (if3d) then
!c$$$         igz = int((ieg-1)/(nx*ny))+1+dz
!c$$$      endif
!c$$$      itmp = mod(ieg-1,nx*ny)
!c$$$      igy = int(itmp/nx)+1+dy
!c$$$      igx = mod(ieg-1,nx)+1+dx
!c$$$
!c$$$      if ( igx .lt. 1 .or. igx .gt. nx .or.
!c$$$     $     igy .lt. 1 .or. igy .gt. ny .or.
!c$$$     $    (if3d .and. (igz .lt. 1 .or. igz .gt. nz))) then
!c$$$         glob_nei = 0
!c$$$      else
!c$$$         glob_nei = (igz-1)*nx*ny + (igy-1)*nx + (igx-1) + 1
!c$$$      endif
!
!      return
!      end
!c-----------------------------------------------------------------------      
!      subroutine ghostbuster
!      implicit none
!      include 'SIZE'
!      include 'INPUT'           ! if3d
!      include 'GHOSTBUSTER'
!
!c     local variables
!      integer iel,inei,ic,nnei,nelt_prev,nelt_nei,ntot_nei
!      save nelt_prev
!      data nelt_prev /1/
!
!c     gs data
!      integer gsh_ext
!
!      if (if3d) then
!         nnei = 27
!      else
!         nnei = 9
!      endif
!
!c     Update lglel_ext (extended list of elements I want)
!      nelt_nei = nelt_ext
!      do iel=nelt_prev,nelt_ext
!         do inei=1,nnei
!            if (ieg_nei(inei,iel).ne.0.and.
!     $           ieg_nei(inei,iel).ne.lglel_ext(iel)) then           ! if not 0 and not yours
!               do ic=1,nelt_nei
!                  if (lglel_ext(ic).eq.ieg_nei(inei,iel)) then ! if already in the list, do not add
!                     exit
!                  endif
!                  if (ic.eq.nelt_nei) then ! add to list
!                     nelt_nei = nelt_nei + 1
!                     lglel_ext(nelt_nei) = ieg_nei(inei,iel)
!                  endif
!               enddo
!            endif
!         enddo
!      enddo
!
!      nelt_prev = nelt_ext
!      nelt_ext = nelt_nei
!
!c     Build gs_setup for extended array
!      call gs_setup_ext(gsh_ext,3)
!
!c     GS operation to extend ieg_nei
!      ntot_nei = nnei*(nelt_ext-nelt_prev)
!      call izero(ieg_nei(1,nelt_prev+1),ntot_nei)
!      call gs_op(gsh_ext,ieg_nei,2,1,0) ! 2->integer, 1->sum, 0->does not matter
!      
!c     Free gs
!      call gs_free(gsh_ext)
!      
!      return
!      end
!c-----------------------------------------------------------------------
!      subroutine gs_setup_ext(gsh_ext,nx)
!      implicit none
!      include 'SIZE'
!      include 'INPUT'           ! if3d
!      include 'GHOSTBUSTER'
!      
!c     argument list
!      integer gsh_ext,nx
!
!c     Local variables
!      integer ny,nz,ntot_ext
!      integer mid,mp,nekcomm
!      common /nekmpi/ mid,mp,nekcomm
!
!c     Build global numbering
!      call build_glo_num_ext(nx)
!
!c     Use global numbering for the handle
!      ny = nx
!      nz = 1
!      if (if3d) nz = nx
!      
!      ntot_ext = nx*ny*nz*nelt_ext
!      call gs_setup(gsh_ext,glo_num_ext,ntot_ext,nekcomm,mp)            
!      
!      return
!      end      
!c-----------------------------------------------------------------------      
!      subroutine build_glo_num_ext(nx)
!      implicit none
!      include 'SIZE'
!      include 'INPUT'           ! if3d
!      include 'GHOSTBUSTER'
!      
!c     argument list
!      integer nx
!      
!c     Local variables
!      integer idx,ix,iy,iz,iel,ny,nz,nxy,nxyz
!
!c     Setup dimensions
!      ny = nx
!      nz = 1
!      if (if3d) nz = nx
!      nxy = nx*ny
!      nxyz = nx*ny*nz
!      
!c     Build global numbering
!      do iel=1,nelt_ext
!         do iz=1,nz
!            do iy=1,ny
!               do ix=1,nx
!                  idx = (iel-1)*nxyz + (iz-1)*nxy + (iy-1)*nx + ix           
!                  glo_num_ext(idx)=int8((lglel_ext(iel)-1)*nxyz
!     $                 + (iz-1)*nxy + (iy-1)*nx + (ix-1) + 1)                        
!               enddo
!            enddo
!         enddo
!      enddo      
!      
!      return
!      end
!c-----------------------------------------------------------------------
