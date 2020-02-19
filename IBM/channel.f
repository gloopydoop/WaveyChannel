#define NUMBER_ELEMENTS_X 10
#define NUMBER_ELEMENTS_Y 10

c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'


      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! FF[XYZ], U[XYZ]
!      include 'ADJOINT'
      include 'INPUT'
      include 'SOLN'
      include 'PARALLEL'
      
	FFX = 0.0
      FFY = 0.0
      FFZ = 0.0
	call IBM_forcing(FFX,FFy,FFZ,ix,iy,iz,ieg,ux,uy,uz)


      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

	QVOL   = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
	include 'LEVELSET'
      include 'MAP2D'
      include 'STATD'


	integer i, n
	integer nelx, nely 
	real u_max_center(lx1*NUMBER_ELEMENTS_X)
	real x_max_center(lx1*NUMBER_ELEMENTS_X), y_max_center(lx1*NUMBER_ELEMENTS_X) 
	real uo(lx1,ly1,lz1,lelv)
	real vo(lx1,ly1,lz1,lelv)
	real wo(lx1,ly1,lz1,lelv)
	real prc(lx2,ly2,lz2,lelv)
	real residu, L2, steady_tol
	real Re_tau_vol_flow, Re_tau_friction, Re_tau_dissipation
	real d_star(lx1,NUMBER_ELEMENTS_X), theta(lx1,NUMBER_ELEMENTS_X)
	real shape_factor(lx1,NUMBER_ELEMENTS_X)
	real tmp_real

	

	steady_tol = 0.00001

	nelx = NUMBER_ELEMENTS_X
	nely = NUMBER_ELEMENTS_Y

!	Frame stuff
c-----------------------------------------------------------------------
!     start framework
	if (ISTEP.eq.0) then
		residu = 99999
      	call frame_start
		call opcopy(uo,vo,wo,vx,vy,vz)
		call set_obj

!		Set up mask!
		
		call filter_init()
		call filter_openclose(IBM_MSKNF)
		call outpost(IBM_MSKNF,IBM_MSKLS,vz,pr,temp,'IC')

      endif

!     monitor simulation
      call frame_monitor

!     finalise framework
      if (ISTEP.eq.NSTEPS.or.LASTEP.eq.1) then
         call frame_end

      endif


! Stuff on every time step

c-----------------------------------------------------------------------
! Convergence check
! 	Check convergence to steady state
	call opsub2(uo,vo,wo,vx,vy,vz)
      call normvc(H1,SEMI,L2,LINF,uo,vo,wo)
	residu = L2/dt/param(55)
	call opcopy(uo,vo,wo,vx,vy,vz)
	write(6,*), residu
	if (residu.lt.steady_tol) then
		!call mntr_set_conv(.true.)
	endif

!	call my version of vol_flow
	if (istep.gt.0) call vol_flow_porous


	call compute_Re_tau_vol_flow(Re_tau_vol_flow)
	call compute_Re_tau_friction(Re_tau_friction)
	call compute_Re_tau_dissipation(Ree_tau_dissipation)

	write(6,*) 'Re_T'
	write(6,*) 'Friction		Dissipation		Vol_flow'
	write(6,*), Re_tau_friction, Ree_tau_dissipation, Re_tau_vol_flow

! Kill

	if (ISTEP.eq.800) then
		call outpost(vx,IBM_MSKNF,vz,pr,temp,'spy')
	endif

! Statistics on last timestep:
!	if (istep.eq.10) then ! (change this when you're ready!)
	write(6,*) 'Residual = ', residu
	if (residu.lt.steady_tol) then
!======================================================================
! Calculate Shape Factor
c-----------------------------------------------------------------------
	call stat_compute
      call stat_gs_sum(3)
	call stat_gs_max(4)
	! remember
	! 1) <1>
	! 2) <u>
	! 3) <uu>
	! 4) max(u)

	! Compute 
	open(unit=1000,file='Max_U_Adam.dat')
	open(unit=999 ,file='H_Adam.dat')
	write(1000,*) 'x_location			U_max'
	write(999 ,*) 'x_location			Shape Factor'
	do il = 1, map2d_lnum
		do i = 1, lx1
			tmp_real = stat_ruavg(i,il,2)/stat_ruavg(i,il,4)
			d_star(i,il) = stat_ruavg(i,il,1) - tmp_real

			tmp_real = stat_ruavg(i,il,4)*stat_ruavg(i,il,4)
			tmp_real = stat_ruavg(i,il,3)/tmp_real
			theta(i,il) = stat_ruavg(i,il,2)/stat_ruavg(i,il,4) - tmp_real
			shape_factor(i,il) = d_star(i,il)/theta(i,il)
			write(1000,*) map2d_xm1(i,1,il), stat_ruavg(i,il,4)
			write(999 ,*) map2d_xm1(i,1,il), shape_factor(i,il)
		enddo
	enddo
	close(1000)
	close(999)
! Calculate pressure drops
c-----------------------------------------------------------------------
	call compute_Re_tau_vol_flow(Re_tau_vol_flow)
	call compute_Re_tau_friction(Re_tau_friction)
	call compute_Re_tau_dissipation(Ree_tau_dissipation)

	call outpost(vx,vy,vz,pr,t,' ')

	write(6,*) 'Re_T'
	write(6,*) 'Friction		Dissipation		Vol_flow'
	write(6,*), Re_tau_friction, Ree_tau_dissipation, Re_tau_vol_flow

! Calculate max velocity line (old way)
c-----------------------------------------------------------------------
! 	Check the max velocity center line
	call planar_max_r(u_max_center,x_max_center,y_max_center,vx)

	open(unit=56,file='Max_U_Harry.dat')
	write(56,*) 'x_location			y_location			U_max '
	do i=1,lx1*NUMBER_ELEMENTS_X
		write(56,*) x_max_center(i), y_max_center(i), u_max_center(i)
	enddo
	close(56)


	call exitt0
	endif
!======================================================================
	
	


      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
c     NOTE ::: This subroutine MAY NOT be called by every process


C     Set boundary conditions

      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
	include 'LEVELSET'

!     argument list

c     velocity
      UX = 0.0
      UY = 0.0
      UZ = 0.0

!     Boundary conditions

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
	! implicit none

	integer ix,iy,iz,eg

C     Set initial conditions

      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
	include 'LEVELSET'

	call IBM_set_wave(ux,uy,uz,ix,iy,iz,ieg)

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

	real h, wave_amp, lam
	common /geom_harry/ h, wave_amp, lam

	ntot = nx1*ny1*nz1*nelt

	lam = 1.1111111110 ! HarryIsTheBomb
	wave_amp = .1111111111 ! HarryIsTheBigDog
	h = 0.5 ! BigDiggityHDizzle

c	param(54) = 0   ! x-direction
c	param(55) = 1 ! HarryDancesLikeAGod


	

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'
	include 'LEVELSET'

      param(55) = 0.0           ! turn off original code
      param(56) = 1.0           ! set ubar

      return
      end
c-----------------------------------------------------------------------

!======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     register modules
      call io_register
      call chkpt_register
	call IBM_register
!      call nseb_register
      call stat_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'SOLN'
!-----------------------------------------------------------------------
!     initialise modules
      call chkpt_init
	call IBM_init
!      call nseb_init
      call stat_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------

      
      return
      end subroutine
!======================================================================

c-----------------------------------------------------------------------
      subroutine planar_max_r(ua,xa,ya,u)
c
c     Compute s-t planar average of quantity u()
c
c 	HARRY: ANY MONEY YOU'VE FUCKED UP THE NELX VS LELX STUFF!
c	WATCH OUT IF YOU GO PARALLEL!
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'ZPER'
c
      real xa(nx1*NUMBER_ELEMENTS_X), ya(nx1*NUMBER_ELEMENTS_X) 
	real ua(nx1*NUMBER_ELEMENTS_X), u(nx1,ny1,nz1,nelv) 
      integer e,k,i,j,table_size, current_index
c
	call rzero(xa,nx1*NUMBER_ELEMENTS_X)
      call rzero(ua,nx1*NUMBER_ELEMENTS_X)
	call rzero(ya,nx1*NUMBER_ELEMENTS_X)

	table_size = 1
c	xa(1) = xm1(i,j,k,e)

c
      do e=1,nelv
	do k=1,nz1
	do i=1,nx1
	do j=1,ny1
		call check_index(current_index,xa,nx1*NUMBER_ELEMENTS_X,table_size,xm1(i,j,k,e))
            if(ua(current_index).lt.u(i,j,k,e)) then
			ua(current_index) = u(i,j,k,e)
			ya(current_index) = ym1(i,j,k,e)
		endif
         enddo
         enddo
         enddo
      enddo
c

      return
      end 


c-----------------------------------------------------------------------
	subroutine check_index(i,vec,max_n,n,x_val)
	
c	Give an array vec() of maximum size max_n but useful size n. 
c	Mapping indexes to unique reals.
c	Find the index ind for a searched value x_val.
c	If the element can't be found, append it to the end of the list
	real vec(max_n)
	real x_val
	integer n, i
	logical finished

	i = 1
	finished = .false.
	do while(i.le.n.and.(.not.finished))
		if(int(vec(i)*10.0**5).eq.int(x_val*10.0**5)) then 
			finished = .true.
		else
			i = i+1
		endif
	enddo

	if (.not.finished) then
		vec(i) = x_val
		n = n + 1
	endif
	
	return
	end
		

!======================================================================
      subroutine compute_Re_tau_dissipation(Re_tau)
      include 'SIZE'
      include 'TOTAL'

	real Re_d, F
	real lam, h, wave_amp
	common /geom_harry/ h, wave_amp, lam
	integer e, ntot, etot
	real sij(lx1*ly1*lz1,ldim,ldim), snrm(lx1*ly1*lz1,lelv)
	real Re_tau, total_dis, dpdx, rho, dnu, u_tau, vx_integrated, xsec, flow_rate, dis_unit_area

	etot = nx1*ny1*nz1
	ntot = etot*nelv
	do e=1,nelv
		call comp_gije(sij,vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),e)
		call comp_sije(sij)
		call mag_tensor_e(snrm(1,e),sij)
		call vsq(snrm(1,e),etot)
		call cmult(snrm(1,e),4.0,etot) ! tak to Philipp about this!!
	enddo
	
	total_dis = glsc2(snrm,bm1,ntot)

	Re_d = 1/param(2)
	flow_rate = 2*h*1
	dpdx = total_dis/flow_rate/lam
	F = total_dis*lam
	tw = dpdx
	u_tau  = sqrt(tw)
	Re_tau = Re_d/2*u_tau
      return
      end

!======================================================================

c-----------------------------------------------------------------------
       subroutine compute_Re_tau_vol_flow(Re_tau)
      include 'SIZE'
      include 'TOTAL'

!	parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)

!      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
!     $                , scale_vf(3)
!      common /cvflow_i/ icvflow,iavflow
!      common /cvflow_c/ chv(3)
!      character*1 chv


! These are from
      parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)
c
	common /cvflow_a/ vxc(kx1,ky1,kz1,lelv)
     $                , vyc(kx1,ky1,kz1,lelv)
     $                , vzc(kx1,ky1,kz1,lelv)
     $                , prc(kx2,ky2,kz2,lelv)
	common /cvflow_i/ icvflow,iavflow

	common /Hdizzle/ saved_scale

	real 			pr_added(kx2,ky2,kz2,lelv)

!	These are mine
	integer ntot2
	real dpdx, Re_tau, tw, rho, dnu, u_tau, base_flow, flow_rate, total_force, Re_d

	real lam, h, wave_amp
	common /geom_harry/ h, wave_amp, lam



	Re_d    = 1/param(2)
	total_force = saved_scale*volvm1
!	call BC_area(wall_area)
	rho    = 1.
	dnu    = param(2)
	tw = total_force/lam

	u_tau  = sqrt(tw/rho)
	Re_tau = u_tau*h/dnu



      return
      end
!======================================================================

c-----------------------------------------------------------------------
      subroutine compute_Re_tau_friction(Re_tau)
      include 'SIZE'
      include 'TOTAL'

	real rho, u_tau, Re_tau, dnu, Re_d
	real x0(3)
	real drag_total, total_wall_area, A1, A2, drag_equiv

	real lam, h, wave_amp
	common /geom_harry/ h, wave_amp, lam
	
c maybe not down here

	
	Re_d    = 1/param(2)
      call rzero(x0,3)               ! torque w.r.t. x0
	call torque_calc(1.0,x0,.false.,.false.) ! wall shear
!	dragx1 is bottom wall
!	dragx2 is top wall
	drag_total = dragx(1) + dragx(2)
!	call BC_area(total_wall_area)
!	A2 = lam ! top wall
!	A1 = total_wall_area - lam ! bottom wall
!	tw = drag_equiv/(lam)
!	tw = dragx(1)/(wall_area-lam) + dragx(2)/lam
	tw = drag_total/lam
	u_tau  = sqrt(tw)
	Re_tau = Re_d/2*u_tau

      return
      end
c-----------------------------------------------------------------------
	subroutine BC_area(BCarea)

      include 'SIZE'
      include 'TOTAL'


	integer nface, nxz
	real length

	nface = 2*ndim
	nxz   = nx1*nz1

	do e=1,nelv
		do f=1,nface
		   if (cbc(f,e,1).eq.'W  ') then
		      BCarea = BCarea + vlsum(area(1,1,f,e),nxz)
		endif
		enddo
	enddo

	return
      end


c-----------------------------------------------------------------------
      subroutine set_obj  ! define objects for surface integrals
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e,f
c
c     Define new objects
c
      nobj = 2			! for Periodic
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F' ! 'F'
         hcode( 2,ii) = 'F' ! 'F'
         hcode( 3,ii) = 'F' ! 'F'
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj
c
      if (maxobj.lt.nobj) write(6,*) 'increase maxobj in SIZEu. rm *.o'
      if (maxobj.lt.nobj) call exitt
c
      nxyz = nx1*ny1*nz1
      do e=1,nelv
      do f=1,2*ndim
         if (cbc(f,e,1).eq.'W  ') then
            iobj = 0
            if (f.eq.1) iobj=1  ! lower wall
            if (f.eq.3) iobj=2  ! upper wall
            if (iobj.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               ieg = lglel(e)
               object(iobj,mem,1) = ieg
               object(iobj,mem,2) = f
c              write(6,1) iobj,mem,f,ieg,e,nid,' OBJ'
    1          format(6i9,a4)
            endif
c
         endif
      enddo
      enddo
c     write(6,*) 'number',(nmember(k),k=1,4)
c
      return
      end
c-----------------------------------------------------------------------
!======================================================================
!> @brief Provide element coordinates and local numbers (user interface)
!! @param[out]  idir              mapping (uniform) direction
!! @param[out]  ctrs              2D element centres
!! @param[out]  cell              local element numberring
!! @param[in]   lctrs1,lctrs2     array sizes
!! @param[out]  nelsort           number of local 3D elements to sort
!! @param[out]  map_xm1, map_ym1  2D coordinates of mapped elements
!! @param[out]  ierr              error flag
      subroutine user_map2d_get(idir,ctrs,cell,lctrs1,lctrs2,nelsort,
     $     map_xm1,map_ym1,ierr)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! [XYZ]C
      include 'GEOM'            ! [XYZ]M1

!     argument list
      integer idir
      integer lctrs1,lctrs2
      real ctrs(lctrs1,lctrs2)  ! 2D element centres  and diagonals 
      integer cell(lctrs2)      ! local element numberring
      integer nelsort           ! number of local 3D elements to sort
      real map_xm1(lx1,lz1,lelt), map_ym1(lx1,lz1,lelt)
      integer ierr              ! error flag

!     local variables
      integer ntot              ! tmp array size for copying
      integer el ,il ,jl        ! loop indexes
      integer nvert             ! vertex number
      real rnvert               ! 1/nvert
      real xmid,ymid            ! 2D element centre
      real xmin,xmax,ymin,ymax  ! to get approximate element diagonal
      integer ifc               ! face number

!     dummy arrays
      real xcoord(8,LELT) ! tmp vertex coordinates

!#define DEBUG
#ifdef DEBUG
!     for testing
      character*3 str1, str2
      integer iunit, ierrl
      ! call number
      integer icalldl
      save icalldl
      data icalldl /0/
#endif

!-----------------------------------------------------------------------
!     initial error flag
      ierr = 0
!     set important parameters
!     uniform direction; should be taken as input parameter
!     x-> 1, y-> 2, z-> 3
      idir = 2
      
!     get element midpoints
!     vertex number
      nvert = 2**NDIM
      rnvert= 1.0/real(nvert)

!     eliminate uniform direction
      ntot = 8*NELV
      call copy(xcoord,XC,ntot) ! copy x


!     set initial number of elements to sort
      nelsort = 0
      call izero(cell,NELT)

!     for every element
      do el=1,NELV
!     element centre
         xmid = xcoord(1,el)
!     element diagonal
         xmin = xmid
         xmax = xmid
         do il=2,nvert
            xmid=xmid+xcoord(il,el)
            xmin = min(xmin,xcoord(il,el))
            xmax = max(xmax,xcoord(il,el))
         enddo
         xmid = xmid*rnvert

!     count elements to sort
            nelsort = nelsort + 1
!     1D position
!     in general this coud involve some curvilinear transform
            ctrs(1,nelsort)=xmid
            ctrs(2,nelsort)=0.0
!     reference distance
            ctrs(3,nelsort)=abs(xmax-xmin)
            if (ctrs(3,nelsort).eq.0.0) then
               ierr = 1
               return
            endif
!     element index
            cell(nelsort) = el
      enddo

!     provide 2D mesh
!     in general this coud involve some curvilinear transform
!     uniform y
      ifc = 1
      do el=1,nelv
         call ftovec(map_xm1(1,1,el),xm1,el,ifc,nx1,ny1,nz1)
      enddo
      il = lx1*lz1*nelv
      call rzero(map_ym1(1,1,el),il)

#ifdef DEBUG
!     testing
      ! to output refinement
      icalldl = icalldl+1
      call io_file_freeid(iunit, ierrl)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalldl
      open(unit=iunit,file='map2d_usr.txt'//str1//'i'//str2)
      
      write(iunit,*) idir, NELV, nelsort
      write(iunit,*) 'Centre coordinates and cells'
      do el=1,nelsort
         write(iunit,*) el, ctrs(:,el), cell(el)
      enddo
      write(iunit,*) 'GLL coordinates'
      do el=1,nelsort
         write(iunit,*) 'Element ', el
         write(iunit,*) 'XM1'
         do il=1,nz1
            write(iunit,*) (map_xm1(jl,il,el),jl=1,nx1)
         enddo
         write(iunit,*) 'YM1'
         do il=1,nz1
            write(iunit,*) (map_ym1(jl,il,el),jl=1,nx1)
         enddo
      enddo
      close(iunit)
#endif
#undef DEBUG
      return
      end subroutine
!=======================================================================
!> @brief Provide velocity, deriv. and vort. in required coordinates
!! @param[out]  lvel             velocity
!! @param[out]  dudx,dvdx,dwdx   velocity derivatives
!! @param[out]  vort             vorticity
      subroutine user_stat_trnsv(lvel,dudx,dvdx,dwdx,vort)
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'               ! if3d

!     argument list
      real lvel(LX1,LY1,LZ1,LELT,2) ! velocity array
      real dudx(LX1,LY1,LZ1,LELT,2) ! velocity derivatives; U
      real dvdx(LX1,LY1,LZ1,LELT,2) ! V
      real dwdx(LX1,LY1,LZ1,LELT,2) ! Z
      real vort(LX1,LY1,LZ1,LELT) ! vorticity

!     local variables
      integer itmp              ! dummy variable
      real rtmp(LX1,LY1,LZ1,LELT)
!-----------------------------------------------------------------------
!     Velocity transformation; simple copy
      itmp = NX1*NY1*NZ1*NELV
      call copy(lvel(1,1,1,1,1),VX,itmp)
      call copy(lvel(1,1,1,1,2),VY,itmp)

!     Derivative transformation
!     No transformation
      call gradm1(dudx(1,1,1,1,1),dudx(1,1,1,1,2),rtmp,lvel(1,1,1,1,1))
      call gradm1(dvdx(1,1,1,1,1),dvdx(1,1,1,1,2),rtmp,lvel(1,1,1,1,2))

!     get vorticity
!     curlz
      call sub3(vort(1,1,1,1),dvdx(1,1,1,1,1),dudx(1,1,1,1,2),itmp)

      return
      end subroutine
!======================================================================


















c-----------------------------------------------------------------------
      subroutine vol_flow_porous
c
c
c     Adust flow volume at end of time step to keep flow rate fixed by
c     adding an appropriate multiple of the linear solution to the Stokes
c     problem arising from a unit forcing in the X-direction.  This assumes
c     that the flow rate in the X-direction is to be fixed (as opposed to Y-
c     or Z-) *and* that the periodic boundary conditions in the X-direction
c     occur at the extreme left and right ends of the mesh.
c
c     pff 6/28/98
c
      include 'SIZE'
      include 'TOTAL'
	include 'LEVELSET'
c
c     Swap the comments on these two lines if you don't want to fix the
c     flow rate for periodic-in-X (or Z) flow problems.
c
      parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)
c
      common /cvflow_a/ vxc(kx1,ky1,kz1,lelv)
     $                , vyc(kx1,ky1,kz1,lelv)
     $                , vzc(kx1,ky1,kz1,lelv)
     $                , prc(kx2,ky2,kz2,lelv)
     $                , vdc(kx1*ky1*kz1*lelv,2)
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_i/ icvflow,iavflow
      common /cvflow_c/ chv(3)
      character*1 chv

	real mask_volume

	real saved_scale
	common /Hdizzle/ saved_scale
c
      real bd_vflow,dt_vflow
      save bd_vflow,dt_vflow
      data bd_vflow,dt_vflow /-99.,-99./

	real v_tmp(kx1,ky1,kz1,lelv), pr_tmp(kx2,ky2,kz2,lelv)

      logical ifcomp

c     Check list:

c     param (55) -- volume flow rate, if nonzero
c     forcing in X? or in Z?


      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if (param(55).eq.0.) return
      if (kx1.eq.1) then
         write(6,*) 'ABORT. Recompile vol_flow with kx1=lx1, etc.'
         call exitt
      endif

      icvflow   = 1                                  ! Default flow dir. = X
      if (param(54).ne.0) icvflow = abs(param(54))
      iavflow   = 0                                  ! Determine flow rate
      if (param(54).lt.0) iavflow = 1                ! from mean velocity
      flow_rate = param(55)

      chv(1) = 'X'
      chv(2) = 'Y'
      chv(3) = 'Z'

c     If either dt or the backwards difference coefficient change,
c     then recompute base flow solution corresponding to unit forcing:

      ifcomp = .false.
      if (dt.ne.dt_vflow.or.bd(1).ne.bd_vflow.or.ifmvbd) ifcomp=.true.
      if (.not.ifcomp) then
         ifcomp=.true.
         do i=1,ntot1
            if (vdiff (i,1,1,1,1).ne.vdc(i,1)) goto 20
            if (vtrans(i,1,1,1,1).ne.vdc(i,2)) goto 20
         enddo
         ifcomp=.false.  ! If here, then vdiff/vtrans unchanged.
   20    continue
      endif
      call gllog(ifcomp,.true.)
      
      call copy(vdc(1,1),vdiff (1,1,1,1,1),ntot1)
      call copy(vdc(1,2),vtrans(1,1,1,1,1),ntot1)
      dt_vflow = dt
      bd_vflow = bd(1)

      if (ifcomp) call compute_vol_soln_porous(vxc,vyc,vzc,prc)

c-------------------------HARRY------------------------------------------------
c 	could change it in the future to zeros in mask!!
      if (icvflow.eq.1) current_flow=glsc2(vx,bm1,ntot1)/domain_length  ! for X
      if (icvflow.eq.2) current_flow=glsc2(vy,bm1,ntot1)/domain_length  ! for Y
      if (icvflow.eq.3) current_flow=glsc2(vz,bm1,ntot1)/domain_length  ! for Z

      if (iavflow.eq.1) then
c-------------------------HARRY------------------------------------------------
	   mask_volume = glsc2(IBM_MSKLS,bm1,ntot1) ! this is the volume of the fluid
c-------------------------HARRY------------------------------------------------
         xsec = mask_volume / domain_length
         flow_rate = param(55)*xsec
      endif

      delta_flow = flow_rate-current_flow

c     Note, this scale factor corresponds to FFX, provided FFX has
c     not also been specified in userf.   If ffx is also specified
c     in userf then the true FFX is given by ffx_userf + scale.

      scale = delta_flow/base_flow
	saved_scale = scale
      scale_vf(icvflow) = scale
      if (nio.eq.0) write(6,1) istep,chv(icvflow)
     $   ,time,scale,delta_flow,current_flow,flow_rate
    1    format(i11,'  Harrys Volflow ',a1,11x,1p5e13.4)

c--------HARRY: THIS IS SO FUCKING SKETCHY, TALK TO PHILIPP!!!-----------
	call copy(v_tmp,vxc,ntot1)	
	call col2(v_tmp,MSKLS,ntot1)
      call add2s2(vx,v_tmp,scale,ntot1)

	call copy(v_tmp,vyc,ntot1)	
	call col2(v_tmp,MSKLS,ntot1)
      call add2s2(vy,vtmp,scale,ntot1)

	call copy(v_tmp,vzc,ntot1)	
	call col2(v_tmp,MSKLS,ntot1)
      call add2s2(vz,vtmp,scale,ntot1)

	call copy(pr_tmp,prc,ntot2)	
	call col2(pr_tmp,MSKLS,ntot2)
      call add2s2(pr,pr_tmp,scale,ntot2)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_vol_soln_porous(vxc,vyc,vzc,prc)
c
c     Compute the solution to the time-dependent Stokes problem
c     with unit forcing, and find associated flow rate.
c
c     pff 2/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
c
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_i/ icvflow,iavflow
      common /cvflow_c/ chv(3)
      character*1 chv
c
      integer icalld
      save    icalld
      data    icalld/0/
c
c
      ntot1 = lx1*ly1*lz1*nelv
      if (icalld.eq.0) then
         icalld=icalld+1
         xlmin = glmin(xm1,ntot1)
         xlmax = glmax(xm1,ntot1)
         ylmin = glmin(ym1,ntot1)          !  for Y!
         ylmax = glmax(ym1,ntot1)
         zlmin = glmin(zm1,ntot1)          !  for Z!
         zlmax = glmax(zm1,ntot1)
c
         if (icvflow.eq.1) domain_length = xlmax - xlmin
         if (icvflow.eq.2) domain_length = ylmax - ylmin
         if (icvflow.eq.3) domain_length = zlmax - zlmin
c
      endif
c
      if (ifsplit) then
c        call plan2_vol_porous(vxc,vyc,vzc,prc)
         call plan4_vol_porous(vxc,vyc,vzc,prc)
      else
         call plan3_vol_porous(vxc,vyc,vzc,prc)
      endif
c
c     Compute base flow rate
c 
      if (icvflow.eq.1) base_flow = glsc2(vxc,bm1,ntot1)/domain_length
      if (icvflow.eq.2) base_flow = glsc2(vyc,bm1,ntot1)/domain_length
      if (icvflow.eq.3) base_flow = glsc2(vzc,bm1,ntot1)/domain_length
c
      if (nio.eq.0 .and. loglevel.gt.2) write(6,1) 
     $   istep,chv(icvflow),base_flow,domain_length,flow_rate
    1    format(i11,'  basflow ',a1,11x,1p3e13.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine plan2_vol_porous(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (classical splitting scheme).
c
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
C
      COMMON /SCRNS/ RESV1 (LX1,LY1,LZ1,LELV)
     $ ,             RESV2 (LX1,LY1,LZ1,LELV)
     $ ,             RESV3 (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
c
      common /cvflow_i/ icvflow,iavflow
C
C
C     Compute pressure 
C
      ntot1  = lx1*ly1*lz1*nelv
c
      if (icvflow.eq.1) then
         call cdtp     (respr,v1mask,rxm2,sxm2,txm2,1)
      elseif (icvflow.eq.2) then
         call cdtp     (respr,v2mask,rxm2,sxm2,txm2,1)
      else
         call cdtp     (respr,v3mask,rxm2,sxm2,txm2,1)
      endif
c
      call ortho    (respr)
c
      call ctolspl  (tolspl,respr)
      call rone     (h1,ntot1)
      call rzero    (h2,ntot1)
c
      call hmholtz  ('PRES',prc,respr,h1,h2,pmask,vmult,
     $                             imesh,tolspl,nmxp,1)
      call ortho    (prc)
C
C     Compute velocity
C
      call opgrad   (resv1,resv2,resv3,prc)
      call opchsgn  (resv1,resv2,resv3)
      call add2col2 (resv1,bm1,v1mask,ntot1)
c
      intype = -1
      call sethlm   (h1,h2,intype)
      call ophinv   (vxc,vyc,vzc,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine plan3_vol_porous(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (PLAN3).
c
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
C
      COMMON /SCRNS/ rw1   (LX1,LY1,LZ1,LELV)
     $ ,             rw2   (LX1,LY1,LZ1,LELV)
     $ ,             rw3   (LX1,LY1,LZ1,LELV)
     $ ,             dv1   (LX1,LY1,LZ1,LELV)
     $ ,             dv2   (LX1,LY1,LZ1,LELV)
     $ ,             dv3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
      COMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)
      common /cvflow_i/ icvflow,iavflow
c
c
c     Compute velocity, 1st part 
c
      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      ifield = 1
c
      if (icvflow.eq.1) then
         call copy     (rw1,bm1,ntot1)
         call rzero    (rw2,ntot1)
         call rzero    (rw3,ntot1)
      elseif (icvflow.eq.2) then
         call rzero    (rw1,ntot1)
         call copy     (rw2,bm1,ntot1)
         call rzero    (rw3,ntot1)
      else
         call rzero    (rw1,ntot1)        ! Z-flow!
         call rzero    (rw2,ntot1)        ! Z-flow!
         call copy     (rw3,bm1,ntot1)    ! Z-flow!
      endif
      intype = -1
      call sethlm   (h1,h2,intype)
      call ophinv   (vxc,vyc,vzc,rw1,rw2,rw3,h1,h2,tolhv,nmxv)
      call ssnormd  (vxc,vyc,vzc)
c
c     Compute pressure  (from "incompr")
c
      intype = 1
      dtinv  = 1./dt
c
      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult   (h2,dtinv,ntot1)
      call invers2 (h2inv,h2,ntot1)
      call opdiv   (respr,vxc,vyc,vzc)
      call chsign  (respr,ntot2)
      call ortho   (respr)
c
c
c     Set istep=0 so that h1/h2 will be re-initialized in eprec
      i_tmp = istep
      istep = 0
      call esolver (respr,h1,h2,h2inv,intype)
      istep = i_tmp
c
      call opgradt (rw1,rw2,rw3,respr)
      call opbinv  (dv1,dv2,dv3,rw1,rw2,rw3,h2inv)
      call opadd2  (vxc,vyc,vzc,dv1,dv2,dv3)
c
      call cmult2  (prc,respr,bd(1),ntot2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine plan4_vol_porous(vxc,vyc,vzc,prc)

c     Compute pressure and velocity using fractional step method.
c     (Tombo splitting scheme).



      include 'SIZE'
      include 'TOTAL'

      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx1,ly1,lz1,lelv)

      common /scrns/ resv1 (lx1,ly1,lz1,lelv)
     $ ,             resv2 (lx1,ly1,lz1,lelv)
     $ ,             resv3 (lx1,ly1,lz1,lelv)
     $ ,             respr (lx1,ly1,lz1,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)

      common /cvflow_i/ icvflow,iavflow

      n = lx1*ly1*lz1*nelv
      call invers2  (h1,vtrans,n)
      call rzero    (h2,       n)

c     Compute pressure 

      if (icvflow.eq.1) call cdtp(respr,h1,rxm2,sxm2,txm2,1)
      if (icvflow.eq.2) call cdtp(respr,h1,rym2,sym2,tym2,1)
      if (icvflow.eq.3) call cdtp(respr,h1,rzm2,szm2,tzm2,1)

      call ortho    (respr)
      call ctolspl  (tolspl,respr)

      call hmholtz  ('PRES',prc,respr,h1,h2,pmask,vmult,
     $                             imesh,tolspl,nmxp,1)
      call ortho    (prc)

C     Compute velocity

      call opgrad   (resv1,resv2,resv3,prc)
      if (ifaxis) call col2 (resv2,omask,n)
      call opchsgn  (resv1,resv2,resv3)

      if (icvflow.eq.1) call add2col2(resv1,v1mask,bm1,n) ! add forcing
      if (icvflow.eq.2) call add2col2(resv2,v2mask,bm1,n)
      if (icvflow.eq.3) call add2col2(resv3,v3mask,bm1,n)


      if (ifexplvis) call split_vis ! split viscosity into exp/imp part

      intype = -1
      call sethlm   (h1,h2,intype)
      call ophinv   (vxc,vyc,vzc,resv1,resv2,resv3,h1,h2,tolhv,nmxv)

      if (ifexplvis) call redo_split_vis ! restore vdiff

      end
c-----------------------------------------------------------------------

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
