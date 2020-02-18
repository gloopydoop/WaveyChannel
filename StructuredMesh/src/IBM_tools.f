!=======================================================================
!     Clio
!     Tools for Immersed boundary fin
!     Parameters used by this set of subroutines:
!     LEVELSET:
!     DIFTIME, DIFSTEP, DIFK - diffusion of the mask boundary
!                              (smooth interface)
!     PHIMAX, PHIMIN - initial conditions for level set, max and min
!     AMPMSK1, AMPMSK2 - immersed boundary coefficients alpha(linear)
!                        and beta(integral)
!     IC_XMAX, IC_YMAX, IC_YMIN - initial geometry of the fin
!     IBM_PROP - coefficients for the governing equations
!=======================================================================
!***********************************************************************

!=======================================================================
!> @brief Register IBM module
!! @ingroup CLIO
!! @note This routine should be called in frame_usr_register
      subroutine IBM_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'LEVELSET'

      ! local variables
      integer lpmid, il
      character*2 str

!-----------------------------------------------------------------------
      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,IBM_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(IBM_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//'FRAME'//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(IBM_id,lpmid,IBM_name,
     $      'Conjugated heat transfer tools')

      ! register and set active section
      call rprm_sec_reg(IBM_sec_id,IBM_id,'_'//adjustl(IBM_name),
     $     'Runtime paramere section for IBM tool module')
      call rprm_sec_set_act(.true.,IBM_sec_id)

      ! register parameters
!-----------------------------------------------------------------------
      call rprm_rp_reg(IBM_AMPMSK_id,IBM_sec_id,'AMPMSK',
     $     'who knows',rpar_real,0,7250,.false.,' ')
      call rprm_rp_reg(IBM_OPT_XMAX_id,IBM_sec_id,'OPT_XMAX',
     $     'who knows',rpar_real,0,-0.25,.false.,' ')
      call rprm_rp_reg(IBM_OPT_YMAX_id,IBM_sec_id,'OPT_YMAX',
     $     'who knows',rpar_real,0,0.25,.false.,' ')
      call rprm_rp_reg(IBM_OPT_YMIN_id,IBM_sec_id,'OPT_YMIN',
     $     'who knows',rpar_real,0,-0.25,.false.,' ')
      call rprm_rp_reg(IBM_FILTER_R_id,IBM_sec_id,'FILTER_R',
     $     'who knows',rpar_real,0,0.02,.false.,' ')

! not sure about prop!!
      call rprm_rp_reg(IBM_P_111,IBM_sec_id,'P_111',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_112,IBM_sec_id,'P_112',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_113,IBM_sec_id,'P_113',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_121,IBM_sec_id,'P_121',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_122,IBM_sec_id,'P_122',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_123,IBM_sec_id,'P_123',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_211,IBM_sec_id,'P_211',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_212,IBM_sec_id,'P_212',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_213,IBM_sec_id,'P_213',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_221,IBM_sec_id,'P_221',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_222,IBM_sec_id,'P_222',
     $     'who knows',rpar_real,0,0,.false.,' ')
      call rprm_rp_reg(IBM_P_223,IBM_sec_id,'P_223',
     $     'who knows',rpar_real,0,0,.false.,' ')
! not sure about prop!!

      call rprm_rp_reg(IBM_KERNEL_id,IBM_sec_id,'KERNEL',
     $     'who knows',rpar_int,2,2,.false.,' ')
      call rprm_rp_reg(IBM_QC_ramp_id,IBM_sec_id,'QC_ramp',
     $     'who knows',rpar_real,0,8,.false.,' ')
      call rprm_rp_reg(IBM_QA_ramp_id,IBM_sec_id,'QA_ramp',
     $     'who knows',rpar_real,0,88,.false.,' ')
      call rprm_rp_reg(IBM_QF_ramp_id,IBM_sec_id,'QF_ramp',
     $     'who knows',rpar_real,0,8,.false.,' ')
      call rprm_rp_reg(IBM_FLTR_type_id,IBM_sec_id,'FLTR_type',
     $     'who knows',rpar_int,2,7250,.false.,' ')
      call rprm_rp_reg(IBM_FLTR_beta_id,IBM_sec_id,'FLTR_beta',
     $     'who knows',rpar_real,0,0.01,.false.,' ')
      call rprm_rp_reg(IBM_P_SIMP_id,IBM_sec_id,'P_SIMP',
     $     'who knows',rpar_int,0,5,.false.,' ')
      call rprm_rp_reg(IBM_MAP_RHO_id,IBM_sec_id,'MAP_RHO',
     $     'who knows',rpar_int,2,7250,.false.,' ')



      ! set initialisation flag
      IBM_ifinit=.false.

      return
      end subroutine
!=======================================================================
!> @brief Initilise IBM tools  module
!! @ingroup CLIO
!! @note This routine should be called in frame_usr_init
      subroutine IBM_init()
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'
      include 'SOLN'
      include 'LEVELSET'

      ! local variables
      integer itmp, il
      real rtmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (IBM_ifinit) then
         call mntr_warn(IBM_id,
     $        'module ['//trim(IBM_name)//'] already initiaised.')
         return
      endif

      ! get runtime parameters
!-----------------------------------------------------------------------
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_AMPMSK_id,rpar_real)
      IBM_AMPMSK = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_OPT_XMAX_id,rpar_real)
      IBM_OPT_XMAX = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_OPT_YMAX_id,rpar_real)
      IBM_OPT_YMAX = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_OPT_YMIN_id,rpar_real)
      IBM_OPT_YMIN = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_FILTER_R_id,rpar_real)
      IBM_FILTER_R = rtmp

! not sure about prop!!
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_111,rpar_real)
      IBM_PROP(1,1,1) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_112,rpar_real)
      IBM_PROP(1,1,2) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_113,rpar_real)
      IBM_PROP(1,1,3) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_121,rpar_real)
      IBM_PROP(1,2,1) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_122,rpar_real)
      IBM_PROP(1,2,2) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_123,rpar_real)
      IBM_PROP(1,2,3) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_111,rpar_real)
      IBM_PROP(2,1,1) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_212,rpar_real)
      IBM_PROP(2,1,2) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_213,rpar_real)
      IBM_PROP(2,1,3) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_221,rpar_real)
      IBM_PROP(2,2,1) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_222,rpar_real)
      IBM_PROP(2,2,2) = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_223,rpar_real)
      IBM_PROP(2,2,3) = rtmp
! not sure about prop!!

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_KERNEL_id,rpar_int)
      IBM_KERNEL = itmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_QC_ramp_id,rpar_real)
      IBM_QC_ramp = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_QA_ramp_id,rpar_real)
      IBM_QA_ramp = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_QF_ramp_id,rpar_real)
      IBM_QF_ramp = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_FLTR_type_id,rpar_int)
      IBM_FLTR_type = itmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_FLTR_beta_id,rpar_real)
      IBM_FLTR_beta = rtmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_P_SIMP_id,rpar_int)
      IBM_P_SIMP = itmp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,IBM_MAP_RHO_id,rpar_int)
      IBM_MAP_RHO = itmp



      ! everything is initialised
      IBM_ifinit=.true.

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup CLIO
!! @return IBM_is_initialised
      logical function IBM_is_initialised()
      implicit none

      include 'SIZE'
      include 'LEVELSET'
!-----------------------------------------------------------------------
      IBM_is_initialised = IBM_ifinit

      return
      end function








!!=======================================================================
!!     read parameters LEVELSET
!      subroutine IBM_param_in(fid)
!      implicit none
!
!      include 'SIZE'            !
!      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
!      include 'LEVELSET'        ! AMPMSK, IC_XMAX, IC_YMAX,
!                                ! IC_YMIN, IBM_PROP, FILTER_R, 
!                                ! OPT_XMAX, OPT_YMAX, OPT_YMIN
!
!!     argument list
!      integer fid               ! file id
!
!!     local variables
!      integer ierr, len, il, ip, lsent, ilen
!      parameter (lsent = (9+2*2*3))
!      parameter (ilen = 4)
!      real rtmp(lsent)
!      integer itmp(ilen)
!
!!     namelists
!      namelist /LEVELSET/ AMPMSK, IC_XMAX, IC_YMAX, IC_YMIN, FLTR_type
!     $     , OPT_XMAX, OPT_YMAX, OPT_YMIN, FILTER_R, IBM_PROP, KERNEL_
!     $     , QF_ramp, QA_ramp, QC_ramp, FLTR_type, FLTR_beta, P_SIMP
!     $     , MAP_RHO
!!-----------------------------------------------------------------------
!!     default values
!      AMPMSK  = 5.0000E+03          !IBM chi
!      OPT_XMAX= 0.25000
!      OPT_YMAX= 0.40000
!      OPT_YMIN=-0.40000
!      FILTER_R= 0.02500             !radius of the modal filter
!      IBM_PROP(1,1,1)= 0.01000      !cdiff fluid vel
!      IBM_PROP(1,1,2)= 1.41443      !ctrans fluid vel
!      IBM_PROP(1,1,3)= 0.00000      !uforce fluid vel
!      IBM_PROP(1,2,1)= 0.01000      !cdiff fluid temp 
!      IBM_PROP(1,2,2)= 1.00000      !ctrans fluid temp
!      IBM_PROP(1,2,3)= 0.00000      !qvol fluid temp
!      IBM_PROP(2,1,1)= 0.01000      !cdiff solid vel
!      IBM_PROP(2,1,2)= 1.41443      !ctrans solid vel
!      IBM_PROP(2,1,3)= 0.00000      !uforce solid vel
!      IBM_PROP(2,2,1)= 77.4908      !cdiff solid temp 
!      IBM_PROP(2,2,2)= 1132.64      !ctrans solid temp
!      IBM_PROP(2,2,3)= 0.00000      !qvol solid temp
!      KERNEL_ = 5                   !1-tophat, 2-b2spline, 3-gaussian,4-sfunction, 5-expwindow
!      QF_ramp = 88
!      QA_ramp = 8
!      QC_ramp = 8
!      FLTR_type=2               !1_average, 1_openclose
!      FLTR_beta=0.5
!      P_SIMP  = 3
!      MAP_RHO = 2
!         
!!     read the file
!      ierr=0
!      if (NID.eq.0) then
!         read(unit=fid,nml=LEVELSET,iostat=ierr)
!      endif
!      call err_chk(ierr,'Error reading LEVELSET parameters.$')
!
!!     broadcast data
!      if (NID.eq.0) then  
!         rtmp(1) = AMPMSK
!         rtmp(2) = OPT_XMAX 
!         rtmp(3) = OPT_YMAX 
!         rtmp(4) = OPT_YMIN
!         rtmp(5) = FILTER_R
!         ip = 5
!         do il=1,3
!            ip = ip +1
!            rtmp(ip)   = IBM_PROP(1,1,il)
!            rtmp(ip+3) = IBM_PROP(1,2,il)
!            rtmp(ip+6) = IBM_PROP(2,1,il)
!            rtmp(ip+9) = IBM_PROP(2,2,il)
!         enddo
!         itmp(1)    = KERNEL_
!         rtmp(18) = QF_ramp
!         rtmp(19) = QA_ramp
!         rtmp(20) = QC_ramp
!         itmp(2)    = FLTR_type
!         rtmp(21) = FLTR_beta
!         itmp(3)  = P_SIMP
!         itmp(4)  = MAP_RHO
!      endif
!      len = lsent*WDSIZE
!      call bcast(rtmp,len)
!      call bcast(itmp,ilen*ISIZE)
!      if (NID.ne.0) then
!         AMPMSK  = rtmp(1)
!         OPT_XMAX= rtmp(2)
!         OPT_YMAX= rtmp(3) 
!         OPT_YMIN= rtmp(4) 
!         FILTER_R= rtmp(5)
!         ip = 5
!         do il=1,3
!            ip = ip +1
!            IBM_PROP(1,1,il) = rtmp(ip) 
!            IBM_PROP(1,2,il) = rtmp(ip+3) 
!            IBM_PROP(2,1,il) = rtmp(ip+6) 
!            IBM_PROP(2,2,il) = rtmp(ip+9) 
!         enddo
!         KERNEL_ = itmp(1)
!         QF_ramp = rtmp(18)
!         QA_ramp = rtmp(19)
!         QC_ramp = rtmp(20)
!         FLTR_type = itmp(2)
!         FLTR_beta = rtmp(21)
!         P_SIMP  = itmp(3)
!         MAP_RHO = itmp(4)
!      endif
!
!      return
!      end
!!***********************************************************************
!!     write parameters LEVELSET
!      subroutine IBM_param_out(fid)
!      implicit none
!
!      include 'SIZE'            ! 
!      include 'LEVELSET'        ! AMPMSK, IC_XMAX, IC_YMAX, IC_YMIN,
!                                ! IBM_PROP, OPT_XMAX, OPT_XMIN,
!                                ! OPT_YMAX, OPT_YMIN, FILTER_R
!
!!     argument list
!      integer fid               ! file id
!
!!     local variables
!      integer ierr
!
!!     namelists
!      namelist /LEVELSET/ AMPMSK, OPT_XMAX, OPT_YMAX, OPT_YMIN,
!     $     FILTER_R, IBM_PROP, KERNEL_, QF_ramp, QA_ramp, QC_ramp
!     $     , FLTR_type, FLTR_beta, P_SIMP, MAP_RHO
!!-----------------------------------------------------------------------
!      ierr=0
!      if (NID.eq.0) then
!         write(unit=fid,nml=LEVELSET,iostat=ierr)
!      endif
!      call err_chk(ierr,'Error writing LEVELSET parameters.$')
!
!      return
!      end
!***********************************************************************
!     set the conjugate heat transfer properties
      subroutine IBM_properties(utrans,udiff,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'            !
      include 'PARALLEL'        ! GLLEL
      include 'TSTEP'           ! IFIELD, ISTEP
      include 'LEVELSET'        ! IBM_PROP, MSKLS
      include 'ADJOINT'

!     argument list
      integer ix,iy,iz,ieg
      real utrans,udiff

!     local variables
      integer iel
      real maskf_v,qf,qc,cl,cu,kl,ku,p_
!-----------------------------------------------------------------------
!     Set the real properties according to the mask distribution
!     Set fake properties for all fields in order to avoid any
!     complain about diff<=0
      if (ISTEP.lt.1) then
         utrans=1.0
         udiff=1.0
      else
         iel = GLLEL(ieg)
         maskf_v = IBM_MSKLS(ix,iy,iz,iel)
         if (IBM_MAP_RHO.eq.1) then
            qc=IBM_QC_ramp
            qf=IBM_QF_ramp
            if (IFIELD.eq.1) then
!     Flow properties
               utrans=IBM_PROP(1,1,2)
               udiff =IBM_PROP(1,1,1)
            elseif (IFIELD.ge.2) then
!     Thermal properties
               cl     = IBM_PROP(1,2,2)
               cu     = IBM_PROP(2,2,2)
               utrans = cl + (cu-cl)*((1-maskf_v)/(1+qc*maskf_v)) 
               IBM_CLITRANS(ix,iy,iz,iel)
     $              = cl - cl*((1-maskf_v)/(1+qc*maskf_v))
               
               kl     = IBM_PROP(1,2,1)
               ku     = IBM_PROP(2,2,1)
               udiff  = kl + (ku-kl)*((1-maskf_v)/(1+qf*maskf_v))
            endif
         elseif(IBM_MAP_RHO.eq.2)then
            p_=IBM_P_SIMP
            if (IFIELD.eq.1) then
!     Flow properties
               utrans=IBM_PROP(1,1,2)
               udiff =IBM_PROP(1,1,1)
            elseif (IFIELD.ge.2) then
!     Thermal properties
               cl     = IBM_PROP(1,2,2)
               cu     = IBM_PROP(2,2,2)
               utrans = cl + (cu-cl)*(1.0-maskf_v)**p_ 
               IBM_CLITRANS(ix,iy,iz,iel)=cl-cl*(1-maskf_v)**p_
               
               kl     = IBM_PROP(1,2,1)
               ku     = IBM_PROP(2,2,1)
               udiff  = kl + (ku-kl)*(1-maskf_v)**p_
            endif
         endif
      endif
      
      return
      end
!***********************************************************************
!     calcualte immersed boundary forcing
      subroutine IBM_forcing(ffx,ffy,ffz,ix,iy,iz,ieg,ux,uy,uz)
      implicit none

      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! VTRANS
!      include 'ADJOINT'         ! IFADJ
      include 'LEVELSET'        ! AMPMSK, MSKLS
      
!     argument list
      real ffx, ffy, ffz, ux, uy, uz
      integer ip,ix,iy,iz,ieg

!     local variables
      integer iel
      real maskf_v, alpha,qa,p_
!-----------------------------------------------------------------------
!     Initialization of the variables ff[xyz] is done in userf
      iel     = GLLEL(ieg)   	
      maskf_v = IBM_MSKLS(ix,iy,iz,iel)
      if(IBM_MAP_RHO.eq.1)then
         qa=IBM_QA_ramp
         alpha   = IBM_AMPMSK * (1.0 - maskf_v)/(1.0 + qa*maskf_v)
      elseif(IBM_MAP_RHO.eq.2)then
         p_=IBM_P_SIMP
         alpha   = IBM_AMPMSK * (1.0 - maskf_v)**p_
      endif         
!      if (.not.IFADJ) then
!     Direct
!     1) force velocity to zero just in the solid (multiply the IBM
!     force with the mask)
         ffx  = ffx - (ux*alpha)/VTRANS(ix,iy,iz,iel,1)
         ffy  = ffy - (uy*alpha)/VTRANS(ix,iy,iz,iel,1)
         if(IF3D) ffz  =   ffz -
     $        (uz*alpha)/VTRANS(ix,iy,iz,iel,1)
!      else
!     Adjoint
!     1) coupling term (gradT)Tadj from the advection of the direct
!     (multiply the ff[xyz] coming from cht_forcing with VTRANS(2)
!     
!     2) force velocity to zero just in the solid (IBM force is self
!     adjoint)
!         ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(iel-1))) 
!         ffx  = ffx*CLITRANS(ix,iy,iz,iel)-
!     $        (VXP(ip,1)*alpha)/VTRANS(ix,iy,iz,iel,1)
!         ffy  = ffy*CLITRANS(ix,iy,iz,iel)-
!     $        (VYP(ip,1)*alpha)/VTRANS(ix,iy,iz,iel,1)
!         if(IF3D) ffz  = ffz*CLITRANS(ix,iy,iz,iel) -
!     $        (VZP(ip,1)*alpha)/VTRANS(ix,iy,iz,iel,1)
!      endif
         
      return
      end
!***********************************************************************
!     impose the boundary conditions
      subroutine IBM_bc(x,y,z,ifield,temp)
      implicit none

      include 'SIZE'            
!      include 'ADJOINT'         ! IFADJ         
      include 'LEVELSET'        ! XBC_MAX
      include 'CONHT'           ! CHRA
      real x, y, z, temp
      integer ifield
!-----------------------------------------------------------------------
!     Dirichlet boundary conditions on the temperature
!	Likely not the same as me!
!      if(.not.IFADJ) then       !non-linear (direct)
         temp=0.0 
         
!      else                      !linear (adjoint)
!         temp=0.0
!         if(x.gt.0.4990)temp = -sqrt(CHRA)
!      endif
      
      return
      end
!***********************************************************************
!     define the initial conditions
      subroutine IBM_ic(x,y,z,ifield,temp,ux,uy,uz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D
      include 'PARALLEL'
      include 'LEVELSET'        ! PHIMIN, PHIMAX, IC_XMAX, IC_YMAX, IC_YMIN
                                ! OPT_XMAX, OPT_XMIN, OPT_YMAX, OPT_YMIN

!     argument list
      real x, y, z, temp, ux, uy, uz, dist,a_,n_
      integer ifield, ix, iy, iz, ieg, iel

!     local variables
!-----------------------------------------------------------------------
!	Likely not the same as me!
      if (ifield.eq.1) then
!     Velocity
		iel=GLLEL(ieg)
		if(IBM_MSKNF(ix,iy,iz,iel).eq.1) then
         		ux= 4.0*y*(1. - y)
		else
			ux=0
		endif
		uy=0.0 
      elseif (ifield.eq.2) then
!     Temperature
         temp=0      
      elseif (ifield.eq.3) then
!     Material distribution
         iel=GLLEL(ieg)
         temp=IBM_MSKNF(ix,iy,iz,iel)
      endif
      return
      end

!***********************************************************************
!     Clear the mask!
      subroutine IBM_clear_mask()
      implicit none

      include 'SIZE'
      include 'INPUT'          
      include 'PARALLEL'
      include 'LEVELSET'        

!     local variables
      integer ix, iy, iz, iel
!-----------------------------------------------------------------------
      do iel=1,NELT            
         do iz=1,NZ1
            do iy=1,NY1
               do ix=1,NX1
			IBM_MSKNF(ix,iy,iz,iel) = 1
			IBM_MSKLS(ix,iy,iz,iel) = 1
               enddo
            enddo
         enddo
      enddo
      return
      end

!***********************************************************************
!     Clear the mask!
      subroutine IBM_testing_mask()
      implicit none

      include 'SIZE'
      include 'INPUT'          
      include 'PARALLEL'
      include 'LEVELSET'
	include 'GEOM'        

!     local variables
      integer ix, iy, iz, iel
!-----------------------------------------------------------------------
      do iel=1,NELT            
         do iz=1,NZ1
            do iy=1,NY1
               do ix=1,NX1
			if(ym1(ix,iy,iz,iel).gt.1.or.ym1(ix,iy,iz,iel).lt.0) then
				IBM_MSKNF(ix,iy,iz,iel) = 0
				!IBM_MSKLS(ix,iy,iz,iel) = 0
			else
				IBM_MSKNF(ix,iy,iz,iel) = 1
				!IBM_MSKLS(ix,iy,iz,iel) = 1
			endif
               enddo
            enddo
         enddo
      enddo
      return
      end
!!***********************************************************************
!!     prepare the material distribution before the next iteration:
!!     if no restart, project the solution between 0-1, if restart,
!!     just copy the field in the commonblock
!      subroutine IBM_definemask
!      implicit none
!
!      include 'SIZE'            ! LX1, LY1, LZ1, LELT
!      include 'SOLN'            ! T
!      include 'INPUT'           ! NPASCAL
!      include 'ADJOINT'         ! DTD[XYZ]
!      include 'CHKPOINT'        ! IFCHKPTRST
!      include 'LEVELSET'        ! MSKLS, MSKNF
!
!!     local variables
!      integer n
!!-----------------------------------------------------------------------
!      n=NELT*NX1*NY1*NZ1
!      
!      if(.not.IFCHKPTRST) then  ! norestart   
!!     Cut values below zero and above one
!         call IBM_setlsmask
!         call copy(MSKNF,T(1,1,1,1,2),n)
!      else                      ! restart
!!     If restart copy the mask into its commonblock
!         call copy (MSKNF,T(1,1,1,1,2),n)
!!     Compute the temperature gradient
!         call gradm1(DTDX,DTDY,DTDZ,T(1,1,1,1,1))
!         NPSCAL = 1
!      endif
!
!      return
!      end
!!***********************************************************************
!!     cut values below zero and above one
!      subroutine IBM_setlsmask
!      implicit none
!
!      include 'SIZE'            ! NX1, NY1, NZ1, NELT
!      INCLUDE 'SOLN'            ! T
!
!!     local variables
!      integer iel, k, j, i,n
!!-----------------------------------------------------------------------
!!     Clio; initialization of the variables ff? to 0.0 is done in userf, 
!!     before calling cht_forcing and IBM_forcing
!      n=NX1*NY1*NZ1*NELT
!      do iel=1,NELT
!         do k=1,NZ1
!            do j=1,NY1
!               do i=1,NX1
!                  if (T(i,j,k,iel,2).le.0.0) then
!                     T(i,j,k,iel,2) = 0.0
!                  elseif (T(i,j,k,iel,2).ge.1.0) then
!                     T(i,j,k,iel,2) = 1.0
!                  endif
!               enddo
!            enddo
!         enddo
!      enddo
!
!      return
!      end
!!***********************************************************************
!!     Overrite the commonblock of the passive scalars for restart
!!     purposes and call restart
!      subroutine IBM_callchekpoint
!      implicit none
!
!      include 'SIZE'            
!      include 'TSTEP'           ! ISTEP
!      include 'SOLN'            ! T
!      include 'INPUT'           ! NPSCAL
!      include 'LEVELSET'        ! MSKLS, MSKNF
!
!      integer n
!!-----------------------------------------------------------------------
!      NPSCAL=3
!      n=NX1*NY1*NZ1*NELT
!      if (ISTEP.gt.0) then
!         call copy (T(1,1,1,1,3),MSKLS,n)
!         call copy (T(1,1,1,1,2),MSKNF,n)
!      endif
!      call checkpoint
!      NPSCAL=1
!      
!      return
!      end
!!***********************************************************************
!!     Set zero velocity as initial condition inside the solid
!      subroutine IBM_zeroVsolADJ
!      implicit none
!
!      include 'SIZE'            
!      include 'SOLN'            ! VXP,VYP,YZP
!      include 'INPUT'           ! IF3D
!      include 'LEVELSET'        ! MSKLS
!
!!     local variables
!      integer n
!!-----------------------------------------------------------------------
!      n=NX1*NY1*NZ1*NELT
!      call col2(VXP,MSKLS,n)  
!      call col2(VYP,MSKLS,n)
!      if (IF3D) call col2(VZP,MSKLS,n)
!      return
!      end
!!***********************************************************************
