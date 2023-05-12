!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision, parameter :: M_sun = 1.9891000d33
  double precision, parameter :: R_sun = 6.9599000d10
  double precision, parameter :: L_sun = 3.8268000d33
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: StefBoltz

  double precision :: cak_Q, cak_a, cak_base, cak_x0, cak_x1
  integer :: it_start_cak
  double precision :: rho_bound, v_inf, Mdot, static_rph
  double precision :: T_bound, R_lowerb, M_star, R_star

  double precision :: kappa_e, L_bound, Gamma_e_bound, F_bound, gradE, E_out
  logical :: fixed_lum, Cak_in_D, fd_correction, kappa_read_table
  logical :: read_cak_table = .true.
  logical :: CAK_zero = .false.
  logical :: hydro_lowerdensity = .false.
  
  integer :: i_v1, i_v2, i_p
  integer :: i_Trad, i_Tgas, i_Opal, i_CAK2, i_fdfac, i_beta_op, i_lambda,  i_gradE 
  !integer :: i_CAK, i_Mdot
  integer :: i_Gamma, i_Lum, i_F1, i_F2, i_edd

  double precision :: sum_time, sum_radius, iter
  double precision :: R_photosphere, R_photosphere_var
  double precision, allocatable :: vr_sumt(:), rho_sumt(:), rho2_sumt(:), vr2_sumt(:), rad_e_sumt(:), sumt(:)
  double precision, allocatable :: rhovr_sumt(:), rho2vr_sumt(:)

  double precision :: Omega = 0.d0
  double precision :: v_rot = 0.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    !> lasy fix for inexplicable pressure
    usr_internal_bc => Fix_pressure

    ! PseudoPlanar correction
    usr_source => PseudoPlanar

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => OPAL_and_CAK

    !> Set maximum value for diffusion coefficient
    usr_special_diffcoef => ceil_diffcoef

    ! Write out energy levels and temperature
    usr_write_analysis => collapse_to_1D

    !> Additional variables
    usr_process_grid => update_extravars

    ! Output routines
    ! usr_aux_output    => specialvar_output
    ! usr_add_aux_names => specialvarnames_output

    ! Timestep for PseudoPlanar
    ! usr_get_dt => get_dt_cak

    ! Refine mesh near base
    usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate()

    i_v1 = var_set_extravar("v1", "v1")
    i_v2 = var_set_extravar("v2", "v2")
    if (rhd_energy) i_p = var_set_extravar("p","p")
    i_Trad = var_set_extravar("Trad", "Trad")
    if (rhd_energy) i_Tgas = var_set_extravar("Tgas", "Tgas")
    i_Opal = var_set_extravar("OPAL", "OPAL")
    i_CAK2 = var_set_extravar("CAK2", "CAK2")
    i_Edd = var_set_extravar("Edd", "Edd")
    i_Gamma = var_set_extravar("Gamma", "Gamma")
    i_Lum = var_set_extravar("Lum", "Lum")
    i_F1 = var_set_extravar("F1", "F1")
    i_F2 = var_set_extravar("F2", "F2")
    i_fdfac = var_set_extravar("fdfac", "fdfac")
    !i_Mdot = var_set_extravar("Mdot", "Mdot")
    !i_CAK = var_set_extravar("CAK", "CAK")  
    i_lambda = var_set_extravar("lambda", "lambda")    
    i_gradE = var_set_extravar("gradE", "gradE")
    i_beta_op = var_set_extravar("beta_op", "beta_op")
    
  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_opal_opacity, only: init_opal
    use mod_cak_opacity, only: init_cak

    use mod_fld

    integer :: i

    !> Initialise Opal
    call init_opal(He_abundance,'Y02800')

    call init_cak('Y09800')

    !> read usr par
    call params_read(par_files)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_lowerb
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*const_mp)
    unit_velocity = v_inf

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    R_lowerb = R_lowerb/unit_length
    R_star = R_star/unit_length
    M_star = M_star/unit_density/unit_length**3
    !T_bound = T_bound/unit_temperature
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    v_rot = v_rot/unit_velocity
    Mdot = Mdot/unit_density/unit_length**3*unit_time

    kappa_e = kappa_e/unit_opacity
    F_bound = F_bound/unit_radflux
    L_bound = L_bound/(unit_radflux*unit_length**2)

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)

    !> Very bad initial guess for gradE using kappa_e
    gradE = -F_bound*3*kappa_e*rho_bound*unit_velocity/const_c

    if (mype ==0) then
     print*, 'L_bound (cgs)', L_bound*(unit_radflux*unit_length**2)
     print*, 'log10(L_bound)', log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
     print*, 'L_bound', L_bound*(unit_radflux*unit_length**2)/L_sun
     print*, ' '
     print*, 'unit_density', unit_density
     print*, 'unit_time', unit_time
     print*, 'unit_pressure', unit_pressure
     print*, 'unit_length', unit_length
     print*, 'unit_numberdensity', unit_numberdensity
     print*, 'unit_radflux', unit_radflux
     print*, 'unit_temperature', unit_temperature
     print*, 'unit_opacity', unit_opacity
     print*, 'unit_velocity', unit_velocity
   endif

    sum_time = 0.d0
    sum_radius= 0.d0
    iter = 0.d0

    allocate(vr_sumt(domain_nx1))
    vr_sumt = 0.d0
    allocate(rho_sumt(domain_nx1))
    rho_sumt = 0.d0
    allocate(rho2_sumt(domain_nx1))
    rho2_sumt = 0.d0
    allocate(vr2_sumt(domain_nx1))
    vr2_sumt = 0.d0
    allocate(rhovr_sumt(domain_nx1))
    rhovr_sumt = 0.d0
    allocate(rho2vr_sumt(domain_nx1))
    rho2vr_sumt = 0.d0
    allocate(rad_e_sumt(domain_nx1))
    rad_e_sumt = 0.d0

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound, kappa_e, &
    T_bound, R_lowerb, R_star, M_star, v_inf, Mdot, Gamma_e_bound, it_start_cak, fixed_lum, &
    Cak_in_D, read_cak_table, CAK_zero, Omega, fd_correction, kappa_read_table, static_rph

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    R_lowerb = R_lowerb*R_sun
    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_bound = Gamma_e_bound * 4.0 * dpi * const_G * M_star * const_c/kappa_e
    F_bound = L_bound/(4*dpi*R_lowerb**2)
    Mdot = Mdot*M_sun/year
    
    v_rot = Omega*dsqrt(const_G*M_star/R_lowerb)
    ! print*, Gamma_e_bound, dpi, const_G, M_star, const_c

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_constants
    use mod_physics

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: vel(ixI^S)
    double precision :: E_gauge
    double precision :: T_in, E_in, rr(ixI^S), bb, temp(ixI^S)

    integer :: ii

    do ii = ixOmin1,ixOmax1
      w(ii,:,rho_) = read_initial_conditions(x(ii,3,1),2)/unit_density
      w(ii,:,mom(1)) = read_initial_conditions(x(ii,3,1),3)/(unit_velocity)
      w(ii,:,mom(2)) = (v_rot + 1.d-1*dsin(10*(x(ii,:,1)-1.d0))*dcos(6*dpi*x(ii,:,2)))
      if (rhd_energy) w(ii,:,e_) = read_initial_conditions(x(ii,3,1),4)/unit_pressure
      w(ii,:,r_e) = read_initial_conditions(x(ii,3,1),5)/unit_pressure
      !> calculating diff coef from the initial values
      w(ii,:,i_diff_mg) = (const_c/unit_velocity)*(1./3.)/(kappa_e/unit_opacity*w(ii,:,rho_))
    enddo
    

    call phys_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions



  function read_initial_conditions(r_in,index) result(var)
    integer, intent(in) :: index
    double precision, intent(in) :: r_in
    double precision :: var

    double precision :: w(1:6), w_mo(1:6), w_po(1:6)
    integer :: ll

    w(:) = 0.d0

    open(unit=1, file='1D_stable.txt')
    read(1,*) w !> first line of data
    do ll = 1,2700
      w_mo = w
      read(1,*) w
      if (w(1) .gt. r_in*unit_length) then
          w_po = w
          goto 8765
       endif
      w_po = w
    enddo


8765 CLOSE(1)
    !print*,ll, r_in*unit_length, w_mo(1), w_po(1), w_mo(2)
    !var = w_mo(index) + (w_po(index) - w_mo(index))/(w_po(1) - w_mo(1))*(r_in - w_mo(1))
     var = w_mo(index)

    if ((var .lt. 1.d-12*unit_density) .and. (index .eq. 2)) then
      print*, '-------------------------------------'

      print*, ll, r_in*unit_length, var
      print*, w_mo
      print*, w_po

      stop

    end if


  end function read_initial_conditions


  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixI^S), Temp(ixI^S)
    double precision :: Temp0, rho0, T_out, n
    double precision :: Local_gradE(ixI^S), F_adv(ixB^S)
    double precision :: Local_tauout(ixI^S)
    double precision :: Local_Tout(ixI^S), Local_Eout(ixI^S)

    double precision :: kappa_out, kappa_in, del_r, kappa_dummy
    double precision :: radflux_out, diff_out
    double precision :: x_outer
    
    double precision :: csr, csg, cst
    
    double precision :: rad_flux(ixI^S,1:ndim)

    integer :: ix^D

    select case (iB)

    case(1)
    !> Density in the GCs assumed by extrapolating
    !> by assuming hydrostatic eqv
      if (hydro_lowerdensity) then 
        do ix1 = ixBmax1-1,ixBmin1,-1
          w(ix1,ixBmin2:ixBmax2,rho_) = dexp(2*dlog(w(ix1+1,ixBmin2:ixBmax2,rho_))&
          - dlog(w(ix1+2,ixBmin2:ixBmax2,rho_)))
        enddo
      else
        w(ixB^S,rho_) = rho_bound
      end if  
            
      !> rho_o*vo_r = rho*v_r, Mass conservation -- for v_r (mom(1))
      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,mom(1)) = w(ix1+1,ixBmin2:ixBmax2,mom(1))&
        *x(ix1+1,ixBmin2:ixBmax2,1)**2/x(ix1,ixBmin2:ixBmax2,1)**2        
        w(ix1,ixBmin2:ixBmax2,mom(2)) = w(ix1+1,ixBmin2:ixBmax2,mom(2))
      enddo
      
      !print*, w(ixB^S, mom(1))
      !print*, rho_bound
      !print*, 1.d-6*rho_bound
      
      !>Stellar rotation
      w(ixB^S,mom(2)) = v_rot*rho_bound
      
      !> Gas sound speed
      !> cs_gas^2 =  gamma*P_gas/rho
      !> Pgas = e_/(gamma-1) 
      csg = sum(w(ixImin1+nghostcells,ixImin2+nghostcells:ixImax2-nghostcells, e_) &
      - w(ixImin1+nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,mom(1))**2/&
      (2.*w(ixImin1+nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,rho_)))&
      /(ixImax2-ixImin2-2*nghostcells)
      csg = rhd_gamma*csg/(rhd_gamma -1.) 
      
      !> Radiation sound speed
      !> cs_rad^2 = 4*r_e/(9*rho)
      csr = sum(w(ixImin1+nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,r_e) &
      /w(ixImin1+nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,rho_))&
      /(ixImax2-ixImin2-2*nghostcells)
      csr= csr*(4./9.)
      
      !> Total sound speed    
      cst = sqrt(csg + csr)
      
      !> The updated floor values are provided by
      !> the total sound speed
      
      !> Floor value for negative inflow
      where(w(ixB^S,mom(1)) .lt. -cst*rho_bound)
       w(ixB^S,mom(1)) = -cst*rho_bound
      endwhere

      !> Ceil value for positive outflow
      where(w(ixB^S,mom(1)) .gt. cst*rho_bound)
       w(ixB^S,mom(1)) = cst*rho_bound
      endwhere

      !> Floor value for negative inflow
      !where(w(ixB^S,mom(1)) .lt. -1.d-6*rho_bound)
      ! w(ixB^S,mom(1)) = -1.d-6*rho_bound
      !endwhere

      !> Ceil value for positive outflow
      !where(w(ixB^S,mom(1)) .gt. 1.d-6*rho_bound)
      ! w(ixB^S,mom(1)) = 1.d-6*rho_bound
      !endwhere
      

      !F_adv(ixBmax1,ixBmin2:ixBmax2) = 4.d0/3.d0*(w(ixBmax1,ixBmin2:ixBmax2,mom(1))/w(ixBmax1,ixBmin2:ixBmax2,rho_))*w(ixBmax1,ixBmin2:ixBmax2,r_e) &
      !           * 4*dpi*xprobmin1**2

      !where (F_adv(ixB^S) .ne. F_adv(ixB^S)) F_adv(ixB^S) = 0.d0
      !where (F_adv(ixB^S) .le. 0.d0) F_adv(ixB^S) = 0.d0

      !> Calculate gradE using the FLD closure, impose rational floor value on gradE
      !do ix1 = ixImin1,ixImax1
      !  do ix2 = ixBmin2,ixBmax2
      !    Local_gradE(ix1,ix2) = -(F_bound-F_adv(ixBmax1,ix2))/w(nghostcells+1,ix2,i_diff_mg)
      !    Local_gradE(ix1,ix2) = max(Local_gradE(ix1,ix2),-200.d0)
      !  enddo
      !enddo
      !gradE = sum(Local_gradE(nghostcells,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)


      !> Extrapolate using gradE, but impose some rational ceil value for Erad near boundary
      !do ix1 = ixBmax1,ixBmin1,-1
      !  do ix2 = ixBmin2,ixBmax2
      !    w(ix1,ix2,r_e) = min(2.d0,w(ix1+2,ix2,r_e)) &
      !    + (x(ix1,ix2,1)-x(ix1+2,ix2,1))*Local_gradE(ix1+1,ix2)
      !  enddo
      !enddo
      
      call get_kappa_OPAL(ixI^L,ixI^L,w,x,kappa)

      !kappa_in = 0.34/unit_opacity
      kappa_in = sum(kappa(ixImin1-nghostcells,ixBmin2:ixBmax2))/((ixBmax2-ixBmin2))

      if (kappa_in .ne. kappa_in) kappa_in = kappa_e
      kappa_in = max(kappa_in,kappa_e)
      kappa_in = min(kappa_in,20*kappa_e)
      
      !> updating gradE with kappa_in
      gradE = -F_bound*3*kappa_in*rho_bound*unit_velocity/const_c    
      
      if (fixed_lum) then
        do ix1 = ixBmax1,ixBmin1,-1
          do ix2 = ixBmin2,ixBmax2
            w(ix1,ix2,r_e) =  w(ix1+1,ix2,r_e) - w(ix1,ix2,rho_)*kappa_in*F_bound*&
            (x(ix1+1,ix2,1)-x(ix1,ix2,1))*unit_velocity/const_c
          enddo
        enddo    
      else
        w(ixB^S,r_e) = (const_rad_a*T_bound**4)/unit_pressure
      end if 
      
      temp(ixB^S) = max(5.d0*T_bound/unit_temperature,&
      (w(ixB^S,r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature)       
       
      if (rhd_energy) then
       w(ixB^S,e_) = w(ixB^S,rho_)*unit_density*temp(ixB^S)*unit_temperature*const_kB/(fld_mu*const_mp)/&
       (rhd_gamma-1.d0)/unit_pressure + half*(w(ixB^S,mom(1))**2+w(ixB^S,mom(2))**2)/w(ixB^S,rho_)
      end if        

   case(2)
      
      if (fixed_lum) then
        call fld_get_radflux(w, x, ixI^L, ixI^L, rad_flux)
        diff_out = sum(w(ixImax1-nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,i_diff_mg))&
        /(ixImax2-ixImin2)
        radflux_out = sum(rad_flux(ixImax1-nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,1))&
        /(ixImax2-ixImin2)/F_bound
        x_outer =  sum(x(ixImax1-nghostcells,ixImin2+nghostcells:ixImax2-nghostcells,1))/&
        (ixImax2-nghostcells-ixImin2-nghostcells + 1.)
        radflux_out =max(1.d-1*F_bound/x_outer**2,radflux_out)
        !Local_Eout(ixB^S) = radflux_out*x_outer**4/(3*diff_out*x(ixB^S,1)**3)
        !E_out = sum(Local_Eout(ixImax1-nghostcells,ixImin2+nghostcells:ixImax2-nghostcells))&
        !/(ixImax2-ixImin2)
        E_out = radflux_out*x_outer/(3*diff_out)
        E_out = max(0.d0, E_out)
        T_out = (E_out/const_rad_a)**0.25/unit_temperature
        T_out = max(5.d2/unit_temperature, T_out)
        E_out = const_rad_a*(T_out*unit_temperature)**4.0/unit_pressure
        w(ixB^S,r_e) = E_out
     
      else
        call get_kappa_OPAL(ixI^L,ixI^L,w,x,kappa)
        kappa_out = sum(kappa(ixImax1-nghostcells,ixBmin2:ixBmax2))/((ixBmax2-ixBmin2))
        if (kappa_out .ne. kappa_out) kappa_out = kappa_e
          kappa_out = max(kappa_out,kappa_e)
          kappa_out = min(kappa_out,20*kappa_e)
          Local_tauout(ixB^S) = kappa_out*w(ixB^S,rho_)*R_lowerb**2/(3*x(ixB^S,1))
          Local_Tout(ixB^S) = (F_bound/StefBoltz*(1.d0/4.d0*Local_tauout(ixB^S)))**0.25d0
          T_out = sum(Local_Tout(ixBmin1,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)
          T_out = max(5.d2/unit_temperature, T_out)
          E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure
          w(ixB^S,r_e) = const_rad_a*(Local_Tout(ixB^S)*unit_temperature)**4.d0/unit_pressure
 
      end if 
      
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine mg_boundary_conditions(iB)
    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB
    integer          :: iigrid, igrid
    double precision :: snd_g, rcv_g, snd_e, rcv_e

    select case (iB)
    case (1)
      !> communicate boundary
      ! do iigrid=1,igridstail; igrid=igrids(iigrid);
      !   if (node(pig1_,igrid) > 1) &
      !     gradE = 0.d0
      ! enddo
      !
      ! snd_g = gradE
      ! call mpi_allreduce(snd_g,rcv_g,1,&
      ! MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      ! gradE = rcv_g
      
      
      if (fixed_lum) then
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
        mg%bc(iB, mg_iphi)%bc_value = gradE
      else 
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = (const_rad_a*T_bound**4)/unit_pressure
      end if  

    case (2)
      !> communicate boundary
      ! do iigrid=1,igridstail; igrid=igrids(iigrid);
      !   if (node(pig1_,igrid) < domain_nx1/block_nx1*node(plevel_,igrid)) &
      !     E_out = bigdouble
      ! enddo
      !
      ! snd_e = E_out
      ! call mpi_allreduce(snd_e,rcv_e,1,&
      ! MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      ! E_out = rcv_e


      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = E_out

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions


  subroutine Fix_pressure(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: pth(ixI^S), mean_p

    !> fix density as well
    !where ((w(ixO^S,rho_) .lt. 1.d-5) .and. (x(ixO^S,1) .lt. 4.d0))
    !  w(ixO^S,rho_) = 1.d-5
    !endwhere

    !> Fix velocity in atmosphere
    where ((x(ixO^S,1) .lt. 1.2d0) .and. (w(ixO^S,mom(1))/w(ixO^S,rho_) .gt. 1.d-4))
      w(ixO^S,mom(1)) = 1.d-4*w(ixO^S,rho_)
    end where


    !if (.not. rhd_energy) &
    !  call mpistop("no energy equation, no pressure fix")

    !call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)

    !if (any(press(ixO^S) .lt. 0.d0)) then
    !  mean_p = max(sum(pth(ixO^S))/(block_nx1*block_nx2),small_pressure)
    !  where (pth(ixO^S) .le. small_pressure)
    !    w(ixO^S,e_) = mean_p/(rhd_gamma - 1) +  0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1)/w(ixO^S,rho_)
    !  end where
    !endif

    !!> Temperature ceil, Tmax 1d6
    !where ((pth(ixO^S)/w(ixO^S,rho_)*unit_temperature) .gt. 1.d6)
    !  w(ixO^S,e_) = 1.d6/unit_temperature*w(ixO^S,rho_)/(rhd_gamma - 1) &
    !      + 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1)/w(ixO^S,rho_)
    !end where

    !where (w(ixO^S,r_e) .gt. 0.5d0)
    !    w(ixO^S,r_e) = 0.5d0
    !end where


  end subroutine Fix_pressure


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,1:ndim)

    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixI^S,:) = 0.d0
    gravity_field(ixI^S,1) = -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ppsource(ixO^S,1:nw)

    double precision :: k_cak(ixO^S), rad_flux(ixO^S,1:ndim)

    call PseudoPlanarSource(ixI^L,ixO^L,wCT,x,ppsource)
    w(ixO^S,rho_) = w(ixO^S,rho_) + qdt*ppsource(ixO^S,rho_) !> OK
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*ppsource(ixO^S,mom(1)) !> OK
    w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt*ppsource(ixO^S,mom(2)) !> OK
    if (rhd_energy) w(ixO^S,e_) = w(ixO^S,e_) + qdt*ppsource(ixO^S,e_) !> OK
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e) !> TROUBLEMAKER

    if (.not. Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixI^L,ixO^L,wCT,x,k_cak)
      else
        call get_kappa_CAK(ixI^L,ixO^L,wCT,x,k_cak)
      endif

      if (fixed_lum) then
        !> Fixed L = L_bound
        w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
          + qdt*wCT(ixO^S,rho_)*L_bound/(4*dpi*x(ixO^S,1)**2)/const_c*k_cak(ixO^S)*unit_velocity
        if (rhd_energy) then
          w(ixO^S,e_) = w(ixO^S,e_) &
            + qdt*wCT(ixO^S,mom(1))*L_bound/(4*dpi*x(ixO^S,1)**2)/const_c*k_cak(ixO^S)*unit_velocity
        endif
      else
        !> Local flux
        call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux)
        w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
          + qdt*wCT(ixO^S,rho_)*rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity
        if (rhd_energy) then
          w(ixO^S,e_) = w(ixO^S,e_) &
            + qdt*wCT(ixO^S,mom(1))*rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity
        endif
      endif

    endif

  end subroutine PseudoPlanar

  subroutine get_dt_cak(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixI^S)
    double precision :: mass

    double precision :: dt_cak
    double precision :: k_cak(ixO^S), rad_flux(ixO^S,1:ndim)

    if (read_cak_table) then
      call get_kappa_CAK2(ixI^L,ixO^L,w,x,k_cak)
    else
      call get_kappa_CAK(ixI^L,ixO^L,w,x,k_cak)
    endif

    if (fixed_lum) then
      !> Fixed L = L_bound
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/(L_bound/(4*dpi*x(ixO^S,1)**2)&
      /const_c*k_cak(ixO^S)*unit_velocity -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length))))
    else
      !> Local flux
      call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/abs(rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity &
      -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length))))
    endif

    dtnew = min(dt_cak, dtnew)

  end subroutine get_dt_cak
  !

  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)

    double precision :: rad_flux(ixO^S,1:ndir)
    double precision :: pth(ixI^S),v(ixO^S,1:ndim)
    double precision :: radius(ixO^S),  pert(ixO^S)
    double precision :: edd(ixO^S,1:ndim,1:ndim)
    integer :: rdir, pdir

    source(ixO^S,1:nw) = zero

    rdir = 1
    pdir = 2

    v(ixO^S,rdir) = w(ixO^S,mom(rdir))/w(ixO^S,rho_)
    v(ixO^S,pdir) = w(ixO^S,mom(pdir))/w(ixO^S,rho_)

    radius(ixO^S) = x(ixO^S,rdir) ! + half*block%dx(ixO^S,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixO^S,rho_) = -two*w(ixO^S,rho_)*v(ixO^S,rdir)/radius(ixO^S)

    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixO^S,mom(rdir)) = - 2*w(ixO^S,rho_)*v(ixO^S,rdir)**two/radius(ixO^S) &
                              + w(ixO^S,rho_)*v(ixO^S,pdir)**two/radius(ixO^S)

    source(ixO^S,mom(pdir)) = - 3*v(ixO^S,rdir)*v(ixO^S,pdir)*w(ixO^S,rho_)/radius(ixO^S)

    !> de/dt = -2 (e + p)v_r/r
    if (rhd_energy) &
    source(ixO^S,e_) = -two*(w(ixO^S,e_)+pth(ixO^S))*v(ixO^S,rdir)/radius(ixO^S)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
      source(ixO^S,r_e) = source(ixO^S,r_e) - two*rad_flux(ixO^S,rdir)/radius(ixO^S)
    endif

    if (rhd_radiation_advection) then
      source(ixO^S,r_e) = source(ixO^S,r_e) - two*w(ixO^S,r_e)*v(ixO^S,rdir)/radius(ixO^S)
    endif

    ! Not sure about this one
    if (rhd_radiation_force) then
      call fld_get_eddington(w, x, ixI^L, ixO^L, edd)
      source(ixO^S,r_e) = source(ixO^S,r_e) + two*v(ixO^S,rdir)*w(ixO^S,r_e)*edd(ixO^S,1,1)/radius(ixO^S)
    endif

  end subroutine PseudoPlanarSource


  subroutine OPAL_and_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: OPAL(ixO^S), CAK(ixO^S)

    !> Get OPAL opacities by reading from table
    call get_kappa_OPAL(ixI^L,ixO^L,w,x,OPAL)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    if (Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixI^L,ixO^L,w,x,CAK)
      else
        call get_kappa_CAK(ixI^L,ixO^L,w,x,CAK)
      endif
    else
      CAK(ixO^S) = 0.d0
    endif

    !> Add OPAL and CAK for total opacity
    kappa(ixO^S) = OPAL(ixO^S) + CAK(ixO^S)

    where(kappa(ixO^S) .ne. kappa(ixO^S))
      kappa(ixO^S) = kappa_e
    endwhere

  end subroutine OPAL_and_CAK


  subroutine get_kappa_OPAL(ixI^L,ixO^L,w,x,kappa)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    integer :: ix^D
    double precision :: Temp(ixI^S)
    double precision :: n, rho0, Temp0

    !> Get OPAL opacities by reading from table
!    if (rhd_energy) then
!      call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
!    else
      call phys_get_trad(w,x,ixI^L,ixO^L,Temp)
!    endif

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,n)
        kappa(ix^D) = n/unit_opacity
    {enddo\ }
    
    where(kappa(ixO^S) .ne. kappa(ixO^S))
      kappa(ixO^S) = kappa_e
    endwhere

    !> Lower limit is electron scattering:
    where(kappa(ixO^S) .lt. kappa_e)
      kappa(ixO^S) = kappa_e
    endwhere


    !> test without opal
    ! kappa(ixO^S) = kappa_e

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    integer :: ix^D
    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S)
    double precision :: xx(ixO^S), alpha(ixO^S)

    if (.not. CAK_zero) then

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradientO(vel,x,ixI^L,ixO^L,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))

    xx(ixO^S) = 1.d0-xprobmin1/x(ixO^S,1)

    alpha(ixO^S) = cak_a

    where (xx(ixO^S) .le. cak_x0)
      alpha(ixO^S) = cak_base
    elsewhere (xx(ixO^S) .le. cak_x1)
      alpha(ixO^S) = cak_base + (cak_a - cak_base)&
      *(xx(ixO^S) - cak_x0)/(cak_x1 - cak_x0)
    endwhere

    kappa(ixO^S) = kappa_e*cak_Q/(1-alpha(ixO^S)) &
    *(gradv(ixO^S)*unit_velocity/(w(ixO^S,rho_)*const_c*cak_Q*kappa_e))**alpha(ixO^S)

    if (x(ixImax1,nghostcells,1) .ge. xprobmax1) then
      kappa(ixOmax1,ixOmin2:ixOmax2) = kappa(ixOmax1-1,ixOmin2:ixOmax2)
    endif

    ! if (it .le. it_start_cak) then
    !   kappa(ixO^S) = kappa(ixO^S)*dexp(-w(ixO^S,rho_)*kappa_e)
    ! endif

    !{do ix^D=ixOmin^D,ixOmax^D\ }
    !    if (xx(ix^D) .lt. cak_x0) then
    !      kappa(ix^D) = min(2*kappa_e,kappa(ix^D))
    !    else if (xx(ix^D) .lt. cak_x1) then
    !      kappa(ix^D) = min( (2 + (20.d0 - 2.d0)/(cak_x1 - cak_x0)*(xx(ix^D)-cak_x0) )*kappa_e, kappa(ix^D))
    !    else
    !      kappa(ix^D) = min(20*kappa_e,kappa(ix^D))
    !    endif
    !
    !  ! if (kappa(ix^D) .gt. 3*kappa_e) &
    !  !   kappa(ix^D) = 3*kappa_e + half*(kappa(ix^D)-3*kappa_e)
    !{enddo\ }

    {do ix^D=ixOmin^D,ixOmax^D\ }
      kappa(ix^D) = min(50*kappa_e,kappa(ix^D))
    {enddo\ }

    else

    !> test with no cak
    kappa(ixO^S) = 0.d0

    endif

  end subroutine get_kappa_CAK



subroutine get_kappa_CAK2(ixI^L,ixO^L,w,x,kappa)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_cak_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: Temp(ixI^S), radius(ixI^S), rho0, temp0, gradv0, kap0, vel0, rad0, fdfac0
    integer :: ix^D

    double precision :: alpha, Qbar, Q0, kappa_e_t, R_photo, kappa_e_table
    double precision :: tau, M_t
    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S), finite_disk(ixI^S), fdfac(ixO^S)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradientO(vel,x,ixI^L,ixO^L,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))
    radius(ixO^S) = x(ixO^S,1)
    
    !!radius = x(ixI^S)

    !> Get CAK opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
    else
      call phys_get_trad(w,x,ixI^L,ixO^L,Temp)
    endif
    
    !> fdfac is the finite disk correction factor
    !> the finite disk correction factor is implemented 
    !> to the CAK opacity
    !> below the photosphere of the star k_line is reduced by
    !> k_line/(1+ alpha) --DD
    
    call get_fdfac(ixI^L,ixO^L,w,x,fdfac) 
    finite_disk(ixO^S)= fdfac(ixO^S) 
       
    R_photo= R_photosphere/unit_length

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        vel0 = vel(ix^D)*unit_velocity  
        rad0 = radius(ix^D)
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix^D)*(unit_velocity/unit_length)
        fdfac0= finite_disk(ix^D)
        
        call set_cak_opacity(rho0,Temp0,gradv0,alpha, Qbar, Q0, kappa_e_t)
        kappa_e_table= kappa_e_t/unit_opacity
        
        if (kappa_read_table) then
          tau = (kappa_e_table*unit_opacity)*rho0*const_c/gradv0
        else
          tau = (kappa_e*unit_opacity)*rho0*const_c/gradv0
        end if   
                        
        M_t = Qbar/(1-alpha)*((1+Q0*tau)**(1-alpha) - 1)/(Q0*tau)
        
        if (kappa_read_table .and. fd_correction) then
          kap0 = (kappa_e_table*unit_opacity)*M_t*fdfac0
        else if (kappa_read_table) then
          kap0 = (kappa_e_table*unit_opacity)*M_t
        else if (fd_correction) then
          kap0 = (kappa_e*unit_opacity)*M_t*fdfac0
        else
          kap0 = (kappa_e*unit_opacity)*M_t
        end if  
        
        
        !if (fd_correction) then
        !  kap0 = (kappa_e*unit_opacity)*M_t*fdfac0 !> DD --finite disk correction factor
        !else
        !  kap0 = (kappa_e*unit_opacity)*M_t
        !end if   
        
        kappa(ix^D) = kap0/unit_opacity
        
        if (kappa(ix^D) .ne. kappa(ix^D)) kappa(ix^D) = 0.d0
        
        if (kappa_read_table) then
          kappa(ix^D) = min(50*kappa_e_table,kappa(ix^D))
        else
          kappa(ix^D) = min(50*kappa_e,kappa(ix^D))
        end if
        
    {enddo\ }


    if (x(ixImax1,nghostcells,1) .ge. xprobmax1) then
      kappa(ixOmax1,:) = kappa(ixOmax1-1,:)
    endif


  end subroutine get_kappa_CAK2

subroutine get_betaop(ixI^L,ixO^L,w,x,beta_op)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_cak_opacity
    use mod_fld
    !> calculating beta for the finite disk
    !> correction factor -- DD
    
    
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: beta_op(ixO^S) 

    double precision :: Temp(ixI^S), radius(ixI^S), rho0, temp0, gradv0, kap0, vel0, rad0
    integer :: ix^D

    double precision :: alpha, Qbar, Q0, kappa_e_t, beta_op0, R_photo
    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S)

    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradientO(vel,x,ixI^L,ixO^L,1,gradv,nghostcells)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))
    radius(ixO^S) = x(ixO^S,1)
    
    !> Get CAK opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
    else
      call phys_get_trad(w,x,ixI^L,ixO^L,Temp)
    endif

    R_photo= R_photosphere/unit_length

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        vel0 = vel(ix^D)*unit_velocity  
        rad0 = radius(ix^D) 
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix^D)*(unit_velocity/unit_length)
        
        call set_cak_opacity(rho0,Temp0,gradv0,alpha,Qbar,Q0,kappa_e_t)

        beta_op0 = (1.-vel0/(gradv0*rad0*unit_length))*(R_photo/rad0)**2        
        
        beta_op(ix^D) = beta_op0
        
    {enddo\ }
    
    if (x(ixImax1,nghostcells,1) .ge. xprobmax1) then
      beta_op(ixOmax1,:) = beta_op(ixOmax1-1,:)
    endif


  end subroutine get_betaop


subroutine get_fdfac(ixI^L,ixO^L,w,x,fdfac)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_cak_opacity
    use mod_fld
    
    !> fdfac is the 1D finite disk correction factor
    !> DD
    
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: fdfac(ixO^S) 

    double precision :: Temp(ixI^S), radius(ixI^S), rho0, temp0, gradv0, kap0, vel0, rad0
    integer :: ix^D

    double precision :: alpha, Qbar, Q0, kappa_e_t, oma, opa, fdfac0, beta_op, R_photo
    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradientO(vel,x,ixI^L,ixO^L,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))
    radius(ixO^S) = x(ixO^S,1)

    !> Get CAK opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
    else
      call phys_get_trad(w,x,ixI^L,ixO^L,Temp)
    endif

    R_photo= R_photosphere/unit_length

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        vel0 = vel(ix^D)*unit_velocity 
        rad0 = radius(ix^D)
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix^D)*(unit_velocity/unit_length)
        
        call set_cak_opacity(rho0,Temp0,gradv0,alpha,Qbar,Q0,kappa_e_t)
        
        opa = 1.+alpha
        oma = 1.-alpha
        
        beta_op = (1.-vel0/(gradv0*rad0*unit_length))*(R_photo/rad0)**2
        
        if (rad0.le.R_photo) then
         fdfac0 = 1./opa
        else if (rad0.gt.R_photo .and. beta_op.ge.1.0) then
         fdfac0 = 1./opa
        else if (rad0.gt.R_photo .and. beta_op.lt.-1.d10) then
         fdfac0 = ((-beta_op)**alpha)/opa
        else if (rad0.gt.R_photo .and. abs(beta_op).gt.1.d-3) then
         fdfac0 = (1.-(1.-beta_op)**opa)/(beta_op*opa)
        else
         fdfac0 = 1.-0.5*alpha*beta_op*(1.+0.3333333*oma*beta_op)
        end if    
                
        fdfac(ix^D) = fdfac0
        
    {enddo\ }
    

    if (x(ixImax1,nghostcells,1) .ge. xprobmax1) then
      fdfac(ixOmax1,:) = fdfac(ixOmax1-1,:)
    endif


  end subroutine get_fdfac



 subroutine ceil_diffcoef(w, wCT, x, ixI^L, ixO^L)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: wCT(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)


    where (w(ixO^S,i_diff_mg) .gt. 5.d4) &
      w(ixO^S,i_diff_mg) = 5.d4

    where (w(ixO^S,i_diff_mg) .lt. 1.d-10) &
      w(ixO^S,i_diff_mg) = 1.d-10


  end subroutine ceil_diffcoef


  subroutine refine_base(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision :: lim_1, lim_2, lim_3, lim_4

    lim_3 = 1.8d0
    lim_2 = 1.9d0
    lim_1 = 2.0d0

    !refine= -1
    !coarsen= -1

    if (all(x(ixG^S,1) < lim_3)) then
      if (level > 4) coarsen=1
      if (level < 4) refine=1
    elseif (all(x(ixG^S,1) < lim_2)) then
      if (level > 3) coarsen=1
      if (level < 3) refine=1
    elseif (all(x(ixG^S,1) < lim_1)) then
      if (level > 2) coarsen=1
      if (level < 2) refine=1
    endif

  end subroutine refine_base


  subroutine collapse_to_1D()
    use mod_global_parameters
    use mod_opal_opacity
    use mod_constants

    integer          :: iigrid, igrid, jj_blk,ii, jj, nbx, i, il, ih, dn_mdot
    integer          :: lvl, ibx
    integer          :: np_mdot, nc, n_int
    double precision :: ratio, sf_mdot, dx_l1, Rphoto0, Rph_dum, avg_photo_radius, prev_time, avg_radius
    double precision :: R_photosphere_dum, R_photosphere_mype, Rph_S, Rph_R
    double precision:: tau(1:domain_nx1), xx(1:domain_nx1), xx_S(1:domain_nx1), xx_R(1:domain_nx1)


    double precision :: radflux(ixM^T,ndim)

    double precision, allocatable :: rp_mdot(:), p_mdot(:)
    double precision, allocatable :: rp_lum(:),p_lum(:)
    double precision, allocatable :: mdot_S(:), mdot_R(:)
    double precision, allocatable :: lum_S(:), lum_R(:)
    integer, allocatable :: jp_mdot(:), number_data_points(:)

    double precision :: rr(1:domain_nx1), rr_S(1:domain_nx1), rr_R(1:domain_nx1)
    double precision :: vr(1:domain_nx1), vr_S(1:domain_nx1), vr_R(1:domain_nx1)
    double precision :: mdot(1:domain_nx1), lum(1:domain_nx1)
    double precision :: rho(1:domain_nx1), rho_S(1:domain_nx1), rho_R(1:domain_nx1)
    double precision :: rho2(1:domain_nx1), rho2_S(1:domain_nx1), rho2_R(1:domain_nx1)
    double precision :: vr2(1:domain_nx1), vr2_S(1:domain_nx1), vr2_R(1:domain_nx1)
    double precision :: rhovr(1:domain_nx1), rhovr_S(1:domain_nx1), rhovr_R(1:domain_nx1)
    double precision :: rho2vr(1:domain_nx1), rho2vr_S(1:domain_nx1), rho2vr_R(1:domain_nx1)
    double precision :: rad_e(1:domain_nx1),rad_e_S(1:domain_nx1),rad_e_R(1:domain_nx1), n, kappa(1:domain_nx1)
    double precision :: t_rad0(1:domain_nx1), kappa_0(1:domain_nx1), temp(1:domain_nx1)

    integer :: lvl_h(1:domain_nx1), lvl_h_S(1:domain_nx1), lvl_h_R(1:domain_nx1)
    integer :: lvl_l(1:domain_nx1), lvl_l_S(1:domain_nx1), lvl_l_R(1:domain_nx1)
        

    ! if (refine_max_level .ne. 1) &
    ! call mpistop("collapse_to_1D doesnt work YET with mpi")

    !> #R_lowerb -1 in simulation
    np_mdot = floor((xprobmax1-xprobmin1)/R_lowerb) + 1
    

    allocate(rp_mdot(1:np_mdot))
    allocate(p_mdot(1:np_mdot))
    allocate(p_lum(1:np_mdot))
    allocate(jp_mdot(1:np_mdot))
    allocate(mdot_S(1:np_mdot))
    allocate(mdot_R(1:np_mdot))
    allocate(lum_S(1:np_mdot))
    allocate(lum_R(1:np_mdot))


    rr = 0.d0
    rp_mdot = 0.d0
    vr = 0.0d0
    rho = 0.d0
    rad_e = 0.d0
    rho2 = 0.d0
    vr2 = 0.d0
    mdot = 0.d0
    p_mdot = 0.d0
    lum = 0.d0
    p_lum = 0.d0
    lvl_h = 0.d0
    lvl_l = 2* refine_max_level
    rhovr = 0.d0
    rho2vr = 0.d0
    xx= 0.d0


    !> Reconstruct radius at level 1
    do jj = 1,domain_nx1
      rr(jj) = xprobmin1 + (jj-0.5d0)*(xprobmax1-xprobmin1)/domain_nx1
    enddo

    !> Choose radii at which to save mdot
    !> This is done at every stellar radii
    do jj = 1,np_mdot
      dx_l1 = (xprobmax1-xprobmin1)/domain_nx1
      dn_mdot = floor(1.d0*domain_nx1/np_mdot+smalldouble)
      rp_mdot(jj) = xprobmin1*jj
    enddo
    rp_mdot(1) = xprobmin1 + dx_l1/2
    rp_mdot(np_mdot) = xprobmin1*np_mdot - dx_l1/2

    !> Find cells that are closest to radii where we want to track mdot
    do jj = 1,domain_nx1
      do ii = 1,np_mdot
        if (rr(jj) - rp_mdot(ii) .le. smalldouble) &
          jp_mdot(ii) = jj
      enddo
    enddo
    !print*, const_rad_a

    !> loop over all grids for this proc
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)

      !> widht in cells per cell of lvl 1: On lvl 1: 1, lvl 2: 2, lvl 3: 4
      nc = 2**(node(plevel_,igrid)-1)

      !> Map the block index to the global index
      !> nbx = nr of blocks between i=0 and current block
      jj_blk = (node(pig1_,igrid)-1)*block_nx1/nc

      if (nc > block_nx1) &
        call mpistop("collapse_to_1D doesnt work, reduce amr")

      !> Calculate radflux in block for luminosity
      call fld_get_radflux(block%w, block%x, ixG^LL, ixM^LL, radflux)

      !> For all cells in the current block, average velocity over lateral direction.
      !> Take into account possible amr
      do ii = 1,block_nx1/nc
        !> relate global index jj to local index ii on grid
        il = nghostcells + 1 + (ii-1)*nc
        ih = nghostcells + 1 + (ii)*nc -1
        jj = jj_blk + ii

        !> Velocity
        vr(jj) = vr(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,mom(1))/block%w(il:ih,ixMlo2:ixMhi2,rho_))/domain_nx2/nc**2 
        !> average value on lvl 1

        !> Density
        rho(jj) = rho(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,rho_)*unit_density)/domain_nx2/nc**2 
        
        !> radiation energy density
        rad_e(jj)= rad_e(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,r_e)*unit_pressure)/domain_nx2/nc**2 

        !> xx
        xx(jj) = xx(jj) &
        + sum(block%x(il:ih,ixMlo2:ixMhi2,1))/domain_nx2/nc**2
        
        !> Density squared
        rho2(jj) = rho2(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,rho_)**2)/domain_nx2/nc**2 
        !> radial velocity squared
        vr2(jj) = vr2(jj) &
        + sum((block%w(il:ih,ixMlo2:ixMhi2,mom(1))/block%w(il:ih,ixMlo2:ixMhi2,rho_))**2)/domain_nx2/nc**2

        !> Density weighted velocity
        rhovr(jj) = rhovr(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,mom(1)))/domain_nx2/nc**2

        !> Density squared weighted velocity
        rho2vr(jj) = rho2vr(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,mom(1))*block%w(il:ih,ixMlo2:ixMhi2,rho_))/domain_nx2/nc**2
        
        !> Mass loss rate
        mdot(jj) = mdot(jj) &
        + 4*dpi*sum(block%x(il:ih,ixMlo2:ixMhi2,1)**2*block%w(il:ih,ixMlo2:ixMhi2,mom(1)))/domain_nx2/nc**2 

        !> Luminosity
        lum(jj) = lum(jj) &
        + 4*dpi*sum(block%x(il:ih,ixMlo2:ixMhi2,1)**2*radflux(il:ih,ixMlo2:ixMhi2,1))/domain_nx2/nc**2 

        !> Highest/Lowest amr level
        lvl_h(jj) = node(plevel_,igrid)
        lvl_l(jj) = node(plevel_,igrid)
      enddo

    enddo
    
    

    !> communicate velocity
    vr_S=vr
    call mpi_allreduce(vr_S,vr_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    vr=vr_R

    !> communicate density
    rho_S=rho
    call mpi_allreduce(rho_S,rho_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rho=rho_R
    
    !> communicate t_rad
    rad_e_S=rad_e
    call mpi_allreduce(rad_e_S,rad_e_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rad_e=rad_e_R
    
    !> communicate radial coordinate -- xx
    xx_S= xx
    call mpi_allreduce(xx_S,xx_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    xx=xx_R

    !> communicate density squared
    rho2_S=rho2
    call mpi_allreduce(rho2_S,rho2_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rho2=rho2_R

    !> communicate vr squared
    vr2_S=vr2
    call mpi_allreduce(vr2_S,vr2_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    vr2=vr2_R

    !> communicate rho*vr
    rhovr_S=rhovr
    call mpi_allreduce(rhovr_S,rhovr_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rhovr=rhovr_R

    !> communicate rho*vr
    rho2vr_S=rho2vr
    call mpi_allreduce(rho2vr_S,rho2vr_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rho2vr=rho2vr_R
    
    sum_time = sum_time + dt
    
      
    !> Calculation of the R_photosphere (i.e. tau = 2./3.)  
    !> for finite disk correction  -- DD
    
    tau = 0.d0
    do ii = domain_nx1-1, 1, -1
      temp(ii) = (rad_e(ii)/const_rad_a)**0.25 
      t_rad0(ii) = max(temp(ii), 1.d4) 
      call set_opal_opacity(rho(ii),t_rad0(ii),n)
      kappa(ii)= n  !>unit_opacity
      kappa_0(ii) = max(kappa(ii), kappa_e*unit_opacity) !> lower limit is electron scattering
      tau(ii) = tau(ii+1) + 0.5*(kappa_0(ii)*rho(ii)+kappa_0(ii+1)*rho(ii+1))*dx_l1*unit_length*(R_lowerb/xx(ii))**2
      if (tau(ii) >= 2./3. .and. tau(ii+1) < 2./3.) then
            Rphoto0 = xx(ii)*unit_length     
      endif      
    enddo    
    Rph_dum = Rphoto0
   
    !> if global_time < 2, takes initial photospheric value, updates every
    !> ~2 global times until some global time(in base.par as static_rph), 
    !> after some gt it takes an average value for the photospheric radius -- DD
     
    if (global_time <= static_rph) then
      sum_radius = (sum_radius + Rph_dum)
      iter= iter + 1.
    else 
      sum_radius = sum_radius
      iter = iter  
    end if   
    avg_radius = sum_radius/iter
      
    if (global_time < 2.) then
      R_photosphere_dum = R_star*unit_length
      R_photosphere_var = R_star*unit_length !> to check the variability of the photosphere
    else if (2.<= global_time .and. global_time < static_rph .and. global_time-prev_time >= 2.) then 
      R_photosphere_dum = Rph_dum
      R_photosphere_var = Rph_dum
      prev_time= global_time
    else if (global_time < static_rph .and. global_time-prev_time < 2.) then
      R_photosphere_dum = R_photosphere_dum
      R_photosphere_var = R_photosphere_var
    else if (global_time >= static_rph) then
      !R_photosphere = avg_photo_radius
      R_photosphere_dum = avg_radius 
    end if   
    
    
    !if (mype==0) then
    !  R_photosphere_mype = R_photosphere_dum
    !end if
    
    !> communicate Rph
    Rph_S=R_photosphere_dum
    call mpi_allreduce(Rph_S,Rph_R,1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    R_photosphere=Rph_R/npe

    
    !R_photosphere = R_photosphere_mype 
    !print*, 'R_phot at', R_photosphere/R_sun
    !print*, 'R_phot / R_core ', R_photosphere/unit_length
    
    !> Only keep mdot in interested radii
    do ii = 1,np_mdot
      jj = jp_mdot(ii)
      p_mdot(ii) = mdot(jj)
      p_lum(ii) = lum(ii)
    enddo

    !> communicate mdot array
    mdot_S=p_mdot
    call mpi_reduce(mdot_S,mdot_R,np_mdot,&
    MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    p_mdot=mdot_R*unit_density*unit_velocity*unit_length**2/M_sun*year

    !> communicate lum array
    lum_S=p_lum
    call mpi_reduce(lum_S,lum_R,np_mdot,&
    MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    p_lum=lum_R*unit_radflux*unit_length**2/L_sun
    

    !> communicate highest amr level
    lvl_h_S=lvl_h
    call mpi_allreduce(lvl_h_S,lvl_h_R,domain_nx1,&
    MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
    lvl_h = lvl_h_R
    lvl_l_S=lvl_l
    call mpi_allreduce(lvl_l_S,lvl_l_R,domain_nx1,&
    MPI_INTEGER,MPI_MIN,icomm,ierrmpi)
    lvl_l = lvl_l_R

    !> integrate over time
    
    
    vr_sumt(1:domain_nx1) = vr_sumt(1:domain_nx1) + vr(1:domain_nx1)*dt
    rho_sumt(1:domain_nx1) = rho_sumt(1:domain_nx1) + rho(1:domain_nx1)*dt
    rho2_sumt(1:domain_nx1) = rho2_sumt(1:domain_nx1) + rho2(1:domain_nx1)*dt
    vr2_sumt(1:domain_nx1) = vr2_sumt(1:domain_nx1) + vr2(1:domain_nx1)*dt
    rhovr_sumt(1:domain_nx1) = rhovr_sumt(1:domain_nx1) + rhovr(1:domain_nx1)*dt
    rho2vr_sumt(1:domain_nx1) = rho2vr_sumt(1:domain_nx1) + rho2vr(1:domain_nx1)*dt
    rad_e_sumt(1:domain_nx1) = rad_e_sumt(1:domain_nx1) + rad_e(1:domain_nx1)*dt
    

    !> Write out average velocity profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_vr',status='replace')
      !write(unitanalysis,*) 'r | vr | <vr>t | vr**2 | <vr**2>t | dispersion | rho*vr | rho**2 * vr| lvl_h | lvl_l'
      do i=1,domain_nx1
        write(unitanalysis,'(8f25.12,2i4)') rr(i), vr(i), vr_sumt(i)/sum_time, vr2(i), vr2_sumt(i)/sum_time,dsqrt( abs(vr2_sumt(i)/sum_time - (vr_sumt(i)/sum_time)**2 )), rhovr_sumt(i)/rho_sumt(i), rho2vr_sumt(i)/rho2_sumt(i) ,lvl_h(i), lvl_l(i)
      enddo
    close(unitanalysis)
    endif

    !> Write out average density profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_rho',status='replace')
      write(unitanalysis,*) 'r | rho | <rho>t | rho**2 | <rho**2>t | <rho**2>t/<rho>t**2'
      do i=1,domain_nx1
        write(unitanalysis,'(6f30.16)') rr(i), rho(i), rho_sumt(i)/sum_time, rho2(i), rho2_sumt(i)/sum_time, sum_time*rho2_sumt(i)/rho_sumt(i)**2
      enddo
    close(unitanalysis)
    endif
    
        !> Write out average temperature profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_t',status='replace')
      write(unitanalysis,*) 'r | t_rad | <t_rad>t |'
      do i=1,domain_nx1
        write(unitanalysis,'(6f30.15)') rr(i), (rad_e(i)/const_rad_a)**0.25, ((rad_e_sumt(i)/sum_time)/const_rad_a)**0.25
      enddo
    close(unitanalysis)
    endif

    

    !> Write out average mdot at different radii
    if (mype==0) then
       if (global_time<smalldouble) then ! if very 1st iteration
         open(unit=unitanalysis,file=trim(base_filename)//'_mdot',status='replace')
         write(unitanalysis,'(a16,*(f7.2))') 'Mdot at t | r= ', rp_mdot
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_mdot
       else
         open(unit=unitanalysis,file=trim(base_filename)//'_mdot',access='append')
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_mdot
       endif
       close(unitanalysis)
     endif

    !> Write out average lum at different radii
    if (mype==0) then
       if (global_time<smalldouble) then ! if very 1st iteration
         open(unit=unitanalysis,file=trim(base_filename)//'_lum',status='replace')
         write(unitanalysis,'(a16,*(f7.2))') 'Lum at t | r= ', rp_mdot
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_lum
       else
         open(unit=unitanalysis,file=trim(base_filename)//'_lum',access='append')
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_lum
       endif
       close(unitanalysis)
     endif
     
     !> Write out photospheric radius with time
    if (mype==0) then
       if (global_time<smalldouble) then ! if very 1st iteration
         open(unit=unitanalysis,file=trim(base_filename)//'_photo',status='replace')
         write(unitanalysis,'(a16,*(f7.2))') 'time | Rph | Rph at all interation'
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, iter, R_photosphere/unit_length, avg_radius/unit_length, R_photosphere_var/unit_length, max(R_lowerb/unit_length,Rph_dum/unit_length)
       else
         open(unit=unitanalysis,file=trim(base_filename)//'_photo',access='append')
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, iter, R_photosphere/unit_length, avg_radius/unit_length, R_photosphere_var/unit_length, max(R_lowerb/unit_length,Rph_dum/unit_length)
       endif
       close(unitanalysis)
     endif


  end subroutine collapse_to_1D


  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                   :: g_rad(ixI^S), big_gamma(ixI^S)
    double precision                   :: g_grav(ixI^S)
    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: kappa(ixO^S), OPAL(ixO^S), CAK(ixO^S), CAK2(ixO^S), fdfac(ixO^S), beta_op(ixO^S)
    double precision                   :: vel(ixI^S), gradv(ixI^S)
    double precision                   :: rad_flux(ixO^S,1:ndim), Lum(ixO^S)
    double precision                   :: pp_rf(ixO^S), lambda(ixO^S), fld_R(ixO^S), gradOE(ixO^S)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)

    if (rhd_energy) call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

    call get_kappa_OPAL(ixI^L,ixO^L,w,x,OPAL)
    call get_kappa_CAK(ixI^L,ixO^L,w,x,CAK)
    call get_kappa_CAK2(ixI^L,ixO^L,w,x,CAK2)
    call get_fdfac(ixI^L,ixO^L,w,x,fdfac)
    call get_betaop(ixI^L,ixO^L,w,x,beta_op)

    g_rad(ixO^S) = (OPAL(ixO^S)+CAK2(ixO^S))*rad_flux(ixO^S,1)/(const_c/unit_velocity)
    g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient(vel,ixI^L,ixO^L,1,gradv)

    pp_rf(ixO^S) = two*rad_flux(ixO^S,1)/x(ixO^S,1)*dt

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)
    
    Lum(ixO^S) = 4*dpi*rad_flux(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun

    call gradientO(w(ixI^S,r_e),x,ixI^L,ixO^L,1,gradOE,nghostcells)

    w(ixO^S,i_v1) = w(ixO^S,mom(1))/w(ixO^S,rho_)
    w(ixO^S,i_v2) = w(ixO^S,mom(2))/w(ixO^S,rho_)
    if (rhd_energy) w(ixO^S,i_p) = (w(ixO^S,e_) - 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_)) &
          *(rhd_gamma - 1)

    w(ixO^S,i_Trad) = Trad(ixO^S)*unit_temperature
    if (rhd_energy) w(ixO^S,i_Tgas) = Tgas(ixO^S)*unit_temperature
    !w(ixO^S,i_Mdot) = 4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2*unit_density*unit_velocity/M_sun*year
    w(ixO^S,i_Opal) = OPAL(ixO^S)/kappa_e
    !w(ixO^S,i_CAK) = CAK(ixO^S)/kappa_e
    w(ixO^S,i_CAK2) = CAK2(ixO^S)/kappa_e
    w(ixO^S,i_lambda) = lambda(ixO^S)
    w(ixO^S,i_Edd) = lambda(ixO^S) + lambda(ixO^S)**2 * fld_R(ixO^S)**2
    w(ixO^S,i_Gamma) = big_gamma(ixO^S)
    w(ixO^S,i_Lum) = Lum(ixO^S)
    w(ixO^S,i_F1) = rad_flux(ixO^S,1)/F_bound
    w(ixO^S,i_F2) = rad_flux(ixO^S,2)/F_bound
    w(ixO^S,i_gradE) = gradOE(ixO^S)
    w(ixO^S,i_fdfac) = fdfac(ixO^S)
    w(ixO^S,i_beta_op) =  beta_op(ixO^S)
    

  end subroutine update_extravars

end module mod_usr
