!=======================================================================
!   fastran plasma state tools

!=======================================================================
! generate statefile
!
subroutine pstool_init(f_ps)

  use plasma_state_mod

  implicit none

  character*(*),intent(in) :: f_ps

  integer, parameter :: nspecmax = 10
  integer, parameter :: lun = 99
  integer, parameter :: iout = 6

  real(kind=rspec), parameter :: twopi = 6.283185307179586_rspec

  !-- label

  character*128 :: global_label = ''
  character*128 :: tokamak_id = ''
  character*128 :: runid = ''
  integer       :: shot_number = 0
  real*8        :: time = 0.0

  !-- species

  integer :: nspec_th     = 1
  integer :: nspec_beam   = 1
  integer :: nspec_fusion = 0
  integer :: nspec_rfmin  = 0
  integer :: nspec_gas    = 1

  integer, dimension(nspecmax) :: z_ion, a_ion
  integer, dimension(nspecmax) :: z_imp, a_imp
  integer, dimension(nspecmax) :: z_beam, a_beam
  integer, dimension(nspecmax) :: z_fusion, a_fusion
  integer, dimension(nspecmax) :: z_gas, a_gas
  integer, dimension(nspecmax) :: z_rfmin, a_rfmin

  !-- grid

  integer :: nrho         = 101
  integer :: nrho_eq      = 0 
  integer :: nth_eq       = 101
  integer :: nrho_eq_geo  = 0
  integer :: nrho_gas     = 0
  integer :: nrho_nbi     = 0
  integer :: nrho_ecrf    = 0
  integer :: nrho_icrf    = 0
  integer :: nrho_fus     = 0
  integer :: nrho_anom    = 0

  !-- namelist

  namelist /inpstool/   &
    global_label,       &
    tokamak_id,         &
    runid,              &
    shot_number,        &
    time,               &
    nspec_th,           &
    nspec_beam,         &
    nspec_fusion,       &
    nspec_rfmin,        &
    nspec_gas,          &
    z_ion,              &
    a_ion,              &
    z_imp,              &
    a_imp,              &
    z_beam,             &
    a_beam,             &
    z_fusion,           & 
    a_fusion,           &
    z_rfmin,            &
    a_rfmin,            &
    z_gas,              &
    a_gas,              &
    nrho,               &
    nrho_eq,            &
    nth_eq,             &
    nrho_eq_geo,        &
    nrho_gas,           &
    nrho_nbi,           &
    nrho_ecrf,          &
    nrho_icrf,          &
    nrho_fus,           &
    nrho_anom

  !-- local

  integer :: ierr
  integer :: i,ii

  !---------------------------------------------------------------------
  ! entry

  write(iout,*) 'Generate plasma state '
  write(iout,*) 'Using plasma state software version: ',trim(ps_version_id)

  ps_debug = 0 ! 0 for quiet operation

  !----------------------------------------------------------------
  ! read inpstool

  open(unit=lun,file='inpstool',status='old',iostat=ierr)
  if(ierr.ne.0) then
     write(iout,*) ' namelist file open failure','inps'
     stop
     return
  endif
  read(lun,inpstool)
  close(lun)

  !---------------------------------------------------------------------
  ! global label

  ps%global_label = global_label
  ps%runid = runid

  ps%tokamak_id = tokamak_id
  ps%shot_number = shot_number

  ps%t0 = time ! current time
  ps%t1 = time ! dummy

  !---------------------------------------------------------------------
  ! dimension

  ! plasma species

  ps%nspec_th = nspec_th ! main ion + impurity
  ps%nspec_rfmin =  nspec_rfmin  ! no minority specie
  ps%nspec_fusion = nspec_fusion  ! fusion specie
  ps%ngsc0 = nspec_gas
  ps%nspec_beam = nspec_beam

  ! plasma component

  ps%nrho = nrho  ! grid for the plasma component

  if (nrho_eq     == 0)  nrho_eq = nrho   
  if (nrho_eq_geo == 0)  nrho_eq_geo = nrho
  if (nrho_gas    == 0)  nrho_gas = nrho 
  if (nrho_nbi    == 0)  nrho_nbi = nrho 
  if (nrho_ecrf   == 0)  nrho_ecrf = nrho
  if (nrho_icrf   == 0)  nrho_icrf = nrho 
  if (nrho_fus    == 0)  nrho_fus = nrho 
  if (nrho_anom   == 0)  nrho_anom = nrho

  ! equilibrium
  ! r and z grid sizes will be reset when first g-eqdsk is read.

  ps%nrho_eq     = nrho_eq
  ps%nth_eq      = nth_eq
  ps%nr          = 0
  ps%nz          = 0
  ps%nrho_eq_geo = nrho_eq_geo

  ! neutral

  ps%nrho_gas = nrho_gas

  ! nbi

  ps%nrho_nbi = nrho_nbi

  ! ech

  ps%nrho_ecrf = nrho_ecrf

  ! ic

  ps%nrho_icrf = nrho_icrf

  ! fusion

  ps%nrho_fus = nrho_fus

  ! anomalous diffusion

  ps%nrho_anom = nrho_anom

  !---------------------------------------------------------------------
  ! initial allocation of state arrays

  call ps_alloc_plasma_state(ierr)
  call ckerr('ps_alloc_plasma_state',ierr)

  !---------------------------------------------------------------------
  ! initialize grid

  do i = 1,ps%nrho
      ps%rho(i) = (i-1.0)/(ps%nrho-1.0)
  enddo

  do i = 1,ps%nrho_eq
      ps%rho_eq(i) = (i-1.0)/(ps%nrho_eq-1.0)
  enddo

  do i=1,ps%nth_eq
    ps%th_eq(i) = twopi*(i-1)/(nth_eq-1.0)
  enddo

  do i = 1,ps%nrho_eq_geo
      ps%rho_eq_geo(i) = (i-1.0)/(ps%nrho_eq_geo-1.0)
  enddo

  do i = 1,ps%nrho_gas
      ps%rho_gas(i) = (i-1.0)/(ps%nrho_gas-1.0)
  enddo

  do i = 1,ps%nrho_nbi
      ps%rho_nbi(i) = (i-1.0)/(ps%nrho_nbi-1.0)
  enddo

  do i = 1,ps%nrho_ecrf
      ps%rho_ecrf(i) = (i-1.0)/(ps%nrho_ecrf-1.0)
  enddo

  do i = 1,ps%nrho_icrf
      ps%rho_icrf(i) = (i-1.0)/(ps%nrho_icrf-1.0)
  enddo

  do i = 1,ps%nrho_fus
      ps%rho_fus(i) = (i-1.0)/(ps%nrho_fus-1.0)
  enddo

  do i = 1,ps%nrho_anom
      ps%rho_anom(i) = (i-1.0)/(ps%nrho_anom-1.0)
  enddo

  !---------------------------------------------------------------------
  ! species

  ! electron

  call ps_species_convert(-1, -1, 0, &
       ps%qatom_s(0), ps%q_s(0), ps%m_s(0), ierr)  ! electron (species 0)
  call ckerr('ps_species_convert',ierr)

  ! thermal ions (main + impurities)

  do i = 1, nspec_th
    call ps_species_convert(z_ion(i), z_ion(i), a_ion(i), &
         ps%qatom_s(i), ps%q_s(i), ps%m_s(i), ierr)
    call ckerr('ps_species_convert',ierr)
  enddo

  ! beam

  do i = 1, nspec_beam
    call ps_species_convert(z_beam(i), z_beam(i), a_beam(i), &
         ps%qatom_snbi(i), ps%q_snbi(i), ps%m_snbi(i), ierr)
    call ckerr('ps_species_convert',ierr)
  enddo

  ! fusion

  do i = 1, nspec_fusion
    call ps_species_convert(z_fusion(i), z_fusion(i), a_fusion(i), &
         ps%qatom_sfus(i), ps%q_sfus(i), ps%m_sfus(i), ierr)
    call ckerr('ps_species_convert',ierr)
  enddo

  ! minority

  do i = 1, nspec_rfmin
    call ps_species_convert(z_rfmin(i), z_rfmin(i), a_rfmin(i), &
         ps%qatom_rfmin(i), ps%q_rfmin(i), ps%m_rfmin(i), ierr)  
    call ckerr('ps_species_convert',ierr)
  enddo

  ! species list

  call ps_label_species(ierr)
  call ckerr('ps_label_species',ierr)

  call ps_merge_species_lists(ierr)
  call ckerr('ps_merge_species_lists',ierr)

  ! neutral (deutron assumed)

  do i = 1, nspec_gas
     if ( z_gas(i) .eq. 1 .and. a_gas(i) .eq. 1 ) then
        ps%gas_atom(i) = 'H'
        ps%gs_name(i) = 'Hrec'
     else if ( z_gas(i) .eq. 1 .and. a_gas(i) .eq. 2 ) then
        ps%gas_atom(i) = 'D'
        ps%gs_name(i) = 'Drec'
     else if ( z_gas(i) .eq. 1 .and. a_gas(i) .eq. 3 ) then
        ps%gas_atom(i) = 'T'
        ps%gs_name(i) = 'Trec'
     else if ( z_gas(i) .eq. 2 .and. a_gas(i) .eq. 3 ) then
        ps%gas_atom(i) = 'He3'
        ps%gs_name(i) = 'He3rec'
     else if ( z_gas(i) .eq. 2 .and. a_gas(i) .eq. 4 ) then
        ps%gas_atom(i) = 'He4'
        ps%gs_name(i) = 'He4rec'
     else
        call ckerr('gas type not supported',99)
     endif
  enddo

  call ps_neutral_species(ierr)
  call ckerr('ps_neutral_species',ierr)

  !  print out

  do ii=0,ps%nspec_all
     write(iout,1) ii,ps%all_type(ii), &
          trim(ps%all_name(ii)),ps%q_all(ii),ps%m_all(ii)
  enddo
1 format(' Specie index & type: ',i2,1x,i2,1x, '"',a,'" charge & mass: ',2(1pe12.5,1x))

  !----------------------------------------------------------------
  ! lock plasma state

  ! write(iout,*) 'locked machine description (mdescr).'
  ! write(iout,*) 'locked shot configuration (sconfig).'

  ! call ps_update_hashCodes(ierr, mdescr_lock=.TRUE., sconfig_lock=.TRUE.)

  !----------------------------------------------------------------
  ! write plasma state

  write(iout,*) 'updating memory'

  call ps_state_memory_update(ierr,eqcheck=.FALSE.)
  call ckerr('ps_state_memory_update',ierr)

  write(iout,*) 'storing plasma state'

  call ps_store_plasma_state(ierr,trim(f_ps))
  call ckerr('ps_store_plasma_state',ierr)


end

!=======================================================================
! generate statefile for geqdsk process
!
subroutine pstool_geqdsk(f_ps, file_geqdsk, bdy_crat)

  use plasma_state_mod

  implicit none

  integer, parameter :: lun = 99
  integer, parameter :: iout = 6

  character*128 :: f_ps ! = 'geqdsk'
  character*128 :: file_geqdsk ! = 'geqdsk'


  real*8 :: bdy_crat != 1.0d-6

  !-- label

  character*128 :: global_label = ''
  character*128 :: tokamak_id = ''
  character*128 :: runid = ''
  integer       :: shot_number = 0
  real*8        :: time = 0.0

  integer :: nrho_eq      = 51
  integer :: nth_eq       = 101
  integer :: nrho_eq_geo  = 51

  !-- namelist

  namelist /inpstool/       &
    global_label,       &
    tokamak_id,         &
    runid,              &
    shot_number,        &
    time,               &
    nrho_eq,            &
    nth_eq,             &
    nrho_eq_geo

  !-- local

  integer :: ierr,i, ii

  !---------------------------------------------------------------------
  ! entry

  write(iout,*) 'Generate plasma state v',trim(ps_version_id)

  ps_debug = 0 ! 0 for quiet operation

  !----------------------------------------------------------------
  ! read inps

  open(unit=lun,file='inpstool',status='old',iostat=ierr)
  if(ierr.ne.0) then
     write(iout,*) ' namelist file open failure','inpstool'
     stop
     return
  endif
  read(lun,inpstool)
  close(lun)

  !---------------------------------------------------------------------
  ! global label

  ps%global_label = global_label
  ps%runid = runid

  ps%tokamak_id = tokamak_id
  ps%shot_number = shot_number

  ps%t0 = time ! current time
  ps%t1 = time ! dummy

  !---------------------------------------------------------------------
  ! dimension

  ps%nrho_eq     = nrho_eq
  ps%nth_eq      = nth_eq
  ps%nr          = 0
  ps%nz          = 0
  ps%nrho_eq_geo = nrho_eq_geo

  !---------------------------------------------------------------------
  ! initial allocation of state arrays

  call ps_alloc_plasma_state(ierr)
  call ckerr('ps_alloc_plasma_state',ierr)

  !---------------------------------------------------------------------
  ! initialize grid

  do i = 1,ps%nrho_eq
      ps%rho_eq(i) = (i-1.0)/(ps%nrho_eq-1.0)
  enddo

  do i=1,ps%nth_eq
    ps%th_eq(i) = ps_twopi*(i-1)/(nth_eq-1.0)
  enddo

  do i = 1,ps%nrho_eq_geo
      ps%rho_eq_geo(i) = (i-1.0)/(ps%nrho_eq_geo-1.0)
  enddo

  !----------------------------------------------------------------
  ! equilibrium from geqdsk

  write(iout,*) 'loading plasma state equilibrium from g file'
  write(iout,*) trim(file_geqdsk)

  call ps_update_equilibrium(ierr, g_filepath=trim(file_geqdsk), &
       bdy_crat = bdy_crat )
  call ckerr('ps_update_equilibrium',ierr)

  !----------------------------------------------------------------
  ! flux surface averages, etc., derived from equilibrium

  write(iout,*) 'computing equilibrium derived profiles'

  write(iout,*) 'flux surface computations'
  call ps_mhdeq_derive('everything',ierr)
  call ckerr('ps_mhdeq_derive',ierr)

  ! write(iout,*) 'additional metric'
  ! allocate(metric(ps%nrho_eq_geo))
  ! call ps_rho_metric ("rho_eq_geo","<SQRT(1-B/BMAX)>S",metric,ierr)
  ! ps%gncfh = metric

  !----------------------------------------------------------------
  ! write plasma state

  write(iout,*) 'updating memory'

  call ps_state_memory_update(ierr,eqcheck=.FALSE.)
  call ckerr('ps_state_memory_update',ierr)

  write(iout,*) 'storing plasma state'

  call ps_store_plasma_state(ierr,'ps.nc')
  call ckerr('ps_store_plasma_state',ierr)


  !----------------------------------------------------------------
  ! lock plasma state

  ! write(iout,*) 'locked machine description (mdescr).'
  ! write(iout,*) 'locked shot configuration (sconfig).'

  ! call ps_update_hashCodes(ierr, mdescr_lock=.TRUE., sconfig_lock=.TRUE.)

  !----------------------------------------------------------------
  ! write plasma state

  write(iout,*) 'updating memory'

  call ps_state_memory_update(ierr,eqcheck=.FALSE.)
  call ckerr('ps_state_memory_update',ierr)

  write(iout,*) 'storing plasma state'

  call ps_store_plasma_state(ierr,f_ps)
  call ckerr('ps_store_plasma_state',ierr)

end


subroutine pstool_load_rf()

  use plasma_state_mod

  implicit none

  integer :: ierr,ii
  character*32 :: icrf_src_name
  real*8 :: power_ic
  integer, parameter :: lun = 99
  integer, parameter :: iout = 6

  namelist/intoric/                           &
     power_ic

  !----------------------------------------------------------------
  ! read intoric

  open(unit=lun,file='intoric',status='old',iostat=ierr)
  if(ierr.ne.0) then
     write(iout,*) ' namelist file open failure','intoric'
     stop
     return
  endif
  read(lun,intoric)
  close(lun)

  call ps_get_plasma_state(ierr,'ps.nc')

  ps%nicrf_src = 1
  call ps_alloc_plasma_state(ierr)

  ps%power_ic(1) = power_ic

  do ii=1,ps%nicrf_src

     if (ii .LT. 10) then
         write(icrf_src_name,'(A,I1)') 'ICRF',ii
     else
         write(icrf_src_name,'(A,I2)') 'ICRF',ii
     endif
     ps%icrf_src_name(ii) = icrf_src_name
  end do

  call ps_store_plasma_state(ierr,'ps.nc')

end

subroutine pstool_load_nbi(f_ps)

  use plasma_state_mod

  implicit none

  character*(*), intent(in) :: f_ps

  integer, parameter :: nbmax = 32
  integer, parameter :: lun = 99
  integer, parameter :: iout = 6

  real(kind=rspec), parameter :: twopi = 6.283185307179586_rspec

  !-- nbi power, ...

  real*8 :: dt_nubeam = 0.050
  integer :: nstep = 1
  integer :: navg = 0

  real*8,dimension(nbmax) :: pinja,einja,ffulla,fhalfa

  !--nbi configuration

  integer                  :: nbeam
  logical,dimension(nbmax) :: nlco
  real*8 ,dimension(nbmax) :: abeama,xzbeama
  integer,dimension(nbmax) :: nbshapa
  real*8 ,dimension(nbmax) :: bmwidra,bmwidza
  real*8 ,dimension(nbmax) :: rtcena,xlbtna,xybsca
  real*8 ,dimension(nbmax) :: divra,divza,foclza,foclra
  integer,dimension(nbmax) :: nbapsha
  real*8 ,dimension(nbmax) :: rapedga,xzpedga,xlbapa,xybapa
  integer,dimension(nbmax) :: nbapsh2
  real*8 ,dimension(nbmax) :: xrapoffa, xzapoffa
  real*8 ,dimension(nbmax) :: rapedg2,xzpedg2,xlbapa2,xrapoff2, xzapoff2
  real*8 ,dimension(nbmax) :: xbzeta

  real*8 :: difb_0,difb_a,difb_in,difb_out

  !--namelist

  namelist/nbi_config/                           &
     nbeam,                                      &
     abeama,xzbeama,                             &
     pinja,einja,ffulla,fhalfa,                  &
     nlco,                                       &
     nbshapa,                                    &
     bmwidra,bmwidza,                            &
     rtcena,xlbapa,xlbtna,                       &
     divra,divza,foclza,foclra,                  &
     nbapsha,                                    &
     xybapa,xybsca,                              &
     rapedga,xzpedga,                            &
     nbapsh2,                                    &
     rapedg2,xzpedg2,xlbapa2,                    &
     xrapoffa, xzapoffa,xrapoff2,xzapoff2,       &
     xbzeta

  namelist/nbi_model/ &
     difb_0,difb_a,difb_in,difb_out

  namelist/nubeam_run/ &
     dt_nubeam, nstep, navg


  !-- local

  character*32 :: nbi_src_name
  real(kind=rspec) :: zdivcon
  integer ierr,ii

  !----------------------------------------------------------------
  ! initialize

  nbeam = 0
  nlco = .true.
  abeama = 0.0; xzbeama = 0.0
  nbshapa = 0
  bmwidra = 0.0; bmwidza = 0.0
  rtcena = 0.0; xlbtna = 0.0; xybsca = 0.0
  divra = 0.0; divza = 0.0; foclza = 0.0; foclra = 0.0
  nbapsha = 0
  rapedga = 0.0; xzpedga = 0.0; xlbapa = 0.0; xybapa = 0.0
  xrapoffa = 0.0; xzapoffa = 0.0
  nbapsh2 = 0
  rapedg2 = 0.0; xzpedg2 = 0.0; xlbapa2 = 0.0
  xrapoff2 = 0.0; xzapoff2 = 0.0;
  xbzeta = 0.0

  !----------------------------------------------------------------
  ! read innubeam

  open(unit=lun,file='innubeam',status='old',iostat=ierr)
  if(ierr.ne.0) then
     write(iout,*) ' namelist file open failure','innubeam'
     stop
     return
  endif
  read(lun,nubeam_run)
  read(lun,nbi_config)
  close(lun)

  !---------------------------------------------------------------------
  ! entry

  call ps_get_plasma_state(ierr,trim(f_ps))

  ps%t0 = 0.0
  ps%t1 = dt_nubeam

  !----------------------------------------------------------------
  ! number of beams and memory allocation

  ps%nbeam = nbeam
  call ps_alloc_plasma_state(ierr)
  call ckerr('ps_alloc_plasma_state',ierr)

  do ii=1,ps%nbeam

     if ( int(xzbeama(ii)) .eq. 1 .and. int(abeama(ii)) .eq. 1 ) then
         ps%nbion(ii) = 'H'
     else if ( int(xzbeama(ii)) .eq. 1 .and. int(abeama(ii)) .eq. 2 ) then
         ps%nbion(ii) = 'D'
     else if ( int(xzbeama(ii)) .eq. 2 .and. int(abeama(ii)) .eq. 3 ) then
         ps%nbion(ii) = 'He3'
     else if ( int(xzbeama(ii)) .eq. 2 .and. int(abeama(ii)) .eq. 4 ) then
         ps%nbion(ii) = 'He4'
     else
        call ckerr('beam type not supported',99)
     endif

     write(iout,'(A,I2,X,A)') 'beam', ii, ps%nbion(ii)

     if (ii .LT. 10) then
         write(nbi_src_name,'(A,I1)') 'NB',ii
     else
         write(nbi_src_name,'(A,I2)') 'NB',ii
     endif
     ps%nbi_src_name(ii) = nbi_src_name
  end do

  !----------------------------------------------------------------
  ! nbi geometry

  do ii=1,ps%nbeam

    ps%srtcen(ii) = 1.0d-2*rtcena(ii)
    if (nlco(ii) .eqv. .false.) then
       ps%srtcen(ii) = -ps%srtcen(ii)   ! (signed: + means ccw momentum inj.)
    endif

    ps%lbsctan(ii) =  1.0d-2*xlbtna(ii) ! dist., sce to tangency pt.
    ps%zbsc(ii)    =  1.0d-2*xybsca(ii) ! height, sce above midplane
    ps%phibsc(ii)  =  xbzeta(ii)        ! toroidal angle of sce

    ps%lbscap(ii)  =  1.0d-2*xlbapa(ii) ! dist., sce to aperture
    ps%zbap(ii)    =  1.0d-2*xybapa(ii) ! height, aperture center

    zdivcon = twopi/(360.0d0*sqrt(2.0d0))
    if ( nbshapa(ii) .eq. 1 ) then
        ps%nbshape(ii)         = 'rectangle'        ! shape of source
        ps%b_halfwidth(ii)     = 1.0d-2*bmwidra(ii) ! half-width
        ps%b_halfheight(ii)    = 1.0d-2*bmwidza(ii) ! half-height
        ps%b_hfocal_length(ii) = 1.0d-2*foclra(ii)  ! horiz. focal len
        ps%b_vfocal_length(ii) = 1.0d-2*foclza(ii)  ! vert. focal len
        ps%b_hdivergence(ii)   = divra(ii)/zdivcon  ! horiz. diverg.
        ps%b_vdivergence(ii)   = divza(ii)/zdivcon  ! vert.  diverg.
    else if ( nbshapa(ii) .eq. 2 ) then
        ps%nbshape(ii)         = 'circle'
        ps%b_halfwidth(ii)     = 1.0d-2*bmwidra(ii)
        ps%b_hfocal_length(ii) = 1.0d-2*foclra(ii)
        ps%b_hdivergence(ii)   = divra(ii)
    else
        write(iout,*) 'nbshapa = 1 or 2'
        ierr = 1
        call ckerr('nbshapa',ierr)
    endif

    if(nbapsha(ii) .eq. 1) then
        ps%nbap_shape(ii)    = 'rectangle'
        ps%ap_halfwidth(ii)  =  1.0d-2*rapedga(ii)
        ps%ap_halfheight(ii) =  1.0d-2*xzpedga(ii)
        ps%ap_horiz_offset(ii)= 1.0d-2*xrapoffa(ii)
        ps%ap_vert_offset(ii) = 1.0d-2*xzapoffa(ii)
    else if(nbapsha(ii) .eq. 2) then
        ps%nbap_shape(ii) = 'circle'
        ps%ap_halfwidth(ii) = 1.0d-2*rapedga(ii)
    else
        write(iout,*) 'nbapsha = 1 or 2'
        ierr = 1
        call ckerr('nbshapa',ierr)
    endif

    if(nbapsh2(ii) .eq. 0) THEN
        ps%nbap2_shape(ii) = 'None'
    else if(nbapsh2(ii) .eq. 1) then
        ps%nbap2_shape(ii)   = 'Rectangle'
        ps%Lbscap2(ii)       = 1.0d-2*xlbapa2(ii)
        ps%ap2_halfwidth(ii) = 1.0d-2*rapedg2(ii)
        ps%ap2_halfheight(ii)= 1.0d-2*xzpedg2(ii)
        ps%ap2_horiz_offset(ii)= 1.0d-2*xrapoff2(ii)
        ps%ap2_vert_offset(ii) = 1.0d-2*xzapoff2(ii)
    else if(nbapsh2(ii) .eq. 2) then
        ps%nbap2_shape(ii)   = 'Circle'
        ps%Lbscap2(ii)       = 1.0d-2*xlbapa2(ii)
        ps%ap2_halfwidth(ii) = 1.0d-2*rapedg2(ii)
    else
        write(iout,*) 'nbapsh2 = 0 or 1 or 2'
        ierr = 1
        call ckerr('nbapsh2',ierr)
    endif

    ps%beam_type(ii) = 'Standard'

  end do

  !----------------------------------------------------------------
  ! nbi power

  do ii=1,ps%nbeam

    ps%power_nbi(ii) = pinja(ii)  ! power on each beam source
    ps%kvolt_nbi(ii) = 1.0d-3*einja(ii)  ! energy of each beam source **kev**
    ps%frac_full(ii) = ffulla(ii)  ! fraction of beam current at full voltage
    ps%frac_half(ii) = fhalfa(ii)  ! fraction of beam current at half voltage

  end do

  do ii=1,ps%nbeam
      ps%einj_max(ii) = 1.25*1.0d-3*einja(ii)  ! max einj of beam
      ps%einj_min(ii) = 20.0  ! min einj of beam
  enddo

  ps%dn0out = 0.5e18  !1.0e6*dn0out, defalut:0.5e18

  !----------------------------------------------------------------
  ! write plasma state

  write(iout,*) 'storing plasma state'

  call ps_store_plasma_state(ierr,f_ps)
  call ckerr('ps_store_plasma_state',ierr)

end

subroutine pstool_load_geqdsk(f_ps, file_geqdsk,bdy_crat)

  use plasma_state_mod

  implicit none

  character*(*), intent(in) :: f_ps
  character*(*), intent(in) :: file_geqdsk
  !character*128 :: file_geqdsk = 'geqdsk'
  real*8 :: bdy_crat != 1.0d-6

  integer, parameter :: iout = 6

  !----------------------------------------------------------------
  ! local

  real(kind=rspec), dimension(:), allocatable :: metric

  integer :: nrho_zc
  real*8, dimension(:), allocatable :: vec_zc,rho_zc
  real*8 :: a,b,c,d
  real*8:: zne_adj,zzne_adj,nion,nimp
  integer :: i,k,ierr

  !----------------------------------------------------------------
  ! read plasma state

  call ps_get_plasma_state(ierr,trim(f_ps))

  !----------------------------------------------------------------
  ! equilibrium from geqdsk

  write(iout,*) 'loading plasma state equilibrium from g file'
  write(iout,*) trim(file_geqdsk)

  call ps_update_equilibrium(ierr, g_filepath=trim(file_geqdsk), &
       bdy_crat = bdy_crat )
  call ckerr('ps_update_equilibrium',ierr)

  !----------------------------------------------------------------
  ! flux surface averages, etc., derived from equilibrium

  write(iout,*) 'computing equilibrium derived profiles'

  write(iout,*) 'flux surface computations'
  call ps_mhdeq_derive('everything',ierr)
  call ckerr('ps_mhdeq_derive',ierr)

  !write(iout,*) 'additional metric'
  !allocate(metric(ps%nrho_eq_geo))
  !call ps_rho_metric ("rho_eq_geo","<SQRT(1-B/BMAX)>S",metric,ierr)
  !ps%gncfh = metric

  !----------------------------------------------------------------
  ! write plasma state

  write(iout,*) 'updating memory'

  call ps_state_memory_update(ierr,eqcheck=.FALSE.)
  call ckerr('ps_state_memory_update',ierr)

  write(iout,*) 'storing plasma state'

  call ps_store_plasma_state(ierr,f_ps)
  call ckerr('ps_store_plasma_state',ierr)

end

subroutine pstool_load_instate()

  use plasma_state_mod

  implicit none

  integer, parameter :: nspecmax = 10
  integer, parameter :: nrhomax  = 201
  integer, parameter :: nbdrymax = 201
  integer, parameter :: lun = 99
  integer, parameter :: iout = 6

  !--instate

  character*128 :: tokamak_id = ''
  real*8  :: ip, b0, r0, rmajor,aminor, kappa, delta
  integer :: nrho
  real*8,dimension(nrhomax)  :: rho
  real*8,dimension(nrhomax)  :: ne,te,ti,omega,zeff
  integer                    :: n_ion, n_imp, n_beam
  integer,dimension(nspecmax):: z_ion, a_ion
  integer,dimension(nspecmax):: z_imp, a_imp
  real*8,dimension(nspecmax) :: f_ion
  real*8,dimension(nspecmax) :: f_imp
  real*8,dimension(nrhomax)  :: j_tot, j_oh, j_bs, j_nb, j_ec, j_ic, j_ext
  real*8,dimension(nrhomax)  :: q, psipol
  real*8,dimension(nrhomax)  :: pe_nb, pe_ec, pe_ic, pe_fus, pe_ionization
  real*8,dimension(nrhomax)  :: pi_nb, pi_ec, pi_ic, pi_fus, pi_ionization, pi_cx
  real*8,dimension(nrhomax)  :: p_rad, p_ohm, p_ei
  real*8,dimension(nrhomax)  :: torque_nb, torque_in
  real*8,dimension(nrhomax)  :: se_nb, se_ionization
  real*8,dimension(nrhomax)  :: si_nb, si_ionization
  real*8,dimension(nrhomax)  :: density_beam
  real*8,dimension(nrhomax)  :: density_alpha
  real*8,dimension(nrhomax)  :: wbeam
  real*8,dimension(nrhomax)  :: walpha
  real*8,dimension(nrhomax)  :: chie,chii
  real*8,dimension(nrhomax)  :: p_eq
  integer                    :: nbdry, nlim
  real*8,dimension(nbdrymax) :: rbdry,zbdry,rlim,zlim

  !--namelist

  namelist /instate/                              &
      ip, b0, r0, rmajor,aminor, kappa, delta,    &
      nrho,                                       &
      rho,                                        &
      ne,te,ti,omega,zeff,                        &
      n_ion, n_imp, n_beam,                       &
      z_ion, a_ion, f_ion,                        &
      z_imp, a_imp, f_imp,                        &
      j_tot, j_oh, j_bs, j_nb, j_ec, j_ic, j_ext, &
      q, psipol,                                  &
      pe_nb, pe_ec, pe_ic, pe_fus, pe_ionization,        &
      pi_nb, pi_ec, pi_ic, pi_fus, pi_ionization, pi_cx, &
      p_rad, p_ohm, p_ei,                         &
      torque_nb, torque_in,                       &
      se_nb, se_ionization,                       &
      si_nb, si_ionization,                       &
      density_beam,                               &
      density_alpha,                              &
      wbeam,                                      &
      walpha,                                     &
      tokamak_id,                                 &
      chie,chii,                                  &
      nbdry, rbdry, zbdry, nlim, rlim, zlim,      &
      p_eq

  !----------------------------------------------------------------
  ! local

  real*8, dimension(:),allocatable :: metric

  integer :: nrho_zc
  real*8, dimension(:), allocatable :: vec_zc,rho_zc
  real*8 :: a,b,c,d
  real*8:: zne_adj,zzne_adj,nion,nimp
  integer :: i,k,ierr

  !----------------------------------------------------------------
  ! read instate

  open(unit=lun,file='instate',status='old',iostat=ierr)
  if(ierr.ne.0) then
     write(iout,*) ' namelist file open failure','instate'
     stop
     return
  endif
  read(lun,instate)
  close(lun)

  !----------------------------------------------------------------
  ! entry

  call ps_get_plasma_state(ierr,'ps.nc')

  !----------------------------------------------------------------
  ! kinetic profiles

  ! zone centered grid (construct local copy)

  nrho_zc = ps%nrho - 1
  allocate(vec_zc(nrho_zc))
  allocate(rho_zc(nrho_zc))

  rho_zc = (ps%rho(1:nrho_zc) + ps%rho(2:ps%nrho))/2

  !  electron density

  vec_zc = 0.5*(ne(1:nrho_zc)+ne(2:ps%nrho))
  ps%ns(:,0) = vec_zc(1:nrho_zc)*1.0e19

  ! ion density

  a=0; b=0; c=0; d=0
  do k = 1, n_imp
     b = b+f_imp(k)*z_imp(k)
     d = d+f_imp(k)*z_imp(k)*z_imp(k)
  end do
  do k = 1, n_ion
     a = a+f_ion(k)*z_ion(k)
     c = c+f_ion(k)*z_ion(k)*z_ion(k)
  end do

  do i=1,nrho_zc

     ! depletion due to RF minority species

     zne_adj = ps%ns(i,0)
     zzne_adj = ps%ns(i,0)*zeff(i)

     ! depletion due to beam ions

     zne_adj = zne_adj - 1.0*density_beam(i)
     zzne_adj = zzne_adj - 1.0**2*density_beam(i)

     ! effective main ion and impurity densities

     nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
     nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

     do k = 1, n_ion
         ps%ns(i,k) = f_ion(k)*nion
     enddo

     do k = 1, n_imp
         ps%ns(i,n_ion+k) = f_imp(k)*nimp
     enddo

     ! total thermal ions

     ps%ni(i) = 0.0

     do k = 1, n_ion+n_imp
        ps%ni(i) = ps%ni(i) +  ps%ns(i,k)
     end do

  enddo

  ! temperatures

  vec_zc = 0.5*(te(1:nrho_zc)+te(2:ps%nrho))
  ps%TS(:,0) = vec_zc(1:nrho_zc)

  vec_zc = 0.5*(ti(1:nrho_zc)+ti(2:ps%nrho))
  do k = 1, n_ion+n_imp
    ps%TS(:,k) = vec_zc(1:nrho_zc)
  end do
  ps%Ti(:) = vec_zc(1:nrho_zc)

  ! rotation

  vec_zc = 0.5*(omega(1:nrho_zc)+omega(2:ps%nrho))
  ps%omegat(:) = vec_zc(1:nrho_zc)

  ! zeff

  vec_zc = 0.5*(zeff(1:nrho_zc)+zeff(2:ps%nrho))
  ps%zeff(:) = vec_zc(1:nrho_zc)
  ps%zeff_th(:) = vec_zc(1:nrho_zc)

  !----------------------------------------------------------------
  ! lock plasma state

  ! write(iout,*) 'locked machine description (mdescr).'
  ! write(iout,*) 'locked shot configuration (sconfig).'

  ! call ps_update_hashCodes(ierr, mdescr_lock=.TRUE., sconfig_lock=.TRUE.)

  write(iout,*) 'updating memory'

  call ps_state_memory_update(ierr,eqcheck=.FALSE.)
  call ckerr('ps_state_memory_update',ierr)

  !----------------------------------------------------------------
  ! write plasma state

  write(iout,*) 'storing plasma state'

  call ps_store_plasma_state(ierr,'ps.nc')
  call ckerr('ps_store_plasma_state',ierr)

end

subroutine pstool_dump_geqdsk()

  use plasma_state_mod

  implicit none

  integer :: ierr

  !----------------------------------------------------------------
  ! entry

  call ps_get_plasma_state(ierr,'ps.nc')

  !----------------------------------------------------------------
  ! entry

  call ps_wr_geqdsk(ierr,'./geqdsk',ps)

end

subroutine pstool_rehash(f_ps)

  use plasma_state_mod
  implicit none
  character*(*),intent(in) :: f_ps
  integer :: ierr
  call ps_get_plasma_state(ierr,trim(f_ps))
  call ps_update_hashCodes(ierr, mdescr_lock=.TRUE., sconfig_lock=.TRUE.)
  call ps_store_plasma_state(ierr,trim(f_ps))
  call ckerr('ps_store_plasma_state',ierr)

end

!======================================================================
! utils
!
subroutine ckerr(sbrtn,ierr)
  integer, intent(in):: ierr
  character*(*), intent(in) :: sbrtn
  if(ierr.ne.0) then
     write(6,*) ' error in '//trim(sbrtn)
     stop
  endif
end subroutine ckerr

!======================================================================
! main
!
program pstool

  implicit none

  character(len=256) :: command, parm1, parm2, parm3
  integer :: ierr,narg
  real*8 :: bdy_crat

  !----------------------------------------------------------------
  ! input argument

  narg = command_argument_count()

  call get_command_argument(1,command)
  call get_command_argument(2,parm1)
  call get_command_argument(3,parm2)
  call get_command_argument(4,parm3)

  !----------------------------------------------------------------
  ! call

!  if (trim(command) .eq. 'init' ) then
!
!     write(*,*) 'pstool: init'
!     call pstool_init(parm1)
!
!  else if (trim(command) .eq. 'geqdsk' ) then
!     read(parm1,'(E15.7)') bdy_crat
!        write(*,'(A,E15.7)') 'geqdsk',bdy_crat
!     call pstool_geqdsk(bdy_crat)
!
!  else if (trim(command) .eq. 'load' ) then
!
!     if (trim(parm1) .eq. 'innubeam' ) then
!        write(*,'(A)') 'load innubeam'
!        call pstool_load_nbi()
!     else if (trim(parm1) .eq. 'geqdsk' ) then
!        read(parm2,'(E15.7)') bdy_crat
!        write(*,'(A,E15.7)') 'load geqdsk',bdy_crat
!        call pstool_load_geqdsk(bdy_crat)
!     else if (trim(parm1) .eq. 'instate' ) then
!        write(*,'(A)') 'load instate'
!        call pstool_load_instate()
!     else if (trim(parm1) .eq. 'intoric' ) then
!        write(*,'(A)') 'load intoric'
!        call pstool_load_rf()
!     endif
!
!  else if (trim(command) .eq. 'dump' ) then
!     if (trim(parm1) .eq. 'geqdsk' ) then
!        write(*,'(A)') 'dump geqdsk'
!        call pstool_dump_geqdsk()
!     endif
!
!  else if (trim(command) .eq. 'rehash' ) then
!     call pstool_rehash(parm1)
!
!  endif

   select case (trim(command))

      case ('init')

         write(*,*) 'pstool: init'
         call pstool_init(parm1)

      case ('rehash')

         write(*,*) 'pstool: rehash'
         call pstool_rehash(parm1)

      case ('geqdsk')

         read(parm3,'(E15.7)') bdy_crat
         write(*,*) 'pstool: load geqdsk ', trim(parm1), trim(parm2), bdy_crat
         call pstool_geqdsk(parm1,parm2,bdy_crat)

      case ('load_geqdsk')

         read(parm3,'(E15.7)') bdy_crat
         write(*,*) 'pstool: load geqdsk ', trim(parm1), trim(parm2), bdy_crat
         call pstool_load_geqdsk(parm1,parm2, bdy_crat)

      case ('load_nbi')

         write(*,*) 'pstool: load nbi'
         call pstool_load_nbi(parm1)

      case default

         write(*,*) 'default'

   end select

end program
