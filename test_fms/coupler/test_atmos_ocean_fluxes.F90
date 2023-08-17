program test_atmos_ocean_fluxes

  use fms_mod, only: fms_init
  use coupler_types_mod
  use atmos_ocean_fluxes_mod
  use field_manager_mod
  use fm_util_mod
  use mpp_mod, only: mpp_error, FATAL

  implicit none

  type(coupler_1d_bc_type) :: gas_fluxes
  type(coupler_1d_bc_type) :: gas_fields_atm
  type(coupler_1d_bc_type) :: gas_fields_ice

  call fms_init
  call test_aof_set_coupler_flux
  call test_atmos_ocean_fluxes_init
  call additional_test

contains
  !--------------------------------------
  subroutine test_atmos_ocean_fluxes_init

    implicit none

    call atmos_ocean_fluxes_init(gas_fluxes, gas_fields_atm, gas_fields_ice, verbosity=5)

  end subroutine test_atmos_ocean_fluxes_init
  !--------------------------------------
  subroutine test_aof_set_coupler_flux

    !> this test sets up the flux fields

    implicit none

    integer, parameter :: n=4

    integer, dimension(n) :: coupler_index, atm_tr_index, nparameter
    character(50), dimension(n) :: fname, flux_type, impl     !< field names

    character(100) :: cresults, thelist
    real :: rresults, rresults2(n)

    real, dimension(n) :: mol_wt, param_array
    integer :: i, success, nn

    fname=["vampires", "eric", "grapes", "coffee"]
    flux_type=["air_sea_gas_flux_generic", &
               "air_sea_gas_flux", &
               "air_sea_deposition", &
               "land_sea_runoff"]
    impl=["ocmip2", "linear", "dry", "river"]
    param_array =[1.0,2.0,3.0,4.0]
    atm_tr_index=[1, 2, 3, 4]
    mol_wt      =[1.0, 2.0, 3.0, 4.0]
    nparameter  =[2, 3, 1, 1]

    call atmos_ocean_type_fluxes_init
    !> set flux fields
    do i=1, n
       coupler_index(1)=aof_set_coupler_flux(name=fname(i), &
                                             flux_type=flux_type(i), &
                                             implementation=impl(i), &
                                             atm_tr_index=atm_tr_index(i), &
                                             param=param_array(1:nparameter(i)), &
                                             mol_wt=mol_wt(i))
    end do

    !> check answers
    do i=1, n
       !> check to see that the flux field has been created
       thelist="coupler_mod/fluxes/"//trim(fname(i))
       if(.not. fm_exists(trim(thelist))) call mpp_error(FATAL, 'coupler')
       !> check to see that the flux_type value has been set correctly
       success=fm_get_value(trim(thelist)//"/flux_type", cresults)
       if(trim(cresults).ne.trim(flux_type(i))) call mpp_error(FATAL,'test_aof_set_coupler_flux ERROR with flux_type')
       !> check to see the implementation value has been set correctly
       success=fm_get_value(trim(thelist)//"/implementation", cresults)
       if(trim(cresults).ne.trim(impl(i))) call mpp_error(FATAL,'test_aof_set_coupler_flux ERROR with implementation')
       !> check to see the param array has been set correctly
       nn=nparameter(i)
       rresults2(1:nn)=fm_util_get_real_array(trim(thelist)//"/param")
       if( any(rresults2(1:nn).ne.param_array(1:nn)) ) call mpp_error(FATAL, 'test_aof_set_coupler_fluxes ERROR with param')
       !> check to see that the mol_wt value has been set correctly
       success=fm_get_value(trim(thelist)//'/mol_wt', rresults)
       if( rresults.ne.mol_wt(i)) call mpp_error(FATAL, 'test_aof_set_coupler_fluxes ERROR with mol_wt')
    end do

  end subroutine test_aof_set_coupler_flux
  !--------------------------------------
  subroutine additional_test

    implicit none

    !write(*,*) 'hereherehere'
    !write(*,*) gas_fluxes%num_bcs
    !write(*,*) fm_util_get_length('/coupler_mod/fluxes/')

  end subroutine additional_test
  !--------------------------------------
end program test_atmos_ocean_fluxes
