!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!! @brief unit test for mpp_global_field
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests the subroutines for mpp_global_field

program test_mpp_global_field

  use platform_mod
  use mpp_mod,         only : mpp_init, mpp_error, FATAL, NOTE, mpp_init_test_read_namelist, mpp_exit
  use mpp_mod,         only : mpp_chksum, mpp_declare_pelist, mpp_pe, mpp_npes, mpp_root_pe, mpp_sync_self
  use mpp_domains_mod, only : domain2D !: variable
  use mpp_domains_mod, only : CENTER, EAST, NORTH, CORNER
  use mpp_domains_mod, only : XUPDATE, YUPDATE
  use mpp_domains_mod, only : mpp_domains_init, mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field

  implicit none

  integer :: pe, npes
  integer :: stdunit=6
  integer :: nx=128, ny=128, nz=40, stackmax=4000000
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: layout(2)

  !> call mpp_init
  call mpp_init(test_level=mpp_init_test_read_namelist)

  !> get pe number and total number of pes
  pe = mpp_pe()
  npes = mpp_npes()

  !--- initialize mpp domains
  call mpp_domains_init()
  call mpp_domains_set_stack_size(stackmax)

  !> part of test_interface
  call test_global_field( 'Non-symmetry' )
  call test_global_field( 'Symmetry center' )
  call test_global_field( 'Symmetry corner' )
  call test_global_field( 'Symmetry east' )
  call test_global_field( 'Symmetry north' )

contains

  subroutine test_global_field( type )

    implicit none

    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, gcheck
    type(domain2D) :: domain
    real, allocatable    :: global1(:,:,:)
    integer              :: ishift, jshift, ni, nj, i, j, position
    integer, allocatable :: pelist(:)
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    integer              :: k

    !--- set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )


    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- determine if an extra point is needed
    ishift = 0; jshift = 0
    position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1; jshift = 1; position=CORNER
    case ('Symmetry east')
       ishift = 1; jshift = 0; position=EAST
    case ('Symmetry north')
       ishift = 0; jshift = 1; position=NORTH
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-whalo:ni+ehalo, 1-shalo:nj+nhalo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    allocate( gcheck(ni, nj, nz) )
    allocate( x (isd:ied,jsd:jed,nz) )

    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !--- test the data on data domain
    gcheck = 0.

    !:MKL id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !:MKL call mpp_clock_end  (id)

    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on data domain' )

    !--- Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !--- will be declared. But for the x-direction global field, mpp_sync_self will
    !--- be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !--- in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !--- deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !--- will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !--- on all pe is needed for those partial pelist. But for y-update, it is ok.
    !--- because the pelist in y-update is not continous.

    allocate(pelist(0:layout(1)-1))
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !xupdate
    gcheck = 0.
    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = XUPDATE, position=position )
    !:MKL call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on data domain' )

    !yupdate
    gcheck = 0.
    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = YUPDATE, position=position )
    !:MKL call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on data domain' )

    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )

    !:MKL call mpp_clock_end(id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, &
                            type//' mpp_global_field on data domain' )


    !--- test the data on compute domain
    gcheck = 0.
    !:MKL id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je, :), gcheck, position=position )
    !:MKL call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on compute domain' )

    !xupdate
    gcheck = 0.
    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = XUPDATE, position=position )
    !:MKL call mpp_clock_end(id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on compute domain' )

    !yupdate
    gcheck = 0.
    !:MKL call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = YUPDATE, position=position )
    !:MKL call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field



  !: duplicate compare_checksums, also in mpp_domains
  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(i8_kind) :: sum1, sum2
    integer :: i, j, k

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    call mpp_sync_self()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_chksum: size of a and b does not match')

    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
                write(*,'(a,i3,a,i3,a,i3,a,i3,a,f20.9,a,f20.9)') trim(string)//" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
       if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       write(stdunit, *)"sum1 =", sum1, mpp_pe()
       write(stdunit, *)"sum2 =", sum2, mpp_pe()
       write(stdunit,'(a,i3,a,i20,a,i20)')" at pe ", mpp_pe(), " sum(a)=", sum1, " sum(b)=", sum2
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums


end program test_mpp_global_field
