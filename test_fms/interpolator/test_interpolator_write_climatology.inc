subroutine write_climatology_file

  use constants_mod, only : PI

  implicit none

  type(FmsNetcdfFile_t) :: fileobj

  integer, parameter :: ntime=3
  integer, parameter :: npfull=3, nphalf=npfull+1

  integer :: i,j,k,l
  integer, allocatable :: pes(:)
  integer :: stringlen
  character(100) :: mystring

  real, dimension(ntime)  :: timedata
  real, dimension(nlonlat)  :: latdata, londata
  real, dimension(nlonlatb) :: latbdata, lonbdata
  real, dimension(npfull) :: pfulldata
  real, dimension(nphalf) :: phalfdata
  real, dimension(nlonlat,nlonlat,npfull,ntime) :: ozonedata

  real :: dx, dy

  allocate(pes(mpp_npes()))
  call mpp_get_current_pelist(pes)

  ! assign time
  timedata=[1, 2, 3]

  !assign latdata
  dx = 1.0
  dy = 1.0

  do i=1, nlonlat
     latdata(i)= real(2*(i-1))
     londata(i)= real(2*(i-1))
  end do

  lonbdata(1)=0.0
  latbdata(1)=0.0
  do i=2,nlonlat
     lonbdata(i)=0.5*(londata(i-1)+londata(i))
     latbdata(i)=0.5*(latdata(i-1)+latdata(i))
  enddo
  lonbdata(nlonlat+1)=real(2*(nlonlat-1))
  latbdata(nlonlat+1)=real(2*(nlonlat-1))

  !p data
  pfulldata=[0.0001, 0.0003, 0.0005]
  phalfdata=[0.0, 0.0002, 0.0004, 0.0005]

  !ozone data
  do i=1, ntime
     do j=1, npfull
        do k=1, nlonlat
           do l=1, nlonlat
              ozonedata(l,k,j,i)=real(j)
           end do
        end do
     end do
  end do


  ! write netcdf file
  if(open_file(fileobj, 'INPUT/'//trim(ncfile), 'overwrite', pelist=pes)) then
     call register_axis(fileobj, 'time',unlimited)
     call register_axis(fileobj, 'lon',  nlonlat)
     call register_axis(fileobj, 'lat',  nlonlat)
     call register_axis(fileobj, 'latb', nlonlatb)
     call register_axis(fileobj, 'lonb', nlonlatb)
     call register_axis(fileobj, 'pfull', npfull)
     call register_axis(fileobj, 'phalf', nphalf)
     call register_axis(fileobj, 'longitude', nlonlat)
     call register_axis(fileobj, 'latitude',  nlonlat)

     !time axes
     call register_field(fileobj, 'time', 'double', dimensions=(/'time'/))
     call register_string_attribute(fileobj, 'time', 'axis', 'T')
     call register_string_attribute(fileobj, 'time', 'calendar', 'noleap')
     call register_string_attribute(fileobj, 'time', 'climatology', '1979-01-01 00:00:00, 1998-01-01 00:00:00')
     call register_string_attribute(fileobj, 'time', 'units', 'days since 1849-01-01 00:00:00')

     !lon axes
     call register_field(fileobj, 'lon', 'float', dimensions=(/'longitude'/))
     call register_string_attribute(fileobj, 'lon', 'edges', 'lonb')
     call register_string_attribute(fileobj, 'lon', 'cartesian_axis', 'X')
     call register_string_attribute(fileobj, 'lon', 'units', 'degrees_east')

     !lonb  axes longitude edges
     call register_field(fileobj, 'lonb', 'float', dimensions=(/'lonb'/))
     call register_string_attribute(fileobj, 'lonb', 'cartesian_axis','X')
     call register_string_attribute(fileobj, 'lonb', 'units', 'degrees_edges')

     !lat
     call register_field(fileobj, 'lat', 'float', dimensions=(/'latitude'/))
     call register_string_attribute(fileobj, 'lat', 'edges', 'latb')
     call register_string_attribute(fileobj, 'lat', 'cartesian_axis', 'Y')
     call register_string_attribute(fileobj, 'lat', 'units', 'degrees_north')

     !latb latitude_edges
     call register_field(fileobj, 'latb', 'float', dimensions=(/'latb'/))
     call register_string_attribute(fileobj, 'latb', 'cartesian_axis', 'Y')
     call register_string_attribute(fileobj, 'latb', 'units', 'degrees_N')

     !pfull layer-mean pressure
     call register_field(fileobj, 'pfull', 'float', dimensions=(/'pfull'/))
     call register_string_attribute(fileobj, 'pfull', 'positive', 'down')
     call register_string_attribute(fileobj, 'pfull', 'cartesian_axis', 'Z')
     call register_string_attribute(fileobj, 'pfull', 'units', 'hPa')

     !phalf half pressure level
     call register_field(fileobj, 'phalf', 'float', dimensions=(/'phalf'/))
     call register_string_attribute(fileobj, 'phalf', 'positive', 'down')
     call register_string_attribute(fileobj, 'phalf', 'cartesian_axis', 'Z')
     call register_string_attribute(fileobj, 'phalf', 'units', 'hPa')

     !ozone ozone mass mixing ratio
     call register_field(fileobj, 'ozone', 'float', dimensions=(/'lon','lat','pfull','time'/))
     call register_string_attribute(fileobj, 'ozone', 'units', 'kg/kg')

     call write_data(fileobj, 'time', timedata)
     call write_data(fileobj, 'lon', londata)
     call write_data(fileobj, 'lonb', lonbdata)
     call write_data(fileobj, 'lat', latdata)
     call write_data(fileobj, 'latb', latbdata)
     call write_data(fileobj, 'pfull', pfulldata)
     call write_data(fileobj, 'phalf', phalfdata)
     call write_data(fileobj, 'ozone', ozonedata)

     call close_file(fileobj)
  end if

end subroutine write_climatology_file

subroutine register_string_attribute(fileobj, variable, attribute, att_value)

  implicit none

  type(FmsNetcdfFile_t) :: fileobj
  character(*) :: variable, attribute, att_value
  integer :: str_len

  str_len=len(trim(att_value))
  call register_variable_attribute(fileobj, trim(variable), trim(attribute), trim(att_value), str_len=str_len)

end subroutine register_string_attribute
