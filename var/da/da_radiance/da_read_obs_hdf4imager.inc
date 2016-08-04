subroutine da_read_obs_hdf4imager(iv,infile_tb)

   !--------------------------------------------------------------------
   ! Purpose: read goes imager data in hdf4 from CIMSS
   !
   ! History: 2/5/2015                 Chun Yang
   !--------------------------------------------------------------------

   implicit none

   type (iv_type),   intent (inout) :: iv
   character(len=*), intent(in)     :: infile_tb

   type (datalink_type), pointer  :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   ! variable for hdf4
   integer, parameter  :: DFACC_READ=1
   integer             :: sfstart,sfn2index,sffattr,sfselect,sfginfo,sfrdata
   integer             :: sfendacc, sfend, sfrnatt
   integer             :: sd_id,sds_id,sds_index,status,attr_index
   integer             :: rank, data_type, n_attrs
   integer             :: dim_sizes(2),start(2),edges(2),stride(2)
   character(len=100)  :: sds_name

   integer(2), allocatable :: latitude_i(:,:)
   integer(2), allocatable :: longitude_i(:,:)
   integer(2), allocatable :: temp_3_75um_nom_i(:,:)
   integer(2), allocatable :: temp_6_7um_nom_i(:,:)
   integer(2), allocatable :: temp_11_0um_nom_i(:,:)
   integer(2), allocatable :: temp_13_3um_nom_i(:,:)
   integer(1), allocatable :: sensor_zenith_angle_i(:,:)
   integer(1), allocatable :: solar_zenith_angle_i(:,:)
   integer(1), allocatable :: sensor_azimuth_angle_i(:,:)
   integer(1), allocatable :: solar_azimuth_angle_i(:,:)

   integer(1), allocatable :: bad_pixel_mask(:,:)
   integer(1), allocatable :: cloud_mask(:,:)
   integer(1), allocatable :: surface_type(:,:)
   integer(2)              :: start_year
   integer(2)              :: start_day

   real(4)                 :: start_time
   real(4)                 :: scale_factor, add_offset
   real(4),    allocatable :: scan_time(:)

   real, allocatable  :: latitude(:,:)
   real, allocatable  :: longitude(:,:)
   real, allocatable  :: sensor_zenith_angle(:,:)
   real, allocatable  :: solar_zenith_angle(:,:)
   real, allocatable  :: solar_azimuth_angle(:,:)
   real, allocatable  :: sensor_azimuth_angle(:,:)
   real, allocatable  :: temp_3_75um_nom(:,:)
   real, allocatable  :: temp_6_7um_nom(:,:)
   real, allocatable  :: temp_11_0um_nom(:,:)
   real, allocatable  :: temp_13_3um_nom(:,:)

   real, parameter :: r360   = 360.0_r_kind
   real, parameter :: tbmin  = 180._r_kind
   real, parameter :: tbmax  = 340._r_kind

   real          :: obs_time
   integer       :: idate5(5)
   integer       :: idate, im, iy, idd, ihh
   integer       :: inst,platform_id,satellite_id,sensor_id
   real          :: tb, crit
   integer       :: ifgat, iout, iobs
   logical       :: outside, outside_all, iuse

   integer       :: i,j,k,l,m,n, numhdf,ihdf,ngoes,landsea_mask
   logical       :: found, head_found, head_allocated
   character(10) :: date
   real          :: dlon_earth,dlat_earth,dlon_earth_deg,dlat_earth_deg

   integer  :: num_imager_local, num_imager_global, num_imager_used, num_imager_thinned
   integer  :: num_imager_used_tmp, num_imager_file
   integer  :: num_imager_local_local, num_imager_global_local, num_imager_file_local

   integer  :: itx, itt

   character(80)    :: filename
   logical          :: fexist
   integer          :: nchan, ifov,iscan,ichannels

  ! Allocatable arrays
   integer,    allocatable     :: ptotal(:),nread(:)
   real,               allocatable     :: in(:), out(:)
   real,       allocatable     :: data_all(:)
   integer,    allocatable     :: cloud_mask_all(:)

   integer, allocatable   :: obstime(:,:)

   integer :: doy,mlen(12),mday(12),mon

   real     :: theta_true
   integer  :: ios, goes_unit
   data  mlen/31,28,31,30,31,30, 31,31,30,31,30,31/

!-----------------------------------------------------------------------

   if ( trace_use ) call da_trace_entry("da_read_obs_hdf4imager")

!---------------------------------------

   platform_id    = 4
   sensor_id      = 22

   if(infile_tb(12:13)=='13') then
      satellite_id = 13
   elseif(infile_tb(12:13)=='15') then
      satellite_id = 15
   else
      write(*,*) 'This character is not coped with presently'
      return
   endif

   inst = 0
   do ngoes = 1, rtminit_nsensor
      if (platform_id  == rtminit_platform(ngoes) &
         .and. sensor_id == rtminit_sensor(ngoes) &
         .and. satellite_id == rtminit_satid(ngoes)) then
            inst = ngoes
      else
         cycle
      end if
   end do 
   if (inst == 0) then
      write(unit=message(1),fmt='(A,I2)') "No observation : goes-",satellite_id
      call da_warning(__FILE__,__LINE__, message(1:1))
      call da_trace_exit("da_read_obs_hdf4imager")
      return
   end if

   numhdf=0
   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0
   num_imager_file     = 0
   num_imager_local    = 0
   num_imager_global   = 0
   num_imager_used     = 0
   num_imager_thinned  = 0

   nchan = iv%instid(inst)%nchan  !4          !fellow the infomation file
   write(unit=stdout,fmt=*)'GOES-IMAGER nchan: ',nchan
   allocate(data_all(1:nchan))
   allocate(cloud_mask_all(1:nchan))

  !-------------------------------------------------------------------------
  ! 1.0 Assign file names and prepare to read imager files
  !-------------------------------------------------------------------------
   filename = trim(infile_tb)//'.hdf'
   inquire (file=filename, exist=fexist)
   if ( fexist )  numhdf=1

   if ( numhdf == 0 ) then
      call da_warning(__FILE__,__LINE__, &
         (/"No valid GOES HDF4 file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf4imager")
      return
   end if

      num_imager_file_local    = 0
      num_imager_local_local   = 0
      num_imager_global_local  = 0

      sd_id = sfstart(filename, DFACC_READ)
      if(sd_id.lt.0)then
         call da_warning(__FILE__,__LINE__, &
            (/"Cannot open HDF4 file "//trim(filename)/))
         return
      endif

!-------------- Read bad_pixel_mask ------------------
      sds_index = sfn2index(sd_id,"bad_pixel_mask")
      sds_id = sfselect(sd_id,sds_index)
      status = sfginfo(sds_id,sds_name,rank,dim_sizes,data_type,n_attrs)

      start(1)=0
      start(2)=0
      edges(1)=dim_sizes(1)
      edges(2)=dim_sizes(2)
      stride(1)=1
      stride(2)=1
      allocate(bad_pixel_mask(dim_sizes(1),dim_sizes(2)))
      status = sfrdata(sds_id,start,stride,edges,bad_pixel_mask)
!------------- Read cloud_mask --------------------
      allocate(cloud_mask(dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"cloud_mask")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,cloud_mask)
!------------- Read surface_type --------------------
      allocate(surface_type(dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"surface_type")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,surface_type)
!------------- Read scan_line_time --------------------
      allocate(scan_time(dim_sizes(2)))
      sds_index = sfn2index(sd_id,"scan_line_time")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,dim_sizes(2),scan_time)
!---------------- Read start_year/day/time from hdf-----------------
      allocate  (obstime(5,dim_sizes(2))) 
      attr_index = sffattr(sd_id,'START_YEAR')
      status = sfrnatt(sd_id,attr_index,start_year)
      attr_index = sffattr(sd_id,'START_DAY')
      status = sfrnatt(sd_id,attr_index,start_day)

      IF ( MOD(start_year,4) == 0 ) THEN
         mlen(2) = 29
         IF ( MOD(start_year,100) == 0 ) THEN
            mlen(2) = 28
         ENDIF
         IF ( MOD(start_year,400) == 0 ) THEN
            mlen(2) = 29
         ENDIF
      ENDIF
      obstime(1,:) = start_year
      m = 0
      do mon=1,12
         mday(mon) = m
         m = m + mlen(mon)
      end do  
      do mon= 2,12
         if (start_day.gt.mday(mon-1).and.start_day.le.mday(mon))then
            obstime(2,:) = mon-1
            obstime(3,:) = start_day - mday(mon-1)
         end if
      end do
      do i = 1,dim_sizes(2)
         obstime(4,i) = int(scan_time(i))
         obstime(5,i) = int((scan_time(i)-int(scan_time(i)))*60.0)
      end do
!------------------latitude--------------------------
      allocate(latitude_i(dim_sizes(1),dim_sizes(2)))
      allocate(latitude  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"latitude")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,latitude_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            latitude(i,j) = latitude_i(i,j)
            if(latitude_i(i,j).gt.-32768) latitude(i,j) =latitude_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(latitude_i)
!------------------longitude--------------------------
      allocate(longitude_i(dim_sizes(1),dim_sizes(2)))
      allocate(longitude  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"longitude")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,longitude_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            longitude(i,j) = longitude_i(i,j)
            if(longitude_i(i,j).gt.-32768) longitude(i,j) =longitude_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(longitude_i)
!------------------temp_3_75um_nom--------------------------
      allocate(temp_3_75um_nom_i(dim_sizes(1),dim_sizes(2)))
      allocate(temp_3_75um_nom  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"temp_3_75um_nom")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,temp_3_75um_nom_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            temp_3_75um_nom(i,j) = temp_3_75um_nom_i(i,j)
            if(temp_3_75um_nom_i(i,j).gt.-32768) temp_3_75um_nom(i,j) =temp_3_75um_nom_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(temp_3_75um_nom_i)
!------------------temp_6_7um_nom--------------------------
      allocate(temp_6_7um_nom_i(dim_sizes(1),dim_sizes(2)))
      allocate(temp_6_7um_nom  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"temp_6_7um_nom")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,temp_6_7um_nom_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            temp_6_7um_nom(i,j) = temp_6_7um_nom_i(i,j)
            if(temp_6_7um_nom_i(i,j).gt.-32768) temp_6_7um_nom(i,j) =temp_6_7um_nom_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(temp_6_7um_nom_i)
!------------------temp_11_0um_nom--------------------------
      allocate(temp_11_0um_nom_i(dim_sizes(1),dim_sizes(2)))
      allocate(temp_11_0um_nom  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"temp_11_0um_nom")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,temp_11_0um_nom_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            temp_11_0um_nom(i,j) = temp_11_0um_nom_i(i,j)
            if(temp_11_0um_nom_i(i,j).gt.-32768)temp_11_0um_nom(i,j) =temp_11_0um_nom_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(temp_11_0um_nom_i)
!------------------temp_13_3um_nom--------------------------
      allocate(temp_13_3um_nom_i(dim_sizes(1),dim_sizes(2)))
      allocate(temp_13_3um_nom  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"temp_13_3um_nom")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,temp_13_3um_nom_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            temp_13_3um_nom(i,j) = temp_13_3um_nom_i(i,j)
            if(temp_13_3um_nom_i(i,j).gt.-32768)temp_13_3um_nom(i,j) =temp_13_3um_nom_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(temp_13_3um_nom_i)
!------------------sensor_zenith_angle--------------------------
      allocate(sensor_zenith_angle_i(dim_sizes(1),dim_sizes(2)))
      allocate(sensor_zenith_angle  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"sensor_zenith_angle")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,sensor_zenith_angle_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            sensor_zenith_angle(i,j) = sensor_zenith_angle_i(i,j)
            if(sensor_zenith_angle_i(i,j).gt.-128) sensor_zenith_angle(i,j) =sensor_zenith_angle_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(sensor_zenith_angle_i)
!------------------solar_zenith_angle--------------------------
      allocate(solar_zenith_angle_i(dim_sizes(1),dim_sizes(2)))
      allocate(solar_zenith_angle  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"solar_zenith_angle")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,solar_zenith_angle_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            solar_zenith_angle(i,j) = solar_zenith_angle_i(i,j)
            if(solar_zenith_angle_i(i,j).gt.-128) solar_zenith_angle(i,j) =solar_zenith_angle_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(solar_zenith_angle_i)
!------------------solar_azimuth_angle--------------------------
      allocate(solar_azimuth_angle_i(dim_sizes(1),dim_sizes(2)))
      allocate(solar_azimuth_angle  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"solar_azimuth_angle")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,solar_azimuth_angle_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            solar_azimuth_angle(i,j) = solar_azimuth_angle_i(i,j)
            if(solar_azimuth_angle_i(i,j).gt.-128) solar_azimuth_angle(i,j) =solar_azimuth_angle_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(solar_azimuth_angle_i)
!------------------sensor_azimuth_angle--------------------------
      allocate(sensor_azimuth_angle_i(dim_sizes(1),dim_sizes(2)))
      allocate(sensor_azimuth_angle  (dim_sizes(1),dim_sizes(2)))
      sds_index = sfn2index(sd_id,"sensor_azimuth_angle")
      sds_id = sfselect(sd_id,sds_index)
      status = sfrdata(sds_id,start,stride,edges,sensor_azimuth_angle_i)
      attr_index = sffattr(sds_id,'scale_factor')
      status = sfrnatt(sds_id,attr_index,scale_factor)
      attr_index = sffattr(sds_id,'add_offset')
      status = sfrnatt(sds_id,attr_index,add_offset)
      do j = 1, dim_sizes(2)
         do i = 1, dim_sizes(1)
            sensor_azimuth_angle(i,j) = sensor_azimuth_angle_i(i,j)
            if(sensor_azimuth_angle_i(i,j).gt.-128) sensor_azimuth_angle(i,j) =sensor_azimuth_angle_i(i,j)*scale_factor + add_offset
         end do
      end do
      deallocate(sensor_azimuth_angle_i)

      status = sfendacc(sds_id)
      status = sfend(sd_id)


    !2.0 Loop to read hdf file and assign information to a sequential structure
    !------------------------------------------------------------------------

! Allocate arrays to hold data
      allocate(head)
      nullify (head % next )
      p => head

! ------------ start scan_loop ------------------
      scan_loop:     do iscan=1, dim_sizes(2)
         do i=1,5
            idate5(i) = obstime(i,iscan)
         end do
         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)
         if ( obs_time < time_slots(0) .or.  &
           obs_time >= time_slots(num_fgat_time) ) cycle scan_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and.  &
               obs_time  < time_slots(ifgat) ) exit
         end do

!-------------------- start fov_loop -------------------
         fov_loop:  do ifov=1,dim_sizes(1)
            num_imager_file = num_imager_file +1
            num_imager_file_local = num_imager_file_local + 1
            info%lat = latitude(ifov,iscan)
            info%lon = longitude(ifov,iscan)

            call da_llxy (info,loc,outside,outside_all)
            if (outside_all) cycle fov_loop

            num_imager_global = num_imager_global + 1
            num_imager_global_local = num_imager_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop

            if (sensor_zenith_angle(ifov,iscan) > 80.0) cycle fov_loop !YC

            num_imager_local = num_imager_local + 1
            num_imager_local_local = num_imager_local_local + 1

            if (surface_type(ifov,iscan) > 0) then
               landsea_mask = 0
            else if (surface_type(ifov,iscan) == 0) then
               landsea_mask = 1
            end if

            write(unit=info%date_char, &
               fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')  &
               idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
               ':', idate5(5),':00'

            info%elv = 0.0

         ! 3.0 Make Thinning
         ! Map obs to thinning grid
         !--------------------------------------------

            if (thinning) then
               dlat_earth = info%lat
               dlon_earth = info%lon
               if (dlon_earth< zero) dlon_earth = dlon_earth + r360
               if (dlon_earth>=r360) dlon_earth = dlon_earth - r360
               dlat_earth_deg = dlat_earth
               dlon_earth_deg = dlon_earth

               dlat_earth = dlat_earth*deg2rad
               dlon_earth = dlon_earth*deg2rad

               crit = 1.
               call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
               if (.not. iuse) then
                  num_imager_thinned = num_imager_thinned+1
                  cycle fov_loop
               end if
            end if

            num_imager_used = num_imager_used + 1

            data_all = missing_r
            tb = temp_3_75um_nom(ifov,iscan)
            if (tb < tbmin .or. tb > tbmax) tb = missing_r
            data_all(1) = tb

            tb = temp_6_7um_nom(ifov,iscan)
            if (tb < tbmin .or. tb > tbmax) tb = missing_r
            data_all(2) = tb 

            tb = temp_11_0um_nom(ifov,iscan)
            if (tb < tbmin .or. tb > tbmax) tb = missing_r
            data_all(3) = tb

            tb = temp_13_3um_nom(ifov,iscan)
            if (tb < tbmin .or. tb > tbmax) tb = missing_r
            data_all(4) = tb

            cloud_mask_all = missing_r
            cloud_mask_all(1) = cloud_mask(ifov,iscan)
            cloud_mask_all(2) = cloud_mask(ifov,iscan)
            cloud_mask_all(3) = cloud_mask(ifov,iscan)
            cloud_mask_all(4) = cloud_mask(ifov,iscan)

     ! 4.0 assign information to sequential radiance structure
     !--------------------------------------------------------------------------
            allocate ( p % tb_inv (1:nchan ))
            allocate ( p % cloud_flag(1:nchan))
            p%info             = info
            p%loc              = loc
            p%landsea_mask     = landsea_mask
            p%scanpos          = ifov
            p%cloud_flag(1:nchan) = cloud_mask_all(1:nchan)
            p%satzen           = sensor_zenith_angle(ifov,iscan)
            p%satazi           = sensor_azimuth_angle(ifov,iscan)
            p%solzen           = solar_zenith_angle(ifov,iscan)
            p%solazi           = solar_azimuth_angle(ifov,iscan)
            p%tb_inv(1:nchan)  = data_all(1:nchan)
            p%sensor_index     = inst
            p%ifgat            = ifgat

            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)
         end do fov_loop
      end do scan_loop

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(filename),' got num_imager_file    : ',num_imager_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(filename),' got num_imager_global  : ',num_imager_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(filename),' got num_imager_local   : ',num_imager_local_local

   deallocate(data_all) ! Deallocate data arrays
   deallocate(cloud_mask_all)


   if (thinning .and. num_imager_global > 0 ) then
#ifdef DM_PARALLEL 
   ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         j = j + thinning_grid(inst,ifgat)%itxmax
      end do 

      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            in(j) = thinning_grid(inst,ifgat)%score_crit(i)
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            if ( ABS(out(j)-thinning_grid(inst,ifgat)%score_crit(i)) > 1.0E-10 ) &
               thinning_grid(inst,ifgat)%ibest_obs(i) = 0
         end do
      end do

      deallocate( in  )
      deallocate( out )
#endif

   ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_imager_used_tmp = num_imager_used
      do j = 1, num_imager_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current % cloud_flag ) !yc
            deallocate ( current )
            num_imager_thinned = num_imager_thinned + 1
            num_imager_used = num_imager_used - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_imager_used
   iv%total_rad_channel = iv%total_rad_channel + num_imager_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_imager_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_imager_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning(__FILE__,__LINE__,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'GOES-IMAGER data counts: '
   write(stdout,fmt='(a,i7)') ' In file: ',num_imager_file
   write(stdout,fmt='(a,i7)') ' Global : ',num_imager_global
   write(stdout,fmt='(a,i7)') ' Local  : ',num_imager_local
   write(stdout,fmt='(a,i7)') ' Used   : ',num_imager_used
   write(stdout,fmt='(a,i7)') ' Thinned: ',num_imager_thinned


   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------

   if (num_imager_used > 0) then
      iv%instid(inst)%num_rad  = num_imager_used
      iv%instid(inst)%info%nlocal = num_imager_used
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
         'Allocating space for radiance innov structure', &
          inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad
      call da_allocate_rad_iv (inst, nchan, iv)
   end if


   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   p => head

   do n = 1, num_imager_used
      i = p%sensor_index
      call da_initialize_rad_iv (i, n, iv, p)

      current => p
      p => p%next
      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current % cloud_flag ) !yc
      deallocate ( current )
   end do

   deallocate ( p )
   deallocate (ptotal)

! Dellocate arrays
   deallocate  (obstime)
   deallocate  (bad_pixel_mask)
   deallocate  (cloud_mask)
   deallocate  (scan_time)
   deallocate  (longitude)
   deallocate  (latitude)
   deallocate  (surface_type)
   deallocate  (temp_13_3um_nom)
   deallocate  (temp_11_0um_nom)
   deallocate  (temp_6_7um_nom)
   deallocate  (temp_3_75um_nom)
   deallocate  (sensor_zenith_angle)
   deallocate  (sensor_azimuth_angle)
   deallocate  (solar_zenith_angle)
   deallocate  (solar_azimuth_angle)

   if ( trace_use ) call da_trace_exit("da_read_obs_hdf4imager")


end subroutine da_read_obs_hdf4imager

