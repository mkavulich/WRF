program da_update_bc

   !-----------------------------------------------------------------------
   ! Purpose: update BC file from wrfvar output.
   ! current version reads only wrf-netcdf file format
   !
   ! Jliu, 08/10/2011:
   !   Introducing 4dvar namelist variables: var4d, multi_inc, da_bdy_file
   !
   !      var4d     = .false.          updating bdy with da_file, wrf_bdy_file,
   !                                   wrf_input, default.
   !                                   true, updating bdy with 4dvar output bdy
   !      multi_inc = .false.          no mulit_inc, default.
   !                                   ture for 4dvar with multi_inc
   !
   !      da_bdy_file = wrfvar_bdyout  default 4dvar output bdy file
   !      
   !   Only updating bdy with wrfvar_bdyout when 
   !         1) var4d     = .true.
   !         2) multi_inc = .false.
   !   Otherwise, updating bdy with da_file, wrf_bdy_file, wrf_input
   !
   ! Y.-R. Guo, 03/18/2008:
   !   1) Fixed the bug for low_bdy_only;
   !   2) Introducing another namelist variable: update_lsm
   !      update_lsm = .true. --- The LSM predicted variables: 
   !                         TSLB, SMOIS, SNOW, SH2O, RHOSN, CANWAT, SNOWH
   !                              will be updated based on wrf_input file
   !                 = .false. -- no updated, default.
   !
   !-----------------------------------------------------------------------

   use da_netcdf_interface, only : da_get_var_3d_real_cdf, &
      da_put_var_3d_real_cdf, da_get_dims_cdf, da_put_var_2d_real_cdf, &
      da_get_var_2d_real_cdf, da_get_var_2d_int_cdf, da_get_bdytimestr_cdf, &
      da_get_times_cdf, da_get_bdyfrq, stderr, stdout, da_put_var_2d_int_cdf

   use da_module_couple_uv, only : da_couple_uv

   implicit none

   include 'netcdf.inc'

   integer, parameter :: max_3d_variables = 20, &
                         max_2d_variables = 25
 
   character(len=512) :: da_file,      &
                         wrf_bdy_file, &
                         wrf_input, da_bdy_file
 
   character(len=20) :: var_pref, var_name, vbt_name

   character(len=20) :: var3d(max_3d_variables), &
                        varsf(max_2d_variables)

   character(len=10), dimension(4) :: bdyname, tenname

   integer           :: ids, ide, jds, jde, kds, kde
   integer           :: num3d, num2d, ndims
   integer           :: time_level, time_id
   integer           :: i,j,k,l,m,n

   integer, dimension(4) :: dims
 
   real, allocatable, dimension(:,:,:) :: tend3d, scnd3d, frst3d, full3d

   real, allocatable, dimension(:,:,:) :: u, v

   real, allocatable, dimension(:,  :) :: mu, mub, msfu, msfv, msfm, &
                                          tend2d, scnd2d, frst2d, full2d

   real, allocatable, dimension(:,  :) :: tsk, tsk_wrfvar
   real, allocatable, dimension(:,:)   :: snow, snowc, snowh

   integer, allocatable, dimension(:,:) :: ivgtyp, full2dint

   character(len=80), allocatable, dimension(:) :: times, &
                                                   thisbdytime, nextbdytime
 
   integer :: east_end, north_end, io_status, cdfid, varid, domain_id, iswater
   integer :: iostatus(4)

   logical :: debug, update_lateral_bdy, update_low_bdy, update_lsm, keep_tsk_wrf
   logical :: keep_snow_wrf, var4d, multi_inc

   real :: bdyfrq

   character(len=512) :: wrfvar_output_file    ! obsolete. Kept for backward compatibility
   logical            :: cycling, low_bdy_only ! obsolete. Kept for backward compatibility

   integer, parameter :: namelist_unit = 7, &
                         ori_unit = 11, &
                         new_unit = 12

   namelist /control_param/ da_file,      &
                            wrf_bdy_file, &
                            wrf_input, domain_id, &
                            debug, update_lateral_bdy, update_low_bdy, update_lsm, &
                            keep_tsk_wrf, keep_snow_wrf, iswater, &
                            wrfvar_output_file, cycling, low_bdy_only, &
                            var4d, multi_inc, da_bdy_file

   da_bdy_file        = 'wrfvar_bdyout'
   da_file            = 'wrfvar_output'
   wrf_bdy_file       = 'wrfbdy_d01'
   wrf_input          = 'wrfinput_d01'
   domain_id          = 1

   debug              = .false. 
   update_lateral_bdy = .true.
   update_low_bdy     = .true.
   update_lsm         = .false.
   keep_tsk_wrf       = .true.
   keep_snow_wrf      = .true.
   iswater            = 16      ! USGS water index: 16, MODIS water index: 17

   var4d              = .false.
   multi_inc          = .false.

   time_id            = 1

   wrfvar_output_file = 'OBSOLETE'
   cycling            = .false.
   low_bdy_only       = .false.

   !---------------------------------------------------------------------
   ! Read namelist
   !---------------------------------------------------------------------
   io_status = 0

   open(unit = namelist_unit, file = 'parame.in', &
          status = 'old' , access = 'sequential', &
          form   = 'formatted', action = 'read', &
          iostat = io_status)

   if (io_status /= 0) then
      write(unit=stdout,fmt=*) 'Error to open namelist file: parame.in.'
      write(unit=stdout,fmt=*) 'Will work for updating lateral boundary only.'
   else
      read(unit=namelist_unit, nml = control_param , iostat = io_status)

      if (io_status /= 0) then
         write(unit=stdout,fmt=*) 'Error to read control_param. Stopped.'
         stop
      end if

      ! deal with the old namelist
      if ( index(wrfvar_output_file, 'OBSOLETE') <= 0 ) then
         ! wrfvar_output_file is set in the user's parame.in
         ! reset the settings
         da_file = wrfvar_output_file
         if ( domain_id > 1 ) then
            low_bdy_only = .true.
         end if
         if ( cycling .and. domain_id == 1 ) then
            update_lateral_bdy = .true.
            update_low_bdy     = .true.
         else
            if ( low_bdy_only ) then
               update_lateral_bdy = .false.
               update_low_bdy     = .true.
            else
               update_lateral_bdy = .true.
               update_low_bdy     = .false.
            end if
         end if
      end if

      if ( var4d ) then
        time_id = 2
        if ( multi_inc ) then
          var4d     = .false.
          time_id = 1
        end if
      end if
        
      WRITE(unit=stdout, fmt='(2a)') &
           'da_bdy_file   = ', trim(da_bdy_file), &
           'da_file       = ', trim(da_file), &
           'wrf_bdy_file  = ', trim(wrf_bdy_file), &
           'wrf_input     = ', trim(wrf_input)

      WRITE(unit=stdout, fmt='(a, L10)')             &
           'update_lateral_bdy = ', update_lateral_bdy, &
           'update_low_bdy     = ', update_low_bdy, &
           'var4d              = ', var4d, &
           'multi_inc          = ', multi_inc

      if ( update_lsm ) keep_snow_wrf = .false.

      close(unit=namelist_unit)
   end if

   ! 3D need update
   num3d=6
   var3d(1)='U'
   var3d(2)='V'
   var3d(3)='W'
   var3d(4)='T'
   var3d(5)='PH'
   var3d(6)='QVAPOR'

   ! 2D need update
   num2d=23
   varsf(1)='MUB'
   varsf(2)='MU'
   varsf(3)='MAPFAC_U'
   varsf(4)='MAPFAC_V'
   varsf(5)='MAPFAC_M'
   varsf(6)='TMN'
   varsf(7)='SST'
   varsf(8)='TSK'
   varsf(9)='VEGFRA'
   varsf(10)='ALBBCK'
   varsf(11)='TSLB'
   varsf(12)='SMOIS'
   varsf(13)='SNOW'
   varsf(14)='SEAICE'
   varsf(15)='SH2O'
   varsf(16)='CANWAT'
   varsf(17)='RHOSN'
   varsf(18)='SNOWH'
   varsf(19)='LANDMASK'
   varsf(20)='IVGTYP'
   varsf(21)='ISLTYP'
   varsf(22)='SNOWC'
   varsf(23)='XLAND'

   if ( domain_id > 1 ) then
      write(unit=stdout, fmt='(a,i2)') 'Nested domain ID=',domain_id
      write(unit=stdout, fmt='(a)') &
        'No wrfbdy file needed, only low boundary need to be updated.'
      if ( update_lateral_bdy ) then
         write(unit=stdout, fmt='(a)') &
            'Re-setting update_lateral_bdy to be false for nested domain.'
         update_lateral_bdy = .false.
      end if
      update_low_bdy     = .true.
   end if

   if ( update_lateral_bdy ) then
   ! First, the boundary times
   call da_get_dims_cdf(wrf_bdy_file, 'Times', dims, ndims, debug)

   if (debug) then
      write(unit=stdout, fmt='(a,i2,2x,a,4i6)') &
           'Times: ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
   end if

   time_level = dims(2)

   if (time_level < 1) then
      write(unit=stdout, fmt='(a,i2/a)') &
           'time_level = ', time_level, &
           'We need at least one time-level BDY.'
      stop 'Wrong BDY file.'
   end if

   allocate(times(dims(2)))
   allocate(thisbdytime(dims(2)))
   allocate(nextbdytime(dims(2)))

   call da_get_times_cdf(wrf_bdy_file, times, dims(2), dims(2), debug)

   call da_get_bdytimestr_cdf(wrf_bdy_file, 'thisbdytime', thisbdytime, dims(2), debug)
   call da_get_bdytimestr_cdf(wrf_bdy_file, 'nextbdytime', nextbdytime, dims(2), debug)

   call da_get_bdyfrq(thisbdytime(1), nextbdytime(1), bdyfrq, debug)

   if (debug) then
      do n=1, dims(2)
         write(unit=stdout, fmt='(3(a, i2, 2a,2x))') &
           '       times(', n, ')=', trim(times(n)), &
           'thisbdytime (', n, ')=', trim(thisbdytime(n)), &
           'nextbdytime (', n, ')=', trim(nextbdytime(n))
      end do
   end if

   end if

   east_end=0
   north_end=0

   cdfid = ncopn(da_file, NCNOWRIT, io_status )

   ! For 2D variables
   ! Get mu, mub, msfu, and msfv

   do n=1,num2d

      io_status = nf_inq_varid(cdfid, trim(varsf(n)), varid)
      if (io_status /= 0 ) then
         print '(/"N=",i2," io_status=",i5,5x,"VAR=",a,a)', &
                   n, io_status, trim(varsf(n)), " does not exist"
         cycle
      endif

      call da_get_dims_cdf( da_file, trim(varsf(n)), dims, &
         ndims, debug)

      select case(trim(varsf(n)))
      case ('MU') ;
         if ( .not. update_lateral_bdy .or. var4d ) cycle
         allocate(mu(dims(1), dims(2)))

         call da_get_var_2d_real_cdf( da_file, &
            trim(varsf(n)), mu, dims(1), dims(2), 1, debug)

         east_end=dims(1)+1
         north_end=dims(2)+1
      case ('MUB') ;
         if ( .not. update_lateral_bdy .or. var4d ) cycle
         allocate(mub(dims(1), dims(2)))

         call da_get_var_2d_real_cdf( da_file, trim(varsf(n)), mub, &
                                   dims(1), dims(2), 1, debug)
      case ('MAPFAC_U') ;
         if ( .not. update_lateral_bdy .or. var4d ) cycle
         allocate(msfu(dims(1), dims(2)))

         call da_get_var_2d_real_cdf( da_file, trim(varsf(n)), msfu, &
                                   dims(1), dims(2), 1, debug)
      case ('MAPFAC_V') ;
         if ( .not. update_lateral_bdy .or. var4d ) cycle
         allocate(msfv(dims(1), dims(2)))

         call da_get_var_2d_real_cdf( da_file, trim(varsf(n)), msfv, &
                                   dims(1), dims(2), 1, debug)
      case ('MAPFAC_M') ;
         if ( .not. update_lateral_bdy .or. var4d ) cycle
         allocate(msfm(dims(1), dims(2)))

         call da_get_var_2d_real_cdf( da_file, trim(varsf(n)), msfm, &
                                   dims(1), dims(2), 1, debug)
      case ('TSK') ;
         if ( .not. update_low_bdy ) cycle

         allocate(tsk(dims(1), dims(2)))
         allocate(tsk_wrfvar(dims(1), dims(2)))
         allocate(ivgtyp(dims(1), dims(2)))

         call da_get_var_2d_real_cdf( wrf_input, trim(varsf(n)), tsk, &
                                   dims(1), dims(2), 1, debug)

         if ( keep_tsk_wrf ) then
            call da_get_var_2d_real_cdf( da_file, trim(varsf(n)), tsk_wrfvar, &
                                      dims(1), dims(2), 1, debug)
            !hcl call da_get_var_2d_int_cdf( da_file, 'IVGTYP', ivgtyp, &
            call da_get_var_2d_int_cdf( wrf_input, 'IVGTYP', ivgtyp, &
                                      dims(1), dims(2), 1, debug)
            ! update TSK.
            do j=1,dims(2)
               do i=1,dims(1)
                  if (ivgtyp(i,j) /= iswater)  tsk(i,j)=tsk_wrfvar(i,j)
               end do
            end do
         end if

            call da_put_var_2d_real_cdf( da_file, trim(varsf(n)), tsk, &
                                      dims(1), dims(2), 1, debug)
            deallocate(tsk)
            deallocate(ivgtyp)
            deallocate(tsk_wrfvar)

         !hcl case ('TMN', 'SST', 'VEGFRA', 'ALBBCK', 'SEAICE') ;
         case ('TMN', 'SST', 'VEGFRA', 'ALBBCK', 'SEAICE', 'LANDMASK', 'XLAND') ;
            if ( .not. update_low_bdy ) cycle

            allocate(full2d(dims(1), dims(2)))

            call da_get_var_2d_real_cdf( wrf_input, trim(varsf(n)), full2d, &
                                      dims(1), dims(2), 1, debug)

            call da_put_var_2d_real_cdf( da_file, trim(varsf(n)), full2d, &
                                      dims(1), dims(2), 1, debug)
            deallocate(full2d)

         case ('IVGTYP', 'ISLTYP') ;  !hcl add
            if ( .not. update_low_bdy ) cycle

            allocate(full2dint(dims(1), dims(2)))

            call da_get_var_2d_int_cdf( wrf_input, trim(varsf(n)), full2dint, &
                                      dims(1), dims(2), 1, debug)

            call da_put_var_2d_int_cdf( da_file, trim(varsf(n)), full2dint, &
                                      dims(1), dims(2), 1, debug)
            deallocate(full2dint)

         case ('SNOW', 'RHOSN', 'SNOWH', 'SNOWC') ;
            if ( (.not. update_lsm) .and. (.not. update_low_bdy) ) cycle
            if ( keep_snow_wrf ) cycle
               allocate(full2d(dims(1), dims(2)))

               call da_get_var_2d_real_cdf( wrf_input, trim(varsf(n)), full2d, &
                                      dims(1), dims(2), 1, debug )

               call da_put_var_2d_real_cdf( da_file, trim(varsf(n)), full2d, &
                                      dims(1), dims(2), 1, debug )
               deallocate(full2d)

         case ('CANWAT') ;
            if ( .not. update_lsm ) cycle
               allocate(full2d(dims(1), dims(2)))

               call da_get_var_2d_real_cdf( wrf_input, trim(varsf(n)), full2d, &
                                      dims(1), dims(2), 1, debug )
!               print *,"sum(full2d^2)=", sum(full2d*full2d)

               call da_put_var_2d_real_cdf( da_file, trim(varsf(n)), full2d, &
                                      dims(1), dims(2), 1, debug )
               deallocate(full2d)

         case ('TSLB', 'SMOIS', 'SH2O') ;
            if( .not. update_lsm ) cycle
               allocate(full3d(dims(1), dims(2), dims(3)))

               call da_get_var_3d_real_cdf( wrf_input, trim(varsf(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )
!               print *,"sum(full3d^2)=", sum(full3d*full3d)

               call da_put_var_3d_real_cdf( da_file, trim(varsf(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )
               deallocate(full3d)

         case default ;
            write(unit=stdout,fmt=*) 'It is impossible here. varsf(n)=', trim(varsf(n))
      end select
   end do

   ! check for snow over water
   iostatus(1) = nf_inq_varid(cdfid, 'IVGTYP', varid)
   iostatus(2) = nf_inq_varid(cdfid, 'SNOW',   varid)
   iostatus(3) = nf_inq_varid(cdfid, 'SNOWC',  varid)
   iostatus(4) = nf_inq_varid(cdfid, 'SNOWH',  varid)
   if ( iostatus(1) == 0 ) then
      allocate(snow(dims(1), dims(2)))
      allocate(snowc(dims(1), dims(2)))
      allocate(snowh(dims(1), dims(2)))
      allocate(ivgtyp(dims(1), dims(2)))
      if ( iostatus(1) == 0 ) then
         call da_get_var_2d_int_cdf( da_file, 'IVGTYP', ivgtyp,    &
                               dims(1), dims(2), 1, debug)
      end if
      if ( iostatus(2) == 0 ) then
         call da_get_var_2d_real_cdf( da_file, 'SNOW',    snow,     &
                                dims(1), dims(2), 1, debug)
      end if
      if ( iostatus(3) == 0 ) then
         call da_get_var_2d_real_cdf( da_file, 'SNOWC',  snowc,     &
                                dims(1), dims(2), 1, debug)
      end if
      if ( iostatus(4) == 0 ) then
         call da_get_var_2d_real_cdf( da_file, 'SNOWH',  snowh,     &
                                dims(1), dims(2), 1, debug)
      end if
      if ( iostatus(2) == 0 ) then
         do j = 1, dims(2)
            do i = 1, dims(1)
               if (ivgtyp(i,j) == iswater)  then
                  if ( snow(i,j) > 0.0 ) then
                     write(unit=stdout,fmt=*) 'Remove snow over water at i, j = ', i, j
                     if ( iostatus(2) == 0 ) snow(i,j)  = 0.0
                     if ( iostatus(3) == 0 ) snowc(i,j) = 0.0
                     if ( iostatus(4) == 0 ) snowh(i,j) = 0.0
                  end if
               end if
            end do
         end do
      end if
      if ( iostatus(2) == 0 ) then
         call da_put_var_2d_real_cdf( da_file, 'SNOW',   snow, &
                                dims(1), dims(2), 1, debug)
      end if
      if ( iostatus(3) == 0 ) then
         call da_put_var_2d_real_cdf( da_file, 'SNOWC',  snowc, &
                                dims(1), dims(2), 1, debug)
      end if
      if ( iostatus(4) == 0 ) then
         call da_put_var_2d_real_cdf( da_file, 'SNOWH',  snowh, &
                                dims(1), dims(2), 1, debug)
      end if
      deallocate(snow)
      deallocate(snowc)
      deallocate(snowh)
      deallocate(ivgtyp)
   end if
   
 if ( update_lateral_bdy ) then

   if ( .not. var4d .and. ( east_end < 1 .or. north_end < 1)) then
      write(unit=stdout, fmt='(a)') 'Wrong data for Boundary.'
      stop
   end if

   if(debug) write(unit=stdout,fmt='(/a/)') 'Processing the lateral boundary condition:'

   ! boundary variables
   bdyname(1)='_BXS'
   bdyname(2)='_BXE'
   bdyname(3)='_BYS'
   bdyname(4)='_BYE'

   ! boundary tendancy variables
   tenname(1)='_BTXS'
   tenname(2)='_BTXE'
   tenname(3)='_BTYS'
   tenname(4)='_BTYE'
   !debug = .true.
   do m=1,4
      var_name='MU' // trim(bdyname(m))
      vbt_name='MU' // trim(tenname(m))

      call da_get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug)

      allocate(frst2d(dims(1), dims(2)))
      allocate(scnd2d(dims(1), dims(2)))
      allocate(tend2d(dims(1), dims(2)))
      if ( .not. var4d ) then
      ! Get variable at second time level
      if (time_level > 1) then
         call da_get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
                                   dims(1), dims(2), 2, debug)
      else
         call da_get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                   dims(1), dims(2), 1, debug)
         call da_get_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), 1, debug)
      end if

      if (debug) then
         write(unit=ori_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
              'No.', m, 'Variable: ', trim(vbt_name), &
              'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

         call da_get_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), 1, debug)

         write(unit=ori_unit, fmt='(a, 10i12)') &
              ' old ', (i, i=1,dims(2))
         do j=1,dims(1)
            write(unit=ori_unit, fmt='(i4, 1x, 10e20.7)') &
                  j, (tend2d(j,i), i=1,dims(2))
         end do
      end if

      ! calculate variable at first time level
      if ( time_level < 2 ) &
         scnd2d = frst2d + tend2d * bdyfrq
      select case(m)
      case (1) ;             ! West boundary
         do l=1,dims(2)
            do j=1,dims(1)
               frst2d(j,l)=mu(l,j)
            end do
         end do
      case (2) ;             ! East boundary
         do l=1,dims(2)
            do j=1,dims(1)
               frst2d(j,l)=mu(east_end-l,j)
            end do
         end do
      case (3) ;             ! South boundary
         do l=1,dims(2)
            do i=1,dims(1)
               frst2d(i,l)=mu(i,l)
            end do
         end do
      case (4) ;             ! North boundary
         do l=1,dims(2)
            do i=1,dims(1)
               frst2d(i,l)=mu(i,north_end-l)
            end do
         end do
      case default ;
         write(unit=stdout,fmt=*) 'It is impossible here. mu, m=', m
      end select
      else
          call da_get_var_2d_real_cdf( da_bdy_file, trim(var_name), frst2d, &
                                   dims(1), dims(2), 1, debug)
          call da_get_var_2d_real_cdf( da_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), 1, debug)
          call da_put_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                   dims(1), dims(2), 1, debug)
          call da_put_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), 1, debug)

        if (time_level > 2 ) then
          call da_get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
                                   dims(1), dims(2), 3, debug)
        else
          if ( time_level > 1 ) then
            call da_get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
                                 dims(1), dims(2), 2, debug)
          else
            cycle
          end if
        end if

        frst2d = frst2d + tend2d * bdyfrq
        if ( time_level < 3 ) then
          call da_put_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), 2, debug)
          scnd2d = scnd2d + tend2d * bdyfrq
        end if
      end if

      ! calculate new tendancy 
      do l=1,dims(2)
         do i=1,dims(1)
            tend2d(i,l)=(scnd2d(i,l)-frst2d(i,l))/bdyfrq
         end do
      end do

      if (debug) then
         write(unit=new_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
              'No.', m, 'Variable: ', trim(vbt_name), &
              'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

         write(unit=new_unit, fmt='(a, 10i12)') &
              ' new ', (i, i=1,dims(2))

         do j=1,dims(1)
            write(unit=new_unit, fmt='(i4, 1x, 10e20.7)') &
                  j, (tend2d(j,i), i=1,dims(2))
         end do
      end if

      ! output new variable at first time level
      call da_put_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                dims(1), dims(2), time_id, debug)
      ! output new tendancy 
      call da_put_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                dims(1), dims(2), time_id, debug)

      deallocate(frst2d)
      deallocate(scnd2d)
      deallocate(tend2d)
   end do

   !---------------------------------------------------------------------
   ! For 3D variables
   if ( .not. var4d ) then
   ! Get U
   call da_get_dims_cdf( da_file, 'U', dims, ndims, debug)

   ! call da_get_att_cdf( da_file, 'U', debug)

   allocate(u(dims(1), dims(2), dims(3)))

   ids=1
   ide=dims(1)-1
   jds=1
   jde=dims(2)
   kds=1
   kde=dims(3)

   call da_get_var_3d_real_cdf( da_file, 'U', u, &
                             dims(1), dims(2), dims(3), 1, debug)

   ! do j=1,dims(2)
   !    write(unit=stdout, fmt='(2(a,i5), a, f12.8)') &
   !       'u(', dims(1), ',', j, ',1)=', u(dims(1),j,1)
   ! end do

   ! Get V
   call da_get_dims_cdf( da_file, 'V', dims, ndims, debug)

   ! call da_get_att_cdf( da_file, 'V', debug)

   allocate(v(dims(1), dims(2), dims(3)))

   call da_get_var_3d_real_cdf( da_file, 'V', v, &
                             dims(1), dims(2), dims(3), 1, debug)

   ! do i=1,dims(1)
   !    write(unit=stdout, fmt='(2(a,i5), a, f12.8)') &
   !       'v(', i, ',', dims(2), ',1)=', v(i,dims(2),1)
   ! end do

   if (debug) then
      write(unit=stdout, fmt='(a,e20.12,4x)') &
           'Before couple Sample u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'Before couple Sample v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   end if

   !---------------------------------------------------------------------
   ! Couple u, v.
   call da_couple_uv ( u, v, mu, mub, msfu, msfv, ids, ide, jds, jde, kds, kde)

   if (debug) then
      write(unit=stdout, fmt='(a,e20.12,4x)') &
           'After  couple Sample u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'After  couple Sample v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   end if
   end if
   !---------------------------------------------------------------------
   !For 3D variables

   do n=1,num3d
      if (debug) write(unit=stdout, fmt='(a, i3, 2a)') 'Processing: var3d(', n, ')=', trim(var3d(n))
      var_pref=trim(var3d(n))
      if ( .not. var4d ) then
      call da_get_dims_cdf( da_file, trim(var3d(n)), dims, ndims, debug)

      allocate(full3d(dims(1), dims(2), dims(3)))

      east_end=dims(1)+1
      north_end=dims(2)+1

      select case(trim(var3d(n)))
      case ('U') ;           ! U
         ! var_pref='R' // trim(var3d(n))
         !var_pref=trim(var3d(n))
         full3d(:,:,:)=u(:,:,:)
      case ('V') ;           ! V 
         ! var_pref='R' // trim(var3d(n))
         !var_pref=trim(var3d(n))
         full3d(:,:,:)=v(:,:,:)
      case ('W') ;
         ! var_pref = 'R' // trim(var3d(n))
         !var_pref = trim(var3d(n))

         call da_get_var_3d_real_cdf( da_file, trim(var3d(n)), &
            full3d, dims(1), dims(2), dims(3), 1, debug)

         if (debug) then
            write(unit=stdout, fmt='(3a,e20.12,4x)') &
                 'Before couple Sample ', trim(var3d(n)), &
                 '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
         end if

         do k=1,dims(3)
            do j=1,dims(2)
               do i=1,dims(1)
                  full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))/msfm(i,j)
               end do
            end do
         end do

         if (debug) then
            write(unit=stdout, fmt='(3a,e20.12,4x)') &
                 'After  couple Sample ', trim(var3d(n)), &
                 '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
         end if
      case ('T', 'PH') ;
         !var_pref=trim(var3d(n))
 
         call da_get_var_3d_real_cdf( da_file, trim(var3d(n)), &
            full3d, dims(1), dims(2), dims(3), 1, debug)

         if (debug) then
            write(unit=stdout, fmt='(3a,e20.12,4x)') &
                 'Before couple Sample ', trim(var3d(n)), &
                 '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
         end if

         do k=1,dims(3)
            do j=1,dims(2)
               do i=1,dims(1)
                  full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
               end do
            end do
         end do

            if (debug) then
               write(unit=stdout, fmt='(3a,e20.12,4x)') &
                    'After  couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            end if
      case ('QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP') ;
         ! var_pref='R' // var3d(n)(1:2)
         ! var_pref=var3d(n)(1:2)
         !var_pref=var3d(n)
 
         call da_get_var_3d_real_cdf( da_file, trim(var3d(n)), &
            full3d, dims(1), dims(2), dims(3), 1, debug)

         if (debug) then
            write(unit=stdout, fmt='(3a,e20.12,4x)') &
                 'Before couple Sample ', trim(var3d(n)), &
                 '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
         end if

         do k=1,dims(3)
            do j=1,dims(2)
               do i=1,dims(1)
                  full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
               end do
            end do
         end do

         if (debug) then
            write(unit=stdout, fmt='(3a,e20.12,4x)') &
                 'After  couple Sample ', trim(var3d(n)), &
                 '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
         end if
      case default ;
         write(unit=stdout,fmt=*) 'It is impossible here. var3d(', n, ')=', trim(var3d(n))
      end select
      end if
      do m=1,4
         var_name=trim(var_pref) // trim(bdyname(m))
         vbt_name=trim(var_pref) // trim(tenname(m))

         if (debug) write(unit=stdout, fmt='(a, i3, 2a)') &
            'Processing: bdyname(', m, ')=', trim(var_name)

         call da_get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug)

         allocate(frst3d(dims(1), dims(2), dims(3)))
         allocate(scnd3d(dims(1), dims(2), dims(3)))
         allocate(tend3d(dims(1), dims(2), dims(3)))
         if ( .not. var4d ) then
         ! Get variable at second time level
         if (time_level > 1) then
            call da_get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
                                      dims(1), dims(2), dims(3), 2, debug)
         else
            call da_get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                      dims(1), dims(2), dims(3), 1, debug)
            call da_get_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                      dims(1), dims(2), dims(3), 1, debug)
         end if

         if (debug) then
            write(unit=ori_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
                 'No.', m, 'Variable: ', trim(vbt_name), &
                 'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

            call da_get_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                      dims(1), dims(2), dims(3), 1, debug)

            write(unit=ori_unit, fmt='(a, 10i12)') &
                 ' old ', (i, i=1,dims(3))
            do j=1,dims(1)
               write(unit=ori_unit, fmt='(i4, 1x, 10e20.7)') &
                     j, (tend3d(j,dims(2)/2,i), i=1,dims(3))
            end do
         end if

         if (time_level < 2) scnd3d = frst3d+tend3d*bdyfrq 
         select case(trim(bdyname(m)))
         case ('_BXS') ;             ! West boundary
            do l=1,dims(3)
            do k=1,dims(2)
            do j=1,dims(1)
               frst3d(j,k,l)=full3d(l,j,k)
            end do
            end do
            end do
         case ('_BXE') ;             ! East boundary
            do l=1,dims(3)
            do k=1,dims(2)
            do j=1,dims(1)
               frst3d(j,k,l)=full3d(east_end-l,j,k)
            end do
            end do
            end do
         case ('_BYS') ;             ! South boundary
            do l=1,dims(3)
            do k=1,dims(2)
            do i=1,dims(1)
               frst3d(i,k,l)=full3d(i,l,k)
            end do
            end do
            end do
         case ('_BYE') ;             ! North boundary
            do l=1,dims(3)
            do k=1,dims(2)
            do i=1,dims(1)
               frst3d(i,k,l)=full3d(i,north_end-l,k)
            end do
            end do
            end do
         case default ;
            write(unit=stdout,fmt=*) 'It is impossible here.'
            write(unit=stdout,fmt=*) 'bdyname(', m, ')=', trim(bdyname(m))
            stop
         end select
         else
           call da_get_var_3d_real_cdf( da_bdy_file, trim(var_name), frst3d, &
                                   dims(1), dims(2), dims(3), 1, debug)
           call da_get_var_3d_real_cdf( da_bdy_file, trim(vbt_name), tend3d, &
                                  dims(1), dims(2), dims(3), 1, debug)

           call da_put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                  dims(1), dims(2), dims(3), 1, debug)
           call da_put_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                  dims(1), dims(2), dims(3), 1, debug)

           if (time_level > 2 ) then
             call da_get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
                                  dims(1), dims(2), dims(3), 3, debug)
           else
             if ( time_level > 1 ) then
               call da_get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
                              dims(1), dims(2), dims(3), 2, debug)
             else
               cycle
             end if
           end if

           frst3d = frst3d + tend3d * bdyfrq
           if ( time_level < 3 ) then 
             call da_get_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                      dims(1), dims(2), dims(3), 2, debug)
             scnd3d = scnd3d + tend3d * bdyfrq 
           end if
        end if

        if (debug)  write(unit=stdout, fmt='(a, i3, 2a)') &
            'cal. tend: bdyname(', m, ')=', trim(vbt_name)

         ! calculate new tendancy 
         do l=1,dims(3)
            do k=1,dims(2)
               do i=1,dims(1)
                  tend3d(i,k,l)=(scnd3d(i,k,l)-frst3d(i,k,l))/bdyfrq
               end do
            end do
         end do

         if (debug) then
            write(unit=new_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
                 'No.', m, 'Variable: ', trim(vbt_name), &
                 'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

            write(unit=new_unit, fmt='(a, 10i12)') &
                 ' new ', (i, i=1,dims(3))

            do j=1,dims(1)
               write(unit=new_unit, fmt='(i4, 1x, 10e20.7)') &
                     j, (tend3d(j,dims(2)/2,i), i=1,dims(3))
            end do
         end if

         ! output new variable at first time level
         call da_put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                dims(1), dims(2), dims(3), time_id, debug)
         call da_put_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                   dims(1), dims(2), dims(3), time_id, debug)

         deallocate(frst3d)
         deallocate(scnd3d)
         deallocate(tend3d)
      end do
      
      if ( .not. var4d ) deallocate(full3d)
   end do
   if ( .not. var4d ) then
   deallocate(mu)
   deallocate(mub)
   deallocate(msfu)
   deallocate(msfv)
   deallocate(u)
   deallocate(v)
   end if
 end if

 write(unit=stdout,fmt=*) &
    '=================================================================='
 if ( update_lateral_bdy ) then
    write(unit=stdout,fmt=*) 'Lateral boundary tendency updated.'
 end if
 if ( update_low_bdy ) then
    write(unit=stdout,fmt=*) 'Low boundary updated with wrf_input fields.'
 end if
 if ( update_lsm ) then
    write(unit=stdout,fmt=*) 'LSM variables updated with wrf_input fields.'
 end if

   if (io_status == 0) &
      write (unit=stdout,fmt=*) "*** Update_bc completed successfully ***"

end program da_update_bc

