subroutine da_check_max_iv_gpseph(iv,it,num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: max_iv check for gpseph obs
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it       ! External iteration.
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer                      :: k,n,ipr
   logical                      :: failed

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_gpseph")

   !----------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !----------------------------------------------------------------------------

   do n = iv%info(gpseph)%n1,iv%info(gpseph)%n2
      if ( iv%info(gpseph)%levels(n) < 1 ) cycle
      do k = 1, iv%info(gpseph)%levels(n)
         !hcl call da_get_print_lvl(iv%gpseph(n)%p(k)%inv,ipr)
         failed=.false.
         if ( iv%gpseph(n)%eph(k)%qc >= obs_qc_pointer )  then
            call da_max_error_qc(it, iv%info(gpseph), n, iv%gpseph(n)%eph(k), max_error_eph, failed)
            if ( iv%info(gpseph)%proc_domain(k,n) ) then
               !hcl num_qcstat_conv(1,gpseph,8,ipr) = num_qcstat_conv(1,gpseph,8,ipr) + 1
               if (failed) then
                  !hcl num_qcstat_conv(2,gpseph,8,ipr) = num_qcstat_conv(2,gpseph,8,ipr) + 1
                  !hcl write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                  !hcl    'gpseph',ob_vars(8),iv%info(gpseph)%lat(k,n),iv%info(gpseph)%lon(k,n), &
                  !hcl    0.01*iv%gpseph(n)%p(k)%inv
               end if
            end if
         end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_gpseph")

end subroutine da_check_max_iv_gpseph
