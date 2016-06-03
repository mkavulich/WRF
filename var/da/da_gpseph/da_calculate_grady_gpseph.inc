subroutine da_calculate_grady_gpseph(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re (residual) vector
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector
   type (y_type),  intent(inout)  :: re          ! Residual vector
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_gpseph")

   do n=1, iv%info(gpseph)%nlocal
      if ( iv%info(gpseph)%levels(n) < 1 ) cycle
      do k=1, iv%info(gpseph)%levels(n)
         !hcl-check
         ! Since we calculate the nonlocal thing globally, we have to exclude the obs located at boundaries.
         if (iv%gpseph(n)%eph(k)%qc < obs_qc_pointer .or. .not. iv%info(gpseph)%proc_domain(1,n) ) then
            re%gpseph(n)%eph(k) = 0.0
            jo_grad_y%gpseph(n)%eph(k) = 0.0
         else
            jo_grad_y%gpseph(n)%eph(k) = -re%gpseph(n)%eph(k) / (iv%gpseph(n)%eph(k)%error * iv%gpseph(n)%eph(k)%error)
         end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_gpseph")

end subroutine da_calculate_grady_gpseph
