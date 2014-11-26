subroutine da_check_xtoy_adjoint_div(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n, k          ! Loop counter.
   real :: typical_div_rms=0.0001
   if (trace_use) call da_trace_entry("da_check_xtoy_adjoint_div")

   do n=iv%info(idxdiv)%n1,iv%info(idxdiv)%n2
      
            adjtest_lhs = adjtest_lhs + &
               (y%div(n)/typical_div_rms)**2   
        

         pertile_lhs = pertile_lhs + &
            (y%div(n)/typical_div_rms)**2   

         y%div(n)= y%div(n) / typical_div_rms ** 2
         
   end do

end subroutine da_check_xtoy_adjoint_div


