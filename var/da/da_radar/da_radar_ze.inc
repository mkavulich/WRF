subroutine da_radar_ze(ra,rs,rg,tc1,ro1,rze)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: ra, rs, rg, tc1, ro1
   real, intent(inout) :: rze

   real                :: zrr, zds, zws, zgr, zerr, zews, zeds, zegr, czr

!  INTRINSIC LOG10

   zrr = 3.63*1.00e+9  ! rainwater
   zws = 4.26*1.00e+11 ! wet snow
   zds = 9.80*1.00e+8  ! dry snow
   zgr = 4.33*1.00e+10 ! graupel

   zerr = zrr*(ro1*ra)**1.75
   zews = zws*(ro1*rs)**1.75
   zeds = zds*(ro1*rs)**1.75
   zegr = zgr*(ro1*rg)**1.75

   if (tc1 .ge. 5.0) then
      rze = zerr
   else if (tc1 .lt. 5.0 .and. tc1 .ge. -5.0) then
      czr = (tc1+5.0)/10.0
      if (tc1 .le. 0.0) then
         rze = czr*zerr + (1.0-czr)*(zeds+zegr)
      else
         rze = czr*zerr + (1.0-czr)*(zews+zegr)
      end if
   else if (tc1.lt.-5.0) then
      rze = zeds + zegr
   end if

end subroutine da_radar_ze

