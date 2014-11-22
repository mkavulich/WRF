subroutine da_radar_ze(ra, rad, rs, rsd, rg, rgd, tc1, ro1, rzed)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: ra, rad, rs, rsd, rg, rgd, tc1, ro1 
   real, intent(inout) :: rzed

   real                :: zrr, zds, zws, zgr, zerrd, zewsd, zedsd, zegrd, czr

!  INTRINSIC LOG10

   zrr = 3.63*1.00e+9   ! rainwater
   zws = 4.26*1.00e+11  ! wet snow
   zds = 9.80*1.00e+8   ! dry snow
   zgr = 4.33*1.00e+10  ! graupel

   zerrd = zrr*1.75*(ro1*ra)**0.75*ro1*rad
   zewsd = zws*1.75*(ro1*rs)**0.75*ro1*rsd
   zedsd = zds*1.75*(ro1*rs)**0.75*ro1*rsd
   zegrd = zgr*1.75*(ro1*rg)**0.75*ro1*rgd

   if (tc1 .ge. 5.0) then
      rzed = zerrd
   else if (tc1 .le. 5.0 .and. tc1 .ge. -5.0) then
      czr = (tc1+5.0)/10.0
      if (tc1.le.0.0) then
         rzed = czr*zerrd + (1.0-czr)*(zedsd+zegrd)
      else
         rzed = czr*zerrd + (1.0-czr)*(zewsd+zegrd)
      end if
   else if (tc1 .LT. -5.0) then
      rzed = zedsd + zegrd
   end if

end subroutine da_radar_ze

