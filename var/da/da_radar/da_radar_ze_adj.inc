subroutine da_radar_ze_adj(ra, rab, rs, rsb, rg, rgb, tc1, ro1, rzeb)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: ra, rs, rg, tc1, ro1, rzeb 
   real, intent(inout) :: rab, rsb, rgb

   real                :: zrr, zds, zws, zgr, zerrb, zewsb, zedsb, zegrb, czr

   zrr = 3.63*1.00e+9  ! rainwater
   zws = 4.26*1.00e+11 ! wet snow
   zds = 9.80*1.00e+8  ! dry snow
   zgr = 4.33*1.00e+10 ! graupel

   zerrb = 0.0
   zewsb = 0.0
   zedsb = 0.0
   zegrb = 0.0

   if (tc1 .ge. 5.0) then
      zerrb = zerrb + rzeb
   else if (tc1 .le. 5.0 .and. tc1 .ge. -5.0) then
      czr = (tc1+5.0)/10.0
      if (tc1.le.0.0) then
         zerrb = zerrb + czr * rzeb
         zedsb = zedsb + (1.0-czr) * rzeb
         zegrb = zegrb + (1.0-czr) * rzeb
      else
         zerrb = zerrb + czr * rzeb
         zewsb = zewsb + (1.0-czr) * rzeb
         zegrb = zegrb + (1.0-czr) * rzeb
      end if
   else if (tc1 .lt. -5.0) then
!!ATTENTION    rze = zeds + zegr
    zedsb = zedsb + rzeb
    zegrb = zegrb + rzeb
   end if

  rgb = 1.75*(ro1*rg)**0.75*zgr*ro1*zegrb
  rsb = 1.75*(ro1*rs)**0.75*zws*ro1*zewsb + 1.75*(ro1*rs)**0.75*zds*ro1*zedsb
  rab = 1.75*(ro1*ra)**0.75*zrr*ro1*zerrb

end subroutine da_radar_ze_adj
