subroutine da_get_satzen ( lat,lon,sate_index,theta_true )

   !---------------------------------------------------------------------------
   ! Purpose: calculate geostationary satellite_zenith_angle
   !
   ! History: original code provided by Jun Li of CIMSS
   !          adapted by Chun Yang
   !---------------------------------------------------------------------------

   implicit none

   real,    intent(in)  :: lat, lon
   integer, intent(in)  :: sate_index
   real,    intent(out) :: theta_true

   real, parameter :: pi = 3.1415926535
   real, parameter :: r_earth = 6370.
   real, parameter :: r_sun = 35800.
   real            :: alat, alon, alon_sat
   real            :: theta, r_tmp, theta_tmp

   alat = lat
   alon = lon

   if (sate_index .eq. 11) then
      alon_sat = -135.*pi/180.
   else if (sate_index .eq. 12) then
      alon_sat = -60.*pi/180.
   else if (sate_index .eq. 13) then
      alon_sat = -75.*pi/180.
   else if (sate_index .eq. 14) then
      alon_sat = -105.*pi/180.
   else if (sate_index .eq. 15) then
      alon_sat = -135.*pi/180.
   else
      call da_error(__FILE__,__LINE__, (/"satellite id has to be 11-15"/))
   end if

   alat = alat*pi/180.
   alon = alon*pi/180.
   theta = abs(alon-alon_sat)
   r_tmp = (2.0*r_earth*sin(theta/2.)-r_earth*(1.0-cos(alat))*sin(theta/2.))**2 &
           +(2.0*r_earth*sin(alat/2.))**2-(r_earth*(1.0-cos(alat))*sin(theta/2.))**2
   r_tmp = sqrt(r_tmp)
   theta_true = 2*asin(r_tmp/r_earth/2.)
   theta_tmp = atan(r_earth*sin(theta_true)/(r_sun+r_earth*(1.0-sin(theta_true))))
   theta_true = (theta_true+theta_tmp)*180./pi

end subroutine da_get_satzen
