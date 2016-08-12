   subroutine da_varbc_tamdar_direct (iv,ob)

   !-------------------------------------------------------!
   !  PURPOSE: Apply bias correction to TAMDAR innovations !
   !                                                       !
   !  Called from da_get_innov_vector.inc                  !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   type (iv_type), intent(inout)  :: iv       ! O-B structure.
   type (y_type),  intent(inout)  :: ob       ! Observation structure.

   real                           :: bc,bias
   integer                        :: isn,iflt,iobs,ipred,ivar,iphase,nobs
   character(len=1)               :: cvar0
   character(len=3)               :: cvar1

   if (trace_use) call da_trace_entry("da_varbc_tamdar_direct")

   !----------------------------------------------------------!

   if (rootproc) &
      write(unit=varbc_tamdar_unit,fmt='(//8X,A/)')'Calculating corrected innovation'

   do ivar = 1, iv%tamdar_varbc(1)%nbcvar
      do iflt = 1, iv%tamdar_varbc(ivar)%ntflight
         do iphase = 1, iv%tamdar_varbc(ivar)%nphase

            nobs = 0
            bias = 0.

            bc = varbc_tamdar_tune_bc * &
                 SUM( iv%tamdar_varbc(ivar)%param(1:iv%tamdar_varbc(ivar)%npred,iflt,iphase) * &
                      iv%tamdar_varbc(ivar)%pred(1:iv%tamdar_varbc(ivar)%npred,iflt,iphase) )

            if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 ) cvar0 = 'T'
            if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 ) cvar0 = 'Q'
            if( iphase == 1 ) cvar1 = 'asc'
            if( iphase == 2 ) cvar1 = 'des'

            if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) >= varbc_tamdar_nobsmin) then
               if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0 ) then
                  do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                     isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)
                     if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(isn)%t(1)%qc >= 0 ) then
                         bias = bias + iv%tamdar(isn)%t(1)%inv
                         nobs = nobs + 1
                     end if
                     if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(isn)%q(1)%qc >= 0 ) then
                         bias = bias + iv%tamdar(isn)%q(1)%inv
                         nobs = nobs + 1
                     end if
                  end do
               end if

               if( wrf_dm_sum_integer(nobs) > 0 ) bias = wrf_dm_sum_real(bias)/wrf_dm_sum_integer(nobs)

               if( bias * iv%tamdar_varbc(ivar)%param(1,iflt,iphase) < 0.0 ) then
                   iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) = 0
                   if (rootproc) &
                      write(unit=varbc_tamdar_unit,fmt='(13X,A,A,A,I5,A,1X,A,A,2F12.7)') &
                           'Variable ',cvar0, ' in',iv%tamdar_varbc(ivar)%tail_id(iflt),' changed to ifuse=0 at', cvar1, &
                           ', bias/param:', bias, iv%tamdar_varbc(ivar)%param(1,iflt,iphase)
               end if

               if ( iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0 ) then
                    do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                       isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)
                       if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(isn)%t(1)%qc >= 0 ) then
                           iv%tamdar(isn)%t(1)%inv = iv%tamdar(isn)%t(1)%inv - bc
                       end if
                       if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(isn)%q(1)%qc >= 0 ) then
                           iv%tamdar(isn)%q(1)%inv = iv%tamdar(isn)%q(1)%inv - bc
                       end if
                    end do
               end if

               if (rootproc .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                   write(unit=varbc_tamdar_unit,fmt='(13X,A,1X,A,A,I5,2X,A,A,2F12.7)') &
                        'Variable',cvar0,' bias corrected for ',iv%tamdar_varbc(ivar)%tail_id(iflt), cvar1, &
                        '  : (BC,OMB bias)', bc, bias
               end if
            end if

         end do
      end do
   end do

   !----------------------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_direct")

   end subroutine da_varbc_tamdar_direct
