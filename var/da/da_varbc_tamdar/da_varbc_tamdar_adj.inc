   subroutine da_varbc_tamdar_adj (cv_size, cv, iv, y)

   !-------------------------------------------------------!
   !  PURPOSE: Apply bias correction to TAMDAR innovations !
   !                                                       !
   !  Called from da_transform_vtoy.inc                    !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   integer, intent(in)           :: cv_size         ! Size of cv array.
   real, intent(inout)           :: cv(1:cv_size)   ! control variables.
   type (iv_type), intent(inout) :: iv              ! O-B structure.
   type (y_type), intent(in)     :: y               ! y = h (xa)

   integer                       :: isn,iflt,iobs,ipred,iphase,ivar,npred
   real, allocatable             :: varbc_param_adj(:)
   real                          :: yh, cv_local

   if (trace_use) call da_trace_entry("da_varbc_tamdar_adj")

   !--------------------------------------------------------!

   npred  = iv%tamdar_varbc(1)%npred

   allocate( varbc_param_adj(npred) )

   do ivar = 1, iv%tamdar_varbc(1)%nbcvar
      do iflt = 1, iv%tamdar_varbc(ivar)%ntflight
         if( iv%tamdar_varbc(ivar)%nobs_sum(iflt,iv%tamdar_varbc(1)%nphase+1) >= varbc_tamdar_nobsmin ) then
            do iphase = 1, iv%tamdar_varbc(1)%nphase

               varbc_param_adj(:) = 0.0

               if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) >= varbc_tamdar_nobsmin) then

                  if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                     do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                        isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)
                        yh = 0
                        if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(isn)%t(1)%qc >= 0 ) yh = y%tamdar(isn)%t(1)
                        if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(isn)%q(1)%qc >= 0 ) yh = y%tamdar(isn)%q(1)
                        varbc_param_adj(1:npred) = varbc_param_adj(1:npred) + &
                                                   yh * iv%tamdar_varbc(ivar)%pred(1:npred,iflt,iphase)
                     end do
                  end if

                  do ipred = 1, npred
                     cv_local = SUM( varbc_param_adj(1:npred) * iv%tamdar_varbc(ivar)%vtox(ipred,1:npred,iflt,iphase) )
                     cv_local = wrf_dm_sum_real(cv_local)

                     if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) &
                        cv(iv%tamdar_varbc(ivar)%index(ipred,iflt,iphase)) = cv( iv%tamdar_varbc(ivar)%index(ipred,iflt,iphase) ) + cv_local
                  end do

               end if
            end do
         end if
      end do
   end do

   deallocate(varbc_param_adj)

   !--------------------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_adj")

   end subroutine da_varbc_tamdar_adj
