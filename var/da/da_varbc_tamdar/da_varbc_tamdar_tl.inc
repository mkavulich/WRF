   subroutine da_varbc_tamdar_tl (cv_size, cv, iv, y)

   !-------------------------------------------------------!
   !  PURPOSE: Apply bias correction to TAMDAR innovations !
   !  METHOD: delta_d    = y - delta_bc                    !
   !          y          = H (delta_x)                     !
   !          delta_bc   = SUM( delta_beta_i Pred_i )      !
   !                                                       !
   !  Called from da_transform_vtoy.inc                    !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   integer, intent(in)           :: cv_size         ! Size of cv array.
   real, intent(in)              :: cv(1:cv_size)   ! control variables.
   type (iv_type), intent(in)    :: iv              ! O-B structure.
   type (y_type), intent(inout)  :: y               ! y = h (xa)

   integer                       :: isn,ivar,iflt,ipred,iobs,iphase
   real, allocatable             :: varbc_param_tl(:)
   real                          :: delta_bc

   if (trace_use) call da_trace_entry("da_varbc_tamdar_tl")

   !----------------------------------------------------------!

   allocate( varbc_param_tl(iv%tamdar_varbc(1)%npred) )

   do ivar = 1, iv%tamdar_varbc(1)%nbcvar
      do iflt = 1, iv%tamdar_varbc(ivar)%ntflight
         if( iv%tamdar_varbc(ivar)%nobs_sum(iflt,iv%tamdar_varbc(ivar)%nphase+1) >= varbc_tamdar_nobsmin ) then
            do iphase = 1,iv%tamdar_varbc(ivar)%nphase

               varbc_param_tl(:) = 0.0

               if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) >= varbc_tamdar_nobsmin) then

                  if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                     do ipred = 1, iv%tamdar_varbc(ivar)%npred
                        varbc_param_tl(ipred) = SUM( cv(iv%tamdar_varbc(ivar)%index(ipred,iflt,iphase)) * &
                                                     iv%tamdar_varbc(ivar)%vtox(ipred,1:iv%tamdar_varbc(ivar)%npred,iflt,iphase) )
                     end do

                     do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                        isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)
                        delta_bc = SUM( varbc_param_tl(1:iv%tamdar_varbc(ivar)%npred)* &
                                        iv%tamdar_varbc(ivar)%pred(1:iv%tamdar_varbc(ivar)%npred,iflt,iphase) )
                        if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(isn)%t(1)%qc >= 0 ) y%tamdar(isn)%t(1) = y%tamdar(isn)%t(1) + delta_bc
                        if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(isn)%q(1)%qc >= 0 ) y%tamdar(isn)%q(1) = y%tamdar(isn)%q(1) + delta_bc
                     end do
                  end if

               end if
            end do
         end if
      end do
   end do

   deallocate(varbc_param_tl)

   !----------------------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_tl")

   end subroutine da_varbc_tamdar_tl
