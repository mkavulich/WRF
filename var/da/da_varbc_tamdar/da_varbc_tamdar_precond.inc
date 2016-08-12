   Subroutine da_varbc_tamdar_precond (iv)

   !-------------------------------------------------------!
   !  PURPOSE: Estimate Hessian for Preconditioning        !
   !                                                       !
   !  Called from da_get_innov_vector.inc                  !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   type (iv_type), intent (inout) :: iv

   integer                        :: i,j,iflt,iobs,obs_sn,ivar,ipred,iphase,npred
   real                           :: hessian_local,bgerr_local,pred_i,pred_j,verr

   real, allocatable              :: hessian(:,:)
   real*8, allocatable            :: eignvec(:,:), eignval(:)

   if (trace_use) call da_trace_entry("da_varbc_tamdar_precond")

   !-------------------------------------------------------------!

   if (rootproc) &
      write(unit=varbc_tamdar_unit,fmt='(//8X,A)') 'Estimating hessian for preconditioning'

   npred = iv%tamdar_varbc(1)%npred

   allocate ( hessian(npred, npred) )
   allocate ( eignvec(npred, npred) )
   allocate ( eignval(npred)        )

   do ivar = 1, iv%tamdar_varbc(1)%nbcvar
      do iflt = 1, iv%tamdar_varbc(ivar)%ntflight
         if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iv%tamdar_varbc(1)%nphase+1) >= varbc_tamdar_nobsmin ) then
            do iphase = 1, iv%tamdar_varbc(1)%nphase

               if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) >= varbc_tamdar_nobsmin) then

                   hessian(:,:) = 0.
                   eignvec(:,:) = 0.
                   eignval(:)   = 0.

                   bgerr_local = 0.0
                   if(iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                      do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                         obs_sn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)

                         verr = 0.0
                         if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(obs_sn)%t(1)%qc >= 0 ) verr = iv%tamdar(obs_sn)%t(1)%error
                         if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(obs_sn)%q(1)%qc >= 0 ) verr = iv%tamdar(obs_sn)%q(1)%error
                         bgerr_local = bgerr_local + verr**2.0/varbc_tamdar_nbgerr
                      end do
                   end if
                   bgerr_local = wrf_dm_sum_real(bgerr_local)

                   do i = 1, npred
                      do j = i, npred
                         hessian_local = 0.0

                         if(iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                            do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                               obs_sn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)
                               pred_i = iv%tamdar_varbc(ivar)%pred(i,iflt,iphase)
                               pred_j = iv%tamdar_varbc(ivar)%pred(j,iflt,iphase)

                               verr = 0.0
                               if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(obs_sn)%t(1)%qc >= 0 ) verr = iv%tamdar(obs_sn)%t(1)%error
                               if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(obs_sn)%q(1)%qc >= 0 ) verr = iv%tamdar(obs_sn)%q(1)%error
                               if( verr > 0. ) &
                                  hessian_local = hessian_local + pred_i*pred_j/verr**2.0
                            end do
                         end if
                         hessian(i,j) = wrf_dm_sum_real(hessian_local)
                         hessian(j,i) = hessian(i,j)
                      end do

                      if( iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) > 0 ) &
                          iv%tamdar_varbc(ivar)%bgerr(i,iflt,iphase) = bgerr_local/iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase)
                      if (iv%tamdar_varbc(ivar)%bgerr(i,iflt,iphase) > 0.) &
                          hessian(i,i)=hessian(i,i)+1.0/iv%tamdar_varbc(ivar)%bgerr(i,iflt,iphase)
                   end do

                   hessian = hessian/varbc_tamdar_factor**2.0

                   if(iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then

                      call da_eof_decomposition(npred, hessian(1:npred,1:npred), &
                                                eignvec(1:npred,1:npred),eignval(1:npred))

                      if (ANY( eignval(1:npred) <= 0 )) then
                          write(unit=stdout,fmt='(A,I8,A,2F12.5)') &
                               'VARBC_TAMDAR: non-positive Hessian for tail_id ', iv%tamdar_varbc(ivar)%tail_id(iflt), &
                               '. Eigenvalues =',eignval(1:npred)

                          do i = 1, npred
                             if (hessian(i,i) > 0) iv%tamdar_varbc(ivar)%vtox(i,i,iflt,iphase) = 1.0/sqrt(hessian(i,i))
                          end do
                      else
!=======================================================
! Preconditioning = inverse square root of Hessian
!
! Hermite real matrix:
!                 A = Q*D*Q^T
!           A^(1/k) = Q*D^(1/K)*Q^T
!=======================================================
                          do i = 1, npred
                          do j = i, npred
                             iv%tamdar_varbc(ivar)%vtox(i,j,iflt,iphase) = SUM( eignvec(i,1:npred)* &
                                                                           sqrt(1.0/eignval(1:npred)) * &
                                                                           eignvec(j,1:npred) )
                             iv%tamdar_varbc(ivar)%vtox(j,i,iflt,iphase) = iv%tamdar_varbc(ivar)%vtox(i,j,iflt,iphase)
                          end do
                          end do
                      end if
                   end if
               end if
            end do
         end if
      end do   
   end do

   deallocate(hessian, eignvec, eignval)

   if (rootproc) &
      write(unit=varbc_tamdar_unit,fmt='(//8X,A)') 'End estimating hessian for preconditioning'

   !-------------------------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_precond")

   end subroutine da_varbc_tamdar_precond
