   subroutine da_varbc_tamdar_update (it, cv_size, cv, iv)

   !-------------------------------------------------------!
   !  PURPOSE: Update and write out parameters             !
   !                                                       !
   !  Called from da_solve.inc                             !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   integer,        intent(in)    :: it          !outer loop counting
   integer,        intent(in)    :: cv_size
   real,           intent(in)    :: cv(cv_size) ! Control variable structure.
   type (iv_type), intent(inout) :: iv          ! Obs. increment structure.

   integer                       :: i,ivar,ipred,iflt,iphase,npred
   integer                       :: iunit,iost
   character(len=filename_len)   :: filename
   character(len=99)             :: fmt_param
   character(len=1)              :: cvar0
   character(len=3)              :: cvar1
   character(len=6)              :: cvar2
   real                          :: sum_tl

   character(len=5), allocatable :: tail_id(:)
   real, allocatable             :: varbc_param_tl(:,:,:),param(:,:,:,:)
   integer, allocatable          :: ifuse(:,:,:)

   if (trace_use) call da_trace_entry("da_varbc_tamdar_update")

   !----------------------------------------------!

   npred = iv%tamdar_varbc(1)%npred

   if (rootproc) write(unit=varbc_tamdar_unit,fmt='(//8X,A)') 'Updating parameters'

   allocate( varbc_param_tl(npred,iv%tamdar_varbc(1)%ntflight,iv%tamdar_varbc(1)%nphase) )

   do ivar = 1, iv%tamdar_varbc(1)%nbcvar
      varbc_param_tl(:,:,:)=0.0
      do iflt = 1, iv%tamdar_varbc(ivar)%ntflight
         if( iv%tamdar_varbc(ivar)%nobs_sum(iflt,iv%tamdar_varbc(1)%nphase+1) >= varbc_tamdar_nobsmin ) then
            do iphase=1, iv%tamdar_varbc(1)%nphase 

               if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) >= varbc_tamdar_nobsmin) then

                  if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                     do ipred = 1, npred
                        varbc_param_tl(ipred,iflt,iphase) = SUM( cv(iv%tamdar_varbc(ivar)%index(ipred,iflt,iphase)) * &
                                                            iv%tamdar_varbc(ivar)%vtox(ipred,1:npred,iflt,iphase) )
                     end do
                  end if

                  do ipred = 1, npred
                     sum_tl = wrf_dm_sum_real( varbc_param_tl(ipred,iflt,iphase) ) * varbc_tamdar_param_sen(ivar)
                     iv%tamdar_varbc(ivar)%param(ipred,iflt,iphase) = iv%tamdar_varbc(ivar)%param(ipred,iflt,iphase) + sum_tl

                     if (rootproc .and. iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0) then
                         if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 ) cvar0 = 'T'
                         if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 ) cvar0 = 'Q'
                         if( iphase == 1 ) cvar1 = 'asc'
                         if( iphase == 2 ) cvar1 = 'des'
                         if( ipred == 1 ) cvar2 = 'param1'
                         if( ipred == 2 ) cvar2 = 'param2'
                         if( ipred == 3 ) cvar2 = 'param3'
                         if( ipred == 4 ) cvar2 = 'param4'

                         write(unit=varbc_tamdar_unit,fmt='(13X,A,A,A,F13.7,A,I5,3X,A,2X,A)') &
                              'Variable  ',cvar0,'  parameter updated : ',sum_tl, ' in ',iv%tamdar_varbc(1)%tail_id(iflt),cvar2,cvar1
                     end if
                  end do
               end if
            end do
         end if
      end do
   end do

   if (.not. rootproc) then
      deallocate(varbc_param_tl)
      if (trace_use) call da_trace_exit("da_varbc_tamdar_update")
      return
   end if

   !-------------------------------!
   !  Write VARBC_TAMDAR.out file  !
   !-------------------------------!

   allocate( tail_id(iv%tamdar_varbc(1)%ntflight) )
   allocate( ifuse  (iv%tamdar_varbc(1)%ntflight,iv%tamdar_varbc(1)%nphase,iv%tamdar_varbc(1)%ntotvar) )
   allocate( param  (iv%tamdar_varbc(1)%maxpred,iv%tamdar_varbc(1)%ntflight,iv%tamdar_varbc(1)%nphase,iv%tamdar_varbc(1)%ntotvar) )

   call da_get_unit(iunit)
   open(unit=iunit,file='VARBC_TAMDAR.in', form='formatted', status='old')

   do i = 1, 11
      read(iunit, *)
   end do
   do iflt = 1, iv%tamdar_varbc(1)%ntflight
      read(iunit, *) tail_id(iflt), & 
                   ( ifuse(iflt,1:iv%tamdar_varbc(1)%nphase,ivar), &
                     ( param(1:iv%tamdar_varbc(1)%maxpred,iflt,iphase,ivar), iphase=1,iv%tamdar_varbc(1)%nphase ), &
                     ivar=1,iv%tamdar_varbc(1)%ntotvar )
   end do

   close(iunit)
   call da_free_unit(iunit)

   call da_get_unit(iunit)

!   if ( it == max_ext_its ) then
      filename = 'VARBC_TAMDAR.out'
!   else
!      write(unit=filename, fmt='(a,i2.2)') 'VARBC_TAMDAR.out_',it
!   end if

   open(unit=iunit,file=trim(adjustl(filename)),form='formatted',iostat=iost,status='replace')
   if (iost /= 0) then
      message(1)="Cannot open TAMDAR bias correction file "//adjustl(filename)
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   write(iunit, '(A52)')'****************************************************'
   write(iunit, '(A52)')'  PARAMETER TABLE FOR TAMDAR VARBC (MEXICO DOMAIN)  '
   write(iunit, '(A52)')'****************************************************'
   write(iunit, '(A52)')'----------------------------------------------------'
   write(iunit, '(A52)')'Aircrafts currently involved in VarBC   (2016-04-01)'
   write(iunit, '(A52)')'----------------------------------------------------'
   write(iunit, '(I3)') iv%tamdar_varbc(1)%ntflight
   write(iunit, '(A52)')'----------------------------------------------------'
   write(iunit, '(A52)')'ID Ifuse(A/D) Parameters[Ascent(4),Descent(4)]|[T/Q]'
   write(iunit, '(A52)')'[Preditors 1-4: Bias, dp/dt, (dp/dt)^2, Mach number]'
   write(iunit, '(A52)')'----------------------------------------------------'

   do iflt = 1, iv%tamdar_varbc(1)%ntflight
      do ivar = 1, iv%tamdar_varbc(1)%nbcvar
         do iphase = 1, iv%tamdar_varbc(1)%nphase
            param(1:npred,iflt,iphase,iv%tamdar_varbc(ivar)%bcvar_sn) = iv%tamdar_varbc(ivar)%param(1:npred,iflt,iphase)
         end do

         if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 ) param(:,iflt,:,2) = param(:,iflt,:,2)*1000.0
      end do

      write(iunit,fmt=trim(adjustl(iv%tamdar_varbc(1)%fmt_param))) &
            tail_id(iflt), &
          ( ifuse(iflt,1:iv%tamdar_varbc(1)%nphase,ivar), &
            ( param(1:iv%tamdar_varbc(1)%maxpred,iflt,iphase,ivar), iphase=1,iv%tamdar_varbc(1)%nphase ), &
            ivar=1,iv%tamdar_varbc(1)%ntotvar )
   end do

   close(iunit)
   call da_free_unit(iunit)

   deallocate(tail_id, ifuse, param, varbc_param_tl)

   if (rootproc) write(unit=varbc_tamdar_unit,fmt='(//8X,A/)') 'VARBC_TAMDAR is done'

   !----------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_update")

   end subroutine da_varbc_tamdar_update
