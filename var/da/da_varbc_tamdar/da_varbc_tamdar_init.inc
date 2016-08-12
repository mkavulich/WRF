   subroutine da_varbc_tamdar_init (iv)

   !-------------------------------------------------------!
   !  PURPOSE: Initialize VarBC from VARBC_TAMDAR.in file  !
   !                                                       !
   !  Called from da_solve.inc                             !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   type (iv_type), intent(inout)  :: iv    ! O-B structure.

   integer, parameter             :: num_varbc_var=2,nphase=2,maxpred=4
   integer                        :: i,iunit,id,iflt,iobs,ivar,iphase,nvar,npred,ntflight
   integer                        :: var_sn(num_varbc_var)
   logical                        :: exist_table
   character(len=120)             :: filename

   character(len=5), allocatable  :: tail_id(:)
   integer, allocatable           :: ifuse(:,:,:)
   real, allocatable              :: param(:,:,:,:)

   if (trace_use) call da_trace_entry("da_varbc_tamdar_init")

   !-------------------------------------------------------!

   if (rootproc) then
      write(unit=varbc_tamdar_unit,fmt='(//8X,A/)')      'VARBC_TAMDAR NAMELIST OPTIONS'
      write(unit=varbc_tamdar_unit,fmt='(13X,A,I8)')     'use_varbc_tamdar :         ', use_varbc_tamdar
      write(unit=varbc_tamdar_unit,fmt='(13X,A,4X,2L2)') 'varbc_tamdar_var :         ', varbc_tamdar_var(:)
      write(unit=varbc_tamdar_unit,fmt='(13X,A,4X,I1)')  'varbc_tamdar_bm  :         ', varbc_tamdar_bm
      write(unit=varbc_tamdar_unit,fmt='(13X,A,F8.3)')   'varbc_tamdar_tiv :         ', varbc_tamdar_tiv
      write(unit=varbc_tamdar_unit,fmt='(13X,A,F8.3)')   'varbc_tamdar_qiv :         ', varbc_tamdar_qiv
      write(unit=varbc_tamdar_unit,fmt='(13X,A,F8.2)')   'varbc_tamdar_tune_bc :     ', varbc_tamdar_tune_bc
      write(unit=varbc_tamdar_unit,fmt='(13X,A,F8.2)')   'varbc_tamdar_factor :      ', varbc_tamdar_factor
      write(unit=varbc_tamdar_unit,fmt='(13X,A,I8)')     'varbc_tamdar_nbgerr :      ', varbc_tamdar_nbgerr
      write(unit=varbc_tamdar_unit,fmt='(13X,A,I8)')     'varbc_tamdar_nobsmin :     ', varbc_tamdar_nobsmin
   end if

   nvar = COUNT( varbc_tamdar_var )

   if( nvar == 0 ) then
     if (rootproc) write(unit=varbc_tamdar_unit,fmt='(/A/)')'Set varbc_tamdar_var false for all variables. VARBC_TAMDAR is switched off'
     use_varbc_tamdar = 0
     return
     if (trace_use) call da_trace_exit("da_varbc_tamdar_init")
   end if

   filename = 'VARBC_TAMDAR.in'
   inquire(file=trim(adjustl(filename)), exist=exist_table)

   if ( exist_table ) then

       if (rootproc) write(unit=varbc_tamdar_unit,fmt='(//8X,A/)')'Reading VARBC_TAMDAR.in and Initialization'

       allocate ( iv%tamdar_varbc(nvar) )

       iv%tamdar_varbc(:)%ntotvar = num_varbc_var
       iv%tamdar_varbc(:)%maxpred = maxpred
       iv%tamdar_varbc(:)%nphase  = nphase
       iv%tamdar_varbc(:)%nbcvar  = nvar

       var_sn(:) = 0
       var_sn = PACK( (/(i,i=1,iv%tamdar_varbc(1)%ntotvar)/), mask=(varbc_tamdar_var) )

       nvar = 0
       do ivar = 1, iv%tamdar_varbc(1)%ntotvar
          if( var_sn(ivar) < 1 ) cycle
          nvar = nvar+1
          iv%tamdar_varbc(nvar)%bcvar_sn = var_sn(ivar)
       end do

       call da_get_unit(iunit)
       open(unit=iunit,file=filename, form='formatted', status='old')
       do i = 1, 6; read(iunit,*); end do; read(iunit,*) ntflight

       if (rootproc) write(unit=varbc_tamdar_unit,fmt='(13X,A,I6)')'Number of aircraft in table :',ntflight

       iv%tamdar_varbc(:)%ntflight = ntflight

       npred = 2
       if( varbc_tamdar_bm == 1 ) npred = 2
       if( varbc_tamdar_bm == 2 ) npred = 3
       if( varbc_tamdar_bm == 3 ) npred = 4

       iv%tamdar_varbc(:)%npred = npred

       do ivar = 1, iv%tamdar_varbc(1)%nbcvar
          allocate ( iv%tamdar_varbc(ivar)%tail_id  (ntflight) )
          allocate ( iv%tamdar_varbc(ivar)%ifuse    (ntflight,nphase) )
          allocate ( iv%tamdar_varbc(ivar)%param    (npred,ntflight,nphase) )
          allocate ( iv%tamdar_varbc(ivar)%nobs     (ntflight,nphase+1) )
          allocate ( iv%tamdar_varbc(ivar)%nobs_sum (ntflight,nphase+1) )

          iv%tamdar_varbc(ivar)%tail_id(:)    = 0
          iv%tamdar_varbc(ivar)%ifuse(:,:)    = 0
          iv%tamdar_varbc(ivar)%param(:,:,:)  = 0.
          iv%tamdar_varbc(ivar)%nobs(:,:)     = 0
          iv%tamdar_varbc(ivar)%nobs_sum(:,:) = 0
       end do

       do i = 1, 4 ; read(iunit, *) ; end do

       write(iv%tamdar_varbc(1)%fmt_param,*)'(A5,1X,2I2,8F7.3,1X,2I2,8F7.3)'

       allocate( tail_id(iv%tamdar_varbc(1)%ntflight) )
       allocate( ifuse  (iv%tamdar_varbc(1)%ntflight,nphase,iv%tamdar_varbc(1)%ntotvar) )
       allocate( param  (maxpred,iv%tamdar_varbc(1)%ntflight,nphase,iv%tamdar_varbc(1)%ntotvar) )

       tail_id = ''
       do iflt = 1, ntflight
          read(iunit, fmt=trim(adjustl(iv%tamdar_varbc(1)%fmt_param))) &
               tail_id(iflt),     &
             ( ifuse(iflt,1:nphase,ivar), &
               ( param(1:maxpred,iflt,iphase,ivar), iphase=1,nphase ), &
               ivar=1,iv%tamdar_varbc(1)%ntotvar &
             )

          if( ANY( iv%tamdar_varbc(:)%bcvar_sn == 2 ) ) param(:,iflt,:,2) = param(:,iflt,:,2)/1000.0

          do ivar = 1, iv%tamdar_varbc(1)%nbcvar
             iv%tamdar_varbc(ivar)%ifuse(iflt,1:nphase)   = ifuse(iflt,1:nphase,iv%tamdar_varbc(ivar)%bcvar_sn)
             iv%tamdar_varbc(ivar)%param(1:npred,iflt,1:nphase) = param(1:npred,iflt,1:nphase,iv%tamdar_varbc(ivar)%bcvar_sn)
          end do

          read(tail_id(iflt),'(I5)') iv%tamdar_varbc(1)%tail_id(iflt)

          do iobs = 1, iv%info(tamdar)%nlocal
             read(iv%info(tamdar)%id(iobs),'(I5)') id
             if( iv%tamdar_varbc(1)%tail_id(iflt) .eq. id ) &
                 iv%tamdar_varbc(1)%nobs(iflt,nphase+1) = iv%tamdar_varbc(1)%nobs(iflt,nphase+1) + 1
          end do

          iv%tamdar_varbc(1)%nobs_sum(iflt,nphase+1) = wrf_dm_sum_integer(iv%tamdar_varbc(1)%nobs(iflt,nphase+1))
       end do

       iv%tamdar_varbc(1)%maxobs = MAXVAL( iv%tamdar_varbc(1)%nobs_sum(:,nphase+1) )

       if (rootproc) &
          write(unit=varbc_tamdar_unit,fmt='(/13X,A,I4)') &
               'Max obs number in aircrafts :  ', iv%tamdar_varbc(1)%maxobs

       do ivar = 1, iv%tamdar_varbc(1)%nbcvar
          allocate ( iv%tamdar_varbc(ivar)%obs_sn(iv%tamdar_varbc(1)%maxobs,ntflight,nphase+1) )
          iv%tamdar_varbc(ivar)%obs_sn(:,:,:) = 0
       end do

       do iflt = 1, ntflight
          i=0
          do iobs = 1, iv%info(tamdar)%nlocal
             read(iv%info(tamdar)%id(iobs),'(I5)') id
             if( iv%tamdar_varbc(1)%tail_id(iflt) .eq. id ) then
                 i=i+1
                 iv%tamdar_varbc(1)%obs_sn(i,iflt,nphase+1) = iobs
             end if
          end do
       end do

       do ivar = 1, iv%tamdar_varbc(1)%nbcvar
          allocate ( iv%tamdar_varbc(ivar)%pred  (npred,ntflight,nphase) )
          allocate ( iv%tamdar_varbc(ivar)%bgerr (npred,ntflight,nphase) )
          allocate ( iv%tamdar_varbc(ivar)%index (npred,ntflight,nphase) )
          allocate ( iv%tamdar_varbc(ivar)%vtox  (npred,npred,ntflight,nphase) )

          iv%tamdar_varbc(ivar)%pred(:,:,:)   = 0.
          iv%tamdar_varbc(ivar)%pred(1,:,:)   = 1.
          iv%tamdar_varbc(ivar)%bgerr(:,:,:)  = 0.
          iv%tamdar_varbc(ivar)%index(:,:,:)  = 0
          iv%tamdar_varbc(ivar)%vtox(:,:,:,:) = 0.

          if( ivar > 1 ) then
              iv%tamdar_varbc(ivar)%maxobs               = iv%tamdar_varbc(1)%maxobs
              iv%tamdar_varbc(ivar)%tail_id(:)           = iv%tamdar_varbc(1)%tail_id(:)
              iv%tamdar_varbc(ivar)%nobs(:,nphase+1)     = iv%tamdar_varbc(1)%nobs(:,nphase+1)
              iv%tamdar_varbc(ivar)%nobs_sum(:,nphase+1) = iv%tamdar_varbc(1)%nobs_sum(:,nphase+1)
              iv%tamdar_varbc(ivar)%obs_sn(:,:,nphase+1) = iv%tamdar_varbc(1)%obs_sn(:,:,nphase+1)
          end if
       end do

       close(iunit)
       call da_free_unit(iunit)

   else

       use_varbc_tamdar = 0
       if (rootproc) write(unit=varbc_tamdar_unit,fmt='(/A/)')'Could not find VARBC_TAMDAR.in file. VARBC_TAMDAR is switched off'

   end if

   deallocate(tail_id, ifuse, param)

   !-------------------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_init")

   end subroutine da_varbc_tamdar_init 
