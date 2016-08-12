   subroutine da_varbc_tamdar_pred(iv, be)

   !-------------------------------------------------------!
   !  PURPOSE: Calculate TAMDAR VarBC predictors           !
   !           pred(1) - 1 (Constant)                      !
   !           pred(2) - Ascent/descent rate of aircraft   !
   !           pred(3) - Square of pred(2)                 !
   !           pred(4) - Mach number                       !
   !                                                       !
   !  Reference: Gao (Feng) et al. 2017                    !                                          ! 
   !                                                       !
   !  Called from da_solve.inc                             !
   !                                                       !
   !  HISTORY: 04/11/2016                   fgao@ucar.edu  !
   !-------------------------------------------------------!

   implicit none

   type (iv_type), intent(inout) :: iv
   type (be_type), intent(inout) :: be

   integer                       :: n,iflt,iobs,ipred,ivar,isn,iqc,iphase,npred,nphase
   integer                       :: ja_start(2), size_ja(2)
   real                          :: dpdt

   integer, allocatable          :: tsec(:),obs_sn(:)
   real   , allocatable          :: pp(:), mach(:)

   character(len=24)             :: tchar,start_sec
   character(len=4)              :: ccyy

   if (trace_use) call da_trace_entry("da_varbc_tamdar_pred")

   !------------------------------------------------------!

   if (rootproc) write(unit=varbc_tamdar_unit,fmt='(//8X,A/)')'Calculating predictors'

   npred  = iv%tamdar_varbc(1)%npred
   nphase = iv%tamdar_varbc(ivar)%nphase

   start_sec = "2004-01-01-01-00-00"

   size_ja(2)  = 0
   ja_start(1) = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_js + be%cv%size_jl
   ja_start(2) = 0

   allocate( pp     (iv%tamdar_varbc(1)%maxobs) )
   allocate( tsec   (iv%tamdar_varbc(1)%maxobs) )
   allocate( obs_sn (iv%tamdar_varbc(1)%maxobs) )
   allocate( mach   (iv%tamdar_varbc(1)%maxobs) )

   !------------------------------------------------------!

   do ivar = 1, iv%tamdar_varbc(1)%nbcvar
      if( ivar == 2 ) ja_start(2) = ja_start(1) + size_ja(1)
      do iflt = 1, iv%tamdar_varbc(ivar)%ntflight
         if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,nphase+1) >= varbc_tamdar_nobsmin) then

            pp(:) = 0.
            tsec(:) = 0
            mach(:) = 0
            obs_sn(:) = 0
            n = 0

            if (iv%tamdar_varbc(ivar)%nobs(iflt,nphase+1) > 1) then
               do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,nphase+1)
                  isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,nphase+1)

                  iqc = 0
                  if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 ) iqc = iv%tamdar(isn)%t(1)%qc
                  if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 ) iqc = iv%tamdar(isn)%q(1)%qc

                  if (iv%tamdar(isn)%p(1) > missing_r .and. iqc >= 0) then
                     n = n+1
                     obs_sn(n) = isn
                     pp(n)     = iv%tamdar(isn)%p(1)
                     tchar     = iv%info(tamdar)%date_char(isn)
                     tsec(n)   = int(da_diff_seconds(tchar, start_sec))
                     mach(n)   = 0.8  ! Estimation. Should be read from ob.ascii
                  end if
               end do
            end if

            if (n > 1) then
               do iobs = 2, n
                  dpdt = (pp(iobs)-pp(iobs-1))/100.0/(tsec(iobs)-tsec(iobs-1))

                  if( abs(pp(iobs)-pp(iobs-1)) .lt. 700.0 .or. &
                      abs(tsec(iobs)-tsec(iobs-1)) .ge. 120 ) then
                      if ( dpdt .ne. 0 ) then
                          dpdt = 0.8 * dpdt/abs(dpdt)
                      else
                          dpdt = -0.8
                      end if
                  end if

                  if( abs(tsec(iobs)-tsec(iobs-1)) .lt. 90 .and. abs(tsec(iobs)-tsec(iobs-1)) .gt. 30 &
                      .and. abs(dpdt) .lt. 0.5 .and. abs(dpdt) .gt. 0.0 ) dpdt = 0.8 * dpdt/abs(dpdt)

                  if (iobs == 2) then
                     if (dpdt <= 0.) then
                        iv%tamdar_varbc(ivar)%nobs(iflt,1) = 2
                        iv%tamdar_varbc(ivar)%obs_sn(1:2,iflt,1) = obs_sn(1:2)
                        iv%tamdar_varbc(ivar)%pred(2,iflt,1) = iv%tamdar_varbc(ivar)%pred(2,iflt,1) + dpdt

                        if( varbc_tamdar_bm .gt. 1 ) iv%tamdar_varbc(ivar)%pred(3,iflt,1) = iv%tamdar_varbc(ivar)%pred(3,iflt,1) + dpdt ** 2.0
                        if( varbc_tamdar_bm .gt. 2 ) iv%tamdar_varbc(ivar)%pred(4,iflt,1) = iv%tamdar_varbc(ivar)%pred(4,iflt,1) + mach(1) + mach(2)
                     else
                        iv%tamdar_varbc(ivar)%nobs(iflt,2) = 2
                        iv%tamdar_varbc(ivar)%obs_sn(1:2,iflt,2) = obs_sn(1:2)
                        iv%tamdar_varbc(ivar)%pred(2,iflt,2) = iv%tamdar_varbc(ivar)%pred(2,iflt,2) + dpdt

                        if( varbc_tamdar_bm .gt. 1 ) iv%tamdar_varbc(ivar)%pred(3,iflt,2) = iv%tamdar_varbc(ivar)%pred(3,iflt,2) + dpdt ** 2.0
                        if( varbc_tamdar_bm .gt. 2 ) iv%tamdar_varbc(ivar)%pred(4,iflt,2) = iv%tamdar_varbc(ivar)%pred(4,iflt,2) + mach(1) + mach(2)
                     end if
                  else
                     if (dpdt <= 0.)then
                        iv%tamdar_varbc(ivar)%nobs(iflt,1) = iv%tamdar_varbc(ivar)%nobs(iflt,1) + 1
                        iv%tamdar_varbc(ivar)%obs_sn(iv%tamdar_varbc(ivar)%nobs(iflt,1),iflt,1) = obs_sn(iobs)
                        iv%tamdar_varbc(ivar)%pred(2,iflt,1) = iv%tamdar_varbc(ivar)%pred(2,iflt,1) + dpdt

                        if( varbc_tamdar_bm .gt. 1 ) iv%tamdar_varbc(ivar)%pred(3,iflt,1) = iv%tamdar_varbc(ivar)%pred(3,iflt,1) + dpdt ** 2.0
                        if( varbc_tamdar_bm .gt. 2 ) iv%tamdar_varbc(ivar)%pred(4,iflt,1) = iv%tamdar_varbc(ivar)%pred(4,iflt,1) + mach(iobs)

                     else
                        iv%tamdar_varbc(ivar)%nobs(iflt,2) = iv%tamdar_varbc(ivar)%nobs(iflt,2) + 1
                        iv%tamdar_varbc(ivar)%obs_sn(iv%tamdar_varbc(ivar)%nobs(iflt,2),iflt,2) = obs_sn(iobs)
                        iv%tamdar_varbc(ivar)%pred(2,iflt,2) = iv%tamdar_varbc(ivar)%pred(2,iflt,2) + dpdt

                        if( varbc_tamdar_bm .gt. 1 ) iv%tamdar_varbc(ivar)%pred(3,iflt,2) = iv%tamdar_varbc(ivar)%pred(3,iflt,2) + dpdt ** 2.0
                        if( varbc_tamdar_bm .gt. 2 ) iv%tamdar_varbc(ivar)%pred(4,iflt,2) = iv%tamdar_varbc(ivar)%pred(4,iflt,2) + mach(iobs)
                     end if
                  end if
               end do
            end if

            do iphase = 1, nphase
               iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) = wrf_dm_sum_integer(iv%tamdar_varbc(ivar)%nobs(iflt,iphase))
               if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) > 0) then
                  do ipred = 2, npred
                     iv%tamdar_varbc(ivar)%pred(ipred,iflt,iphase) = wrf_dm_sum_real(iv%tamdar_varbc(ivar)%pred(ipred,iflt,iphase)) / &
                                                                     iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase)
                  end do
               end if

               if (iv%tamdar_varbc(ivar)%nobs_sum(iflt,iphase) >= varbc_tamdar_nobsmin) then 
                  if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0 .and. iv%tamdar_varbc(ivar)%ifuse(iflt,iphase) > 0 ) then
                     do ipred = 1, npred
                        size_ja(ivar) = size_ja(ivar) + 1
                        iv%tamdar_varbc(ivar)%index(ipred,iflt,iphase) = ja_start(ivar) + size_ja(ivar)
                        iv%tamdar_varbc(ivar)%vtox(ipred,ipred,iflt,iphase) = 1.0
                     end do
                  end if
               else
                  if (iv%tamdar_varbc(ivar)%nobs(iflt,iphase) > 0) then
                     do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,iphase)
                        isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,iphase)

                        if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(isn)%t(1)%qc >= 0 ) iv%tamdar(isn)%t(1)%qc = -55
                        if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(isn)%q(1)%qc >= 0 ) iv%tamdar(isn)%q(1)%qc = -55
                     end do
                  end if
               end if

            end do

         else

            iv%tamdar_varbc(ivar)%ifuse(iflt,:) = 0

            if (iv%tamdar_varbc(ivar)%nobs(iflt,nphase+1) > 0) then
               do iobs = 1, iv%tamdar_varbc(ivar)%nobs(iflt,nphase+1)
                  isn = iv%tamdar_varbc(ivar)%obs_sn(iobs,iflt,nphase+1)

                  if( iv%tamdar_varbc(ivar)%bcvar_sn == 1 .and. iv%tamdar(isn)%t(1)%qc >= 0 ) iv%tamdar(isn)%t(1)%qc = -55
                  if( iv%tamdar_varbc(ivar)%bcvar_sn == 2 .and. iv%tamdar(isn)%q(1)%qc >= 0 ) iv%tamdar(isn)%q(1)%qc = -55
               end do
            end if

            if (rootproc .and. iv%tamdar_varbc(ivar)%nobs_sum(iflt,nphase+1) > 0) &
               write(unit=varbc_tamdar_unit,fmt='(13X,A,I5,A,I5)') &
                    'Measurement from ', iv%tamdar_varbc(ivar)%tail_id(iflt),' is insufficient: ',iv%tamdar_varbc(ivar)%nobs_sum(iflt,nphase+1)

         end if
      end do
   end do

   be % cv % size_ja(1) = size_ja(1)
   be % cv % size_ja(2) = size_ja(2)

   if (rootproc .and. ANY(iv%tamdar_varbc(1)%nobs_sum(:,nphase+1) >= varbc_tamdar_nobsmin)) &
      write(unit=varbc_tamdar_unit,fmt='(/13X,A5,A7,A16,A10)') &
           ' ID  ','  Ifuse ',' Ascent/Descent',' Nsum(a/d) '

   do iflt = 1, iv%tamdar_varbc(1)%ntflight
      if (rootproc .and. iv%tamdar_varbc(1)%nobs_sum(iflt,nphase+1) >= varbc_tamdar_nobsmin) then
         write(unit=varbc_tamdar_unit,fmt='(13X,I5,1X,2(1X,I2),2F8.3,2I5,1X,2(1X,I2),2F8.3,2I5)') &
               iv%tamdar_varbc(1)%tail_id(iflt), &
             ( iv%tamdar_varbc(ivar)%ifuse(iflt,1:nphase), &
               iv%tamdar_varbc(ivar)%pred(2,iflt,1:nphase), &
               iv%tamdar_varbc(ivar)%nobs_sum(iflt,1:nphase), ivar=1,iv%tamdar_varbc(1)%nbcvar )
      end if
   end do

   deallocate(pp,tsec,obs_sn,mach)

   !------------------------------------------------------!

   if (trace_use) call da_trace_exit("da_varbc_tamdar_pred")

   end subroutine da_varbc_tamdar_pred
