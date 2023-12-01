!>\file GFS_rad_time_vary.fv3.F90
#ifdef USE_MKL_MERSENNE
include 'mkl_vsl.f90'
#endif
!!  Contains code related to GFS radiation suite setup (radiation part of time_vary_step)
   module GFS_rad_time_vary

      implicit none

      private

      public GFS_rad_time_vary_timestep_init

      contains

!>\defgroup mod_GFS_rad_time_vary GFS Radiation Time Update
!! This module contains code related to GFS radiation setup.
!> @{
!> \section arg_table_GFS_rad_time_vary_timestep_init Argument Table
!! \htmlinclude GFS_rad_time_vary_timestep_init.html
!!
      subroutine GFS_rad_time_vary_timestep_init (lrseeds, rseeds,                     &
              lslwr, lsswr, isubc_lw, isubc_sw, icsdsw, icsdlw, cnx, cny, isc, jsc,    &
              imap, jmap, sec, kdt, imp_physics, imp_physics_zhao_carr, ipsd0, ipsdlim,&
              ps_2delt, ps_1delt, t_2delt, t_1delt, qv_2delt, qv_1delt, t, qv, ps,     &
              errmsg, errflg)

         use mersenne_twister,          only: random_setseed, random_index, random_stat
         use machine,                   only: kind_phys
         use radcons,                   only: qmin, con_100
#ifdef USE_MKL_MERSENNE
         use mkl_vsl_type
         use mkl_vsl
#endif
         implicit none

         ! Interface variables
         logical,                intent(in)    :: lrseeds
         integer,                intent(in)    :: rseeds(:,:)
         integer,                intent(in)    :: isubc_lw, isubc_sw, cnx, cny, isc, jsc, kdt
         integer,                intent(in)    :: imp_physics, imp_physics_zhao_carr, ipsd0, ipsdlim
         logical,                intent(in)    :: lslwr, lsswr
         integer,                intent(inout) :: icsdsw(:), icsdlw(:)
         integer,                intent(in)    :: imap(:), jmap(:)
         real(kind_phys),        intent(in)    :: sec
         real(kind_phys),        intent(inout) :: ps_2delt(:)
         real(kind_phys),        intent(inout) :: ps_1delt(:)
         real(kind_phys),        intent(inout) :: t_2delt(:,:)
         real(kind_phys),        intent(inout) :: t_1delt(:,:)
         real(kind_phys),        intent(inout) :: qv_2delt(:,:)
         real(kind_phys),        intent(inout) :: qv_1delt(:,:)
         real(kind_phys),        intent(in)    :: t(:,:), qv(:,:), ps(:)
         character(len=*),       intent(out)   :: errmsg
         integer,                intent(out)   :: errflg

         ! Local variables
         type (random_stat) :: stat
         integer :: ix, j, i, ipseed
         integer :: numrdm(cnx*cny*2)
#ifdef USE_MKL_MERSENNE
         real(kind=4) :: subDomain(size(jmap))
         type (VSL_STREAM_STATE) :: stream
         integer(kind=8) SkipCount
         integer(kind=4) err
#endif
         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (lsswr .or. lslwr) then

           !--- call to GFS_radupdate_timestep_init is now in GFS_rrtmg_setup_timestep_init

           !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
           if ((isubc_lw==2) .or. (isubc_sw==2)) then
             !NRL If random seeds supplied by NEPTUNE
             if(lrseeds) then
               do ix=1,size(jmap)
                 icsdsw(ix) = rseeds(ix,1)
                 icsdlw(ix) = rseeds(ix,2)
               enddo
             else
               ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
#ifdef USE_MKL_MERSENNE
               ! Create a new stream using the new seed
               err=vslnewstream( stream, VSL_BRNG_MT19937,  ipseed )

               ! Skip ahead to local ranks section of the output.
               ! Avoids having every rank compute random numbers for the whole grid.
               SkipCount = int((imap(1)+isc-1+(jmap(1)+jsc-2)*cnx)-1,8)
               err=vslskipaheadstream( stream,  SkipCount)

               ! Compute local ranks section of the output
               err=vsrnguniform( VSL_RNG_METHOD_UNIFORM_STD, stream, size(icsdsw), subDomain, real(0.0,4), real(ipsdlim,4) )
               icsdsw=int(subDomain)

               !do ix=1,size(jmap)
               !  j = jmap(ix)
               !  i = imap(ix)
               !  write(6,'("RNG: ",i12)') (i+isc-1 + (j+jsc-2)*cnx)
               !enddo

               ! Skip ahead again
               err=vslskipaheadstream( stream,  int(cnx*cny-size(icsdsw),8))
               err=vsrnguniform( VSL_RNG_METHOD_UNIFORM_STD, stream, size(icsdlw), subDomain, real(0.0,4), real(ipsdlim,4) )
               icsdlw=int(subDomain)
               !write(6,'("RNG: ",4i12)') ipseed, SkipCount, cnx*cny, size(icsdlw)
               !write(6,'("RNG: ",9i12)') ipseed, size(jmap), isc, jsc, imap(1), jmap(1), imap(size(jmap)), jmap(size(jmap)), cnx*cny

               ! Delete the stream
               err=vsldeletestream( stream )
#else
               call random_setseed (ipseed, stat)
               call random_index (ipsdlim, numrdm, stat)

               do ix=1,size(jmap)
                 j = jmap(ix)
                 i = imap(ix)
                 !--- for testing purposes, replace numrdm with '100'
                 icsdsw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx)
                 icsdlw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx + cnx*cny)
               enddo
#endif
             end if !lrseeds
           endif  ! isubc_lw and isubc_sw

           if (imp_physics == imp_physics_zhao_carr) then
             if (kdt == 1) then
               t_2delt  = t
               t_1delt  = t
               qv_2delt = qv
               qv_1delt = qv
               ps_2delt = ps
               ps_1delt = ps
             endif
           endif

         endif

      end subroutine GFS_rad_time_vary_timestep_init
!> @}

   end module GFS_rad_time_vary
