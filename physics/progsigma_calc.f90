!>\file progsigma
!! This file contains the subroutine that calculates the prognostic
!! updraft area fraction that is used for closure computations in 
!! saSAS deep and shallow convection.

!>\ingroup samfdeepcnv
!! This subroutine computes a prognostic updraft area fraction
!! used in the closure computations in the samfdeepcnv.f scheme
!>\ingroup samfshalcnv
!! This subroutine computes a prognostic updraft area fracftion
!! used in the closure computations in the samfshalcnv. scheme
!!\section progsigma General Algorithm 
!> @{ 

      subroutine progsigma_calc (im,km,flag_init,flag_restart,flag_deep, &
           del,tmf,qmicro,dbyo1,zdqca,omega_u,zeta,hvap,delt,            &
           qgrs_dsave,q,kbcon1,ktcon,cnvflg,gdx,                         &
           do_ca, ca_closure, ca_entr, ca_trigger, nthresh, ca_deep,     &
           ca_turb,ca_micro,ca_shal,ca_rad,convcount,ca1,ca2,ca3,ca4,    &
           sigmain,sigmaout,sigmab,errmsg,errflg)
!                                                           
!                                                                                                                                             
      use machine,  only : kind_phys
      use funcphys, only : fpvs

      implicit none

!     intent in
      integer, intent(in)  :: im,km,kbcon1(im),ktcon(im)
      real,    intent(in)  :: hvap,delt
      real,    intent(in)  :: qgrs_dsave(im,km), q(im,km),del(im,km),    &
           qmicro(im,km),tmf(im,km),dbyo1(im,km),zdqca(im,km),           &
           omega_u(im,km),zeta(im,km),gdx(im)
      logical, intent(in)  :: flag_init,flag_restart,flag_deep,cnvflg(im)
      real(kind=kind_phys), intent(in) :: nthresh
      real(kind=kind_phys), intent(in) :: ca_deep(im)
      real(kind=kind_phys), intent(out):: ca_turb(im),                   &
           ca_micro(im),ca_rad(im),ca_shal(im),convcount(im),ca1(im),    &
           ca2(im),ca3(im),ca4(im)
      logical, intent(in)  :: do_ca,ca_closure,ca_entr,ca_trigger

      real(kind=kind_phys), intent(in) :: sigmain(im,km)

!     intent out
      real(kind=kind_phys), intent(out) :: sigmaout(im,km)
      real(kind=kind_phys), intent(out) :: sigmab(im)
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

!     Local variables
      integer              :: i,k,km1
      real(kind=kind_phys) :: termA(im),termB(im),termC(im),termD(im),   &
                          mcons(im),zfdqa(im),zform(im,km),              &
                          qadv(im,km),sigmamax(im)                         
                          

      real(kind=kind_phys) :: gcvalmx,ZEPS7,ZZ,ZCVG,mcon,buy2,   &
                          zfdqb,dtdyn,dxlim,rmulacvg,dp,tem,     &
                          alpha,DEN
      integer :: inbu(im,km)  

     !Parameters
     gcvalmx = 0.1
     rmulacvg=10.
     ZEPS7=1.E-11
     km1=km-1
     alpha=7000.

     !Initialization 2D
     do k = 1,km
        do i = 1,im
           sigmaout(i,k)=0.
           inbu(i,k)=0
           zform(i,k)=0. 
        enddo
     enddo
     
     !Initialization 1D
     do i=1,im
         sigmab(i)=0.
         sigmamax(i)=0.95
         termA(i)=0.
         termB(i)=0.
         termC(i)=0.
         termD(i)=0.
         zfdqa(i)=0.
         mcons(i)=0.
      enddo

      !Temporary Initialization output:
       do i = 1,im
         if(flag_deep)then
            !ca_turb(i)=0.
            ca_shal(i)=0.
         endif
         if(.not. flag_deep)then
            ca_rad(i)=0.
            convcount(i)=0.
            ca1(i)=0.
         endif
       enddo

      !Initial computations, place maximum sigmain in sigmab

      do k=2,km
       do i=1,im
          if(flag_init .and. .not. flag_restart)then
             if(cnvflg(i))then
                sigmab(i)=0.03
             endif
          else
             if(cnvflg(i))then
                !if(sigmain(i,k)<1.E-5)then
                !   sigmain(i,k)=0.
                !endif
               if(sigmain(i,k)>sigmab(i))then
                   sigmab(i)=sigmain(i,k)
               endif
             endif
          endif
       enddo
      enddo

      do i=1,im
         if(sigmab(i) < 1.E-5)then !after advection
            sigmab(i)=0.                                                                                                             
         endif
      enddo
           
      !Initial computations, sigmamax
      do i=1,im
         sigmamax(i)=alpha/gdx(i)
         sigmamax(i)=MIN(0.95,sigmamax(i))
      enddo

      !Initial computations, dynamic q-tendency
      do k = 1,km
         do i = 1,im
            if(flag_init .and. .not.flag_restart)then
               qadv(i,k)=0.
            else
               qadv(i,k)=(q(i,k) - qgrs_dsave(i,k))/delt
            endif
         enddo
      enddo
      
      !compute termD "The vertical integral of the latent heat convergence is limited to the                                        
      !buoyant layers with positive moisture convergence (accumulated from the surface).                                                       
      !Lowest level:                                                                                                               
       do i = 1,im
          dp = 1000. * del(i,1)
          mcons(i)=(hvap*(qadv(i,1)+tmf(i,1)+qmicro(i,1))*dp)
       enddo
      !Levels above:
       do k = 2,km1
          do i = 1,im
             dp = 1000. * del(i,k)
             if(cnvflg(i))then
                mcon = (hvap*(qadv(i,k)+tmf(i,k)+qmicro(i,k))*dp)
                buy2 = termD(i)+mcon+mcons(i)
!               Do the integral over buoyant layers with positive mcon acc from surface
                if(k > kbcon1(i) .and. k < ktcon(i) .and. buy2 > 0.)then
                   inbu(i,k)=1
                endif
                inbu(i,k-1)=MAX(inbu(i,k-1),inbu(i,k))
                termD(i) = termD(i) + float(inbu(i,k-1))*mcons(i)
                mcons(i)=mcon
             endif
          enddo
       enddo

       !termA
       do k = 2,km1
          do i = 1,im
             dp = 1000. * del(i,k)
             if(cnvflg(i))then
                tem=(sigmab(i)*zeta(i,k)*float(inbu(i,k))*dbyo1(i,k))*dp
                termA(i)=termA(i)+tem
             endif
          enddo
       enddo

       !termB                                                                                                             
       do k = 2,km1
          do i = 1,im
             dp = 1000. * del(i,k)
             if(cnvflg(i))then
                tem=(dbyo1(i,k)*float(inbu(i,k)))*dp
                termB(i)=termB(i)+tem
             endif
          enddo
       enddo
       
      !termC
       do k = 2,km1
          do i = 1,im
             if(cnvflg(i))then
                dp = 1000. * del(i,k)
                zform(i,k)=-1.0*float(inbu(i,k))*(omega_u(i,k)*delt)
                zfdqb=0.5*((zform(i,k)*zdqca(i,k)))
                termC(i)=termC(i)+(float(inbu(i,k))*   &
                     (zfdqb+zfdqa(i))*hvap*zeta(i,k))
                zfdqa(i)=zfdqb
             endif
         enddo
      enddo

      !sigmab
       do i = 1,im                                                                                                                           
          if(cnvflg(i))then

             DEN=MIN(termC(i)+termB(i),1.E8) !1.E8
             !DEN=MAX(termC(i)+termB(i),1.E7) !1.E7

             ZCVG=termD(i)*delt

             ZZ=MAX(0.0,SIGN(1.0,termA(i)))            &
                  *MAX(0.0,SIGN(1.0,termB(i)))         &
                  *MAX(0.0,SIGN(1.0,termC(i)-ZEPS7))   


             ZCVG=MAX(0.0,ZCVG)
             
             if(flag_init)then
                sigmab(i)=0.03
             else
                sigmab(i)=(ZZ*(termA(i)+ZCVG))/(DEN+(1.0-ZZ))
             endif

             if(sigmab(i)>0.)then
                sigmab(i)=MIN(sigmab(i),sigmamax(i))  
                sigmab(i)=MAX(sigmab(i),0.01)
             endif
             
             if(flag_deep)then
                !ca_turb(i)=ZCVG
                ca_shal(i)=termC(i)
             else
                ca_rad(i)=ZCVG
                ca1(i)=termC(i)
             endif
             !ca3(i)=sigmab(i)

          endif!cnvflg
       enddo

       do k=1,km
          do i=1,im
             if(cnvflg(i))then
                sigmaout(i,k)=sigmab(i)
             endif
          enddo
       enddo
     
     end subroutine progsigma_calc
!> @}                            
!! @} 



