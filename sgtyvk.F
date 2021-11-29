************************************************************************
*     -----------------
      SUBROUTINE SGTYVK
*     -----------------
*
*     (Purpose)
*       Check if the photon is Reflected on TYVEK or Not
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1995.04.27 ; First version by Y.Koshio
*       2008.01.17 ; Implement new tyvek reflectivity model R.Wendell (rvw)
*                    The old model is left in for reference and is switchable
*                    using the "use_new_model" below                    
*                    Please also see, sgtyvk_util.F and sgtyvk_odseg.F
*
************************************************************************

#include "geant321/gconst.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geotnk.h"      
#include "dsparma.h"

      real photon_pos(3), photon_dir(3)
      real vrefl(3),tyvknorm(3)
      real threfl,phrefl
      PARAMETER (dstep=10.)
      integer ierr
      real thnorm,phnorm
      real refldir(3)
      real dir1
      logical do_2pi
      real inc_dir(3)
      real inc_dir2(3)
      real e1(3),e2(3),e3(3)
      real Ie1(3),Ie2(3),Ie3(3) ! incident basis
      real refl(3), A(3), B(3), inner, mA, mB
      integer angbin
      real angle
      real gausprob
      real thI, thR, phiI, phiR
      real myDummy
! Specify our options and such....
      integer use_new_model
      real odrad_mid
      real odz_mid

      real angle_bins(10)
      data angle_bins/0.,10.,20.,30.,40.,50.,60.,70.,80.,90./ 

      real gauss_prob_new(10)
!  the new model A.Chavarria 2007, Gausian/ (Gauss+Lambert)
!     binned in degrees, 0, 10, ...80, 90. However, 90 degrees was not measured
!     data gauss_prob_new/0.48,0.56,0.46,0.27,0.34,0.35,0.39,0.49,0.56,0.54/ !air
      data gauss_prob_new/0.60,0.78,0.66,0.59,0.52,0.57,0.63,0.67,0.74,0.54/ !Water


!  the old, double isotropic model
c Roughly 1/(cos(theta)**1.75 dependence) for prob to use Gaussian vs Lambert
c Can do something more sophisticated later..
      real gauss_prob_old(10)
      data gauss_prob_old/0.02,0.021,0.022,0.026,0.03,0.038,0.053,0.083,0.16,0.54/


! Specify our printing options 0: old model, 1: new model , rvw 
      use_new_model = 1

c-J.G. added following r and z
      r = sqrt((vect(1)**2)+vect(2)**2)
      z = vect(3)

      ! midpoint radiues between inner tank and outer tank
      odrad_mid = ((rintk +rmed) + rtktk )/2.
      odz_mid   = ((zpintk+zmed) + zptktk)/2.

      alamda = 1237.8e-18 / gekin

        rref = 0.0
        ! top and bottom tyvek surfaces 
        ! for the top end cap 
         if (z.ge.1855. .and.  r .le. (rintk+rmed) ) then
            rref = tyvksg(alamda)*dstyvkrtop(1)
            ! top surface
            if( z .gt. odz_mid ) then 
               rref = tyvksg(alamda)*dstyvkrtop(2)
            endif
         endif

        ! top and bottom tyvek surfaces 
        ! for the bottom end cap
         if (z.le.-1855.  .and.  r .le. (rintk+rmed) ) then
            rref = tyvksg(alamda)*dstyvkrbot(1)
            if( z .lt. -1. * odz_mid ) then 
                !bottom surface 
                rref = tyvksg(alamda)*dstyvkrbot(2)
            endif
         endif


         ! inner surface of the barrell section
         if ((z.ge.-1855.).and.(z.le.1855.)) then
            rref = tyvksg(alamda)*dstyvkrwal(1)
            ! outer surface of the barrell mid section
            if( r .gt. odrad_mid ) then 
                rref = tyvksg(alamda)*dstyvkrwal(2)
            endif
         endif

         ! outer surface of the entire barrell section as 
         ! well as top/bottom surface outside of the segmentation
         if ( r .gt. (rintk+rmed) ) then
              rref = tyvksg(alamda)*dstyvkrwal(2)
         endif

      if(  rref .gt. 1.00 .or. (rref.lt.1.0e-5.and.tyvksg(alamda).gt. 1.0e-5) ) then 
        write(6,*) 'rref error: ', rref,
     +              r, z,  (rintk+rmed) ,  tyvksg(alamda),dstyvkrwal(2) 
        stop 
      endif 

C      rnref= rndm()
      CALL GRNDM(RNREF,1)

*     Check if the photon is reflected?
c      print*,'REF',rnref,rref
      if(rnref.le.rref) then
         ngphot=1
         do 10 i=1,3
            xphot(i,ngphot) = vect(i)-dstep*vect(i+3)
            photon_pos(i) = xphot(i,ngphot)
            photon_dir(i) = vect(i+3)
 10      continue

         call get_tyvknorm(photon_pos,photon_dir,tyvknorm,ierr)

c Current meaning of dstyvkm parameter is Lambert/Gaussian fraction
c   This could change in the future

c Get basis for Tyvek normal

         call vbasis(tyvknorm,e1,e2,e3)


         ! old apdetsim routine for SK-I pathed in 20071221 rvw
        if ( SK_GEOMETRY.EQ.SK_I .and. use_new_model .eq. 0) THEN
 
          CALL GRNDM(RDUMMY,1)
          ph=RDUMMY*twopi
          CALL GRNDM(RDUMMY,1)
          th=RDUMMY*twopi
        
        else !non SK-I geometry  

         CALL GRNDM(RDUMMY,1)

! Allow component of 2pi isotropic for old model
! the new model has no isotropic component and will
! automatically enter the Gaussian vs. Lambert testing
         if ( (rdummy.lt.DSTYVKM .and. use_new_model .eq. 0)
     +       .or.
     +        use_new_model .eq. 1 
     +      ) then    

c Use Lambert + Gaussian reflection
c Ref: Hasenbalg and Ravignani

c Transform photon direction

            call vtransform(photon_dir,e1,e2,e3,inc_dir)
               
c Probability of Lambert vs Gaussian (a function of incident angle)

            call GRNDM(RDUMMY,1)
c Find the probability bin according to the angle

            if( use_new_model .eq. 1 ) then
              if( inc_dir(3)  .gt. 1.00 ) then
                 inc_dir(3) = 1.00
              elseif( inc_dir(3) .lt. -1.00 )then
                 inc_dir(3)  = -1.00
              endif

              angle = acos(-inc_dir(3))*180./pi
              angbin = nint(angle/10.) + 1  
 
              if( angbin .ge. 10 ) angbin = 9
              ! fortran will round, so move the interpolation
              if( angle - angle_bins(angbin) .lt. 0 ) then
                  angbin = angbin - 1
              endif

             ! Lets do some interpolation
             ! should only be a minor change
              gausprob = (gauss_prob_new( angbin+1)-gauss_prob_new(angbin))
     +                 / ( angle_bins(angbin+1)-angle_bins(angbin) ) 
              gausprob = gausprob * (angle-angle_bins(angbin)) + 
     +                   gauss_prob_new(angbin) 

            else
                angbin = nint(acos(-inc_dir(3))*2./pi*10.)
            endif  ! end of test on which model is being used


            if ( (angbin.gt.10.and. use_new_model.eq.0).or. angbin.le.0 ) angbin=1
            if (  angbin.gt.10 .and. use_new_model.eq.1   ) angbin=10  ! A.Chavarria didn't measure incidence > 80 degrees

            ! 20091222 introduce scaling factor via dstyvkm parameter
            ! aids in the tuning of ODrmean and related variables.
            ! introduced with the idea in mind that the SK water 
            ! is likely more pure than that used for the reflectivity measurement
            if ( (rdummy.gt.gauss_prob_old (angbin).and. use_new_model.eq.0 )
     +           .or.
     +           (rdummy.gt.dstyvkm*gausprob .and. use_new_model.eq.1 )
     +          ) then 

               call sglambert(vrefl)
           
            else

c Gaussian reflection

               call sggauss_refl(inc_dir,vrefl)

            endif

            
         else   ! End of Lambert/Gaussian case only if using the old model

            call sgisotropic(vrefl)

         endif
         
* Get reflected photon direction 
                                                   ! inc_dir is in Tyvek normal coordinates 
                                                   ! vrefl is in Tyvek normal coordinates
         call vuntransform(vrefl,e1,e2,e3,refldir) ! here refldir is in SK coordinates


* Theta and phi of reflected photon in SK frame
* Here th=0 is straight up (th is angle wrt z axis)  

         call sgget_angles(refldir,th,ph)

        endif ! end of SK-I geometry and old model test

         cosph=cos(ph)
         sinph=sin(ph)
         costh=cos(th)
         sinth=sin(th)

         xphot(4,ngphot)=sinth*cosph
         xphot(5,ngphot)=sinth*sinph
         xphot(6,ngphot)=costh

c         dotprod =    xphot(4,ngphot)*tyvknorm(1) +
c     &         xphot(5,ngphot)*tyvknorm(2) +
c     &         xphot(6,ngphot)*tyvknorm(3)
c         write(6,99)'sgtyvk',dotprod, acos(-inc_dir(3)),
c     &            tyvknorm(1),tyvknorm(2),tyvknorm(3),
c     &               xphot(4,ngphot),xphot(5,ngphot),xphot(6,ngphot)
c 99      format(a6,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4)
c         write(6,*)'---------------------------------------'

         xphot(8,ngphot)=costh*cosph
         xphot(9,ngphot)=costh*sinph
         xphot(10,ngphot)=-sinth

         xphot(7,ngphot)=gekin
         xphot(11,ngphot)=tofg
         call gskpho(0)
         ngphot=0
         
c      print*,(xphot(ii,1),ii=1,6)
      endif  ! End of reflected case


*      close(22)


      return
      end
   
************************************************************************

      subroutine get_tyvknorm(pos,dir,tyvknorm,ierr)

      implicit none

#include "geotnk.h"      

* Return normal to nearest Tyvek wall for given input position
* There might be a better way to do this using Geant infrastructure...

      real pos(3),dir(3),tyvknorm(3)
      real norm
      integer i
      integer norf,ic,in
      real hpos(3)
      integer ierr

* Find the intersection with the cylinder where the photon is headed


* First check the inner cylinder
      norf=0
      call crstnksk(pos,dir,zpintk+zmed,rintk+rmed,norf,ic,in,hpos)

      if (ic.ne.0) then

* The photon hits the inner cylinder

        if (ic.eq.1) then

* Crosses inner barrel            

            tyvknorm(1) = pos(1) 
            tyvknorm(2) = pos(2)
            tyvknorm(3) = 0.

        else if (ic.eq.2) then

* Crosses inner top

            tyvknorm(1) = 0.
            tyvknorm(2) = 0.
            tyvknorm(3) = 1.

         else if (ic.eq.3) then

* Crosses inner bottom

            tyvknorm(1) = 0.
            tyvknorm(2) = 0.
            tyvknorm(3) = -1.

         endif

      else

* Photon did not hit inner cylinder...look for intersection with outer cyl

         call crstnksk(pos,dir,zptktk,rtktk,norf,ic,in,hpos)

         if (ic.eq.1) then

* Crosses outer barrel            

            tyvknorm(1) = -pos(1)
            tyvknorm(2) = -pos(2)
            tyvknorm(3) = 0.

         else if (ic.eq.2) then

* Crosses outer top

            tyvknorm(1) = 0.
            tyvknorm(2) = 0.
            tyvknorm(3) = -1.

         else if (ic.eq.3) then

* Crosses outer bottom

            tyvknorm(1) = 0.
            tyvknorm(2) = 0.
            tyvknorm(3) = 1.

         else

            write(6,*)'Error!  No intersection'
            ierr = 1
            return

        endif
      endif

* Normalize

      norm = sqrt (tyvknorm(1)**2 + tyvknorm(2)**2 + tyvknorm(3)**2)

      do i=1,3
         tyvknorm(i) = tyvknorm(i)/norm
      enddo
      ierr = 0

      return 
      end


