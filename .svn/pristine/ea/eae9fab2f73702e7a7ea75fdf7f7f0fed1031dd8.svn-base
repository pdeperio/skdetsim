************************************************************************
*     -----------------
      SUBROUTINE SGTYVK_ODSEGM
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
*       1995.04.27 ; First version by Y.Koshio,
*       
*       2005.08    ;  F.Dufour, Modified for SK-III for OD segmentation
*       2008.01.17 ;  Add in new tyvek model of reflectivity as measured by 
*                  ;  A. Chavarria 2007. More details in sgtyvk.F R.Wendell (rvw) 
*       2020.05    ;  T. Wester updated reflection calculation to be
*                     consistent with sgtyvk.F
*
************************************************************************

#include "geant321/gconst.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "skdetsim_parm.h"
#include "skhead.h"
#include "dsparma.h"

      real rnref,rntrans
     
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
      real e1(3),e2(3),e3(3)
      integer angbin
      real gausprob

      real prob_transmitted

c Roughly 1/(cos(theta)**1.75 dependence) for prob to use Gaussian vs Lambert
c Can do something more sophisticated later..
      real gauss_prob(10)
! the old model:
!      data gauss_prob/0.02,0.021,0.022,0.026,0.03,0.038,0.053,0.083,0.16,0.54/
! A.Chavarria model, Gaus/(Gaus+Lambert) fraction binned in degrees
      data gauss_prob/0.60,0.78,0.66,0.59,0.52,0.57,0.63,0.67,0.74,0.54/ !Water
      real angle_bins(10)
      data angle_bins/0.,10.,20.,30.,40.,50.,60.,70.,80.,90./ 

c-J.G. added following r and z
      r = sqrt((vect(1)**2)+vect(2)**2)
      z = vect(3)

      alamda = 1237.8e-18 / gekin


c-J.G.      rref = tyvksg(alamda)*DSTYVKR
      if (z.ge.1855.) then
         rref = tyvksg(alamda)*DSTYVKRTOP(1)*DSSEGRTOP
      endif
      if (z.le.-1855.) then
         rref = tyvksg(alamda)*DSTYVKRBOT(1)*DSSEGRBOT
      endif
c Not relevant for segmentation Tyvek
      if ((z.gt.-1855.).and.(z.lt.1855.)) then
         rref = tyvksg(alamda)*DSTYVKRWAL(1)
      endif
      
C      rnref= rndm()
      CALL GRNDM(RNREF,1)
      CALL GRNDM(RNTRANS,1)


      if (z.gt.0) then
         prob_transmitted = PROBTT
      else
         prob_transmitted = PROBTB
      endif

      if(RNTRANS.gt.prob_transmitted) then 

c photon is not transmitted:
         
*     Check if the photon is reflected?
c     write(6,*)'rnref,rref',rnref,rref
         if(rnref.le.rref) then
            ngphot=1
            do 10 i=1,3
               xphot(i,ngphot) = vect(i)-dstep*vect(i+3)
               photon_pos(i) = xphot(i,ngphot)
               photon_dir(i) = vect(i+3)
 10         continue
         
            call get_tyvknorm_odsegm(photon_pos,photon_dir,tyvknorm,ierr)
            if (ierr.ne.0)return
         
c     Current meaning of dstyvkm parameter is Lambert/Gaussian fraction
c     This could change in the future
         
c     Get basis for Tyvek normal
            
            call vbasis(tyvknorm,e1,e2,e3)
         
            CALL GRNDM(RDUMMY,1)
         
c     Allow component of 2pi isotropic
         
!            if (rdummy.lt.DSTYVKM) then

c     Use Lambert + Gaussian reflection

c     Ref: Hasenbalg and Ravignani
            
c     Transform photon direction
            
               call vtransform(photon_dir,e1,e2,e3,inc_dir)
            
c     Probability of Lambert vs Gaussian (a function of incident angle)
            
               call GRNDM(RDUMMY,1)

c Find the probability bin according to the angle

!  new model is binned in degrees, rvw
                angbin = nint(acos(-inc_dir(3))*180./pi/10.) + 1
!  old model 
!               angbin = nint(acos(-inc_dir(3))*2./pi*10.)
               
               if (angbin.gt.9  ) angbin=9   ! A.Chavarria didn't measure incidence > 80 degrees
               if (angbin.gt.10.or.angbin.le.0) angbin=1

               if (SK_GEOMETRY.ge.SK_V) then
c                  Follow sgtyvk.F to be consistent
                   ! Lets do some interpolation
                   ! should only be a minor change
                    gausprob = (gauss_prob( angbin+1)-gauss_prob(angbin))
     +                       / ( angle_bins(angbin+1)-angle_bins(angbin))
                    gausprob = gausprob * (angle-angle_bins(angbin)) + 
     +                         gauss_prob(angbin) 

                   if (rdummy.gt.dstyvkm*gausprob) then 
c                      Lambert reflection, independent of incident angle
                       call sglambert(vrefl)
                   else
c                     Gaussian reflection
                      call sggauss_refl(inc_dir,vrefl)
                   endif
               else
c                  <= SK I-IV, use simplified segmentation tyvek model                   
                   if (rdummy.gt.gauss_prob(angbin)) then 
c                      Lambert reflection, independent of incident angle
                       call sglambert(vrefl)
                   else
c                      Gaussian reflection
                       call sggauss_refl(inc_dir,vrefl)
                   endif
               endif



!! this is remant of old model               
!            else                ! End of Lambert/Gaussian case
!               
!               call sgisotropic(vrefl)
!               
!         endif
         
*     Get reflected photon direction 
         
         call vuntransform(vrefl,e1,e2,e3,refldir)
         
            
*     Theta and phi of reflected photon in SK frame
*     Here th=0 is straight up (th is angle wrt z axis)
            
            call sgget_angles(refldir,th,ph)
            
            cosph=cos(ph)
            sinph=sin(ph)
            costh=cos(th)
            sinth=sin(th)
            
            xphot(4,ngphot)=sinth*cosph
            xphot(5,ngphot)=sinth*sinph
            xphot(6,ngphot)=costh
            
c     dotprod =    xphot(4,ngphot)*tyvknorm(1) +
c     &         xphot(5,ngphot)*tyvknorm(2) +
c     &         xphot(6,ngphot)*tyvknorm(3)
c     write(6,99)'sgtyvk',dotprod, acos(-inc_dir(3)),
c     &            tyvknorm(1),tyvknorm(2),tyvknorm(3),
c     &               xphot(4,ngphot),xphot(5,ngphot),xphot(6,ngphot)
c     99      format(a6,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4)
c     write(6,*)'---------------------------------------'
            
            xphot(8,ngphot)=costh*cosph
            xphot(9,ngphot)=costh*sinph
            xphot(10,ngphot)=-sinth
            
            xphot(7,ngphot)=gekin
            xphot(11,ngphot)=tofg
c            write(*,*)"In reflecetd"
c            do i=1,11
c               write(*,*) "xphot",i, xphot(i,ngphot)
c            enddo
           
            call gskpho(0)
            ngphot=0
            
         endif ! End of reflected case
         
      else         ! Photon is transmitted
         
c         print*,"Transmitted light"
c     
         
c     IJKL = IJKL + 1
c     IF ((IJKL .LE. 10).OR.(MOD(IJKL,1000) .EQ. 0)) THEN
c     WRITE(*,*) IJKL," GRNDM RDUMMY =", RDUMMY
c     write(*,*)"before probt", probt
c     END IF
         
c     if (rdummy.lt.PROBT) then
c            write(*,*)"before probt", probt,rdummy
         
         
         ngphot=1
         do 11 i=1,3
            xphot(i,ngphot) = vect(i)+dstep*vect(i+3)
            photon_pos(i) = xphot(i,ngphot)
            photon_dir(i) = vect(i+3)
 11      continue
            
            
         call sgget_angles(photon_dir,th,ph)
            
         cosph=cos(ph)
         sinph=sin(ph)
         costh=cos(th)
         sinth=sin(th)
            
         xphot(4,ngphot)=sinth*cosph
         xphot(5,ngphot)=sinth*sinph
         xphot(6,ngphot)=costh
         
         xphot(8,ngphot)=costh*cosph
         xphot(9,ngphot)=costh*sinph
         xphot(10,ngphot)=-sinth
         
         xphot(7,ngphot)=gekin
         xphot(11,ngphot)=tofg
         
c         write(*,*)"In transmitted"
c         do i=1,11
c            write(*,*) "xphot",i, xphot(i,ngphot)    
c         enddo
         
         call gskpho(0)
         
         
         
         ngphot=0
         
         
      endif
      
      
c     print*,(xphot(ii,1),ii=1,6)
                           
      
      return

      end




**********************************************************

  
      subroutine get_tyvknorm_odsegm(pos,dir,tyvknorm,ierr)

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

      ierr=0
* Find the intersection with the cylinder where the photon is headed

* For OD segmentation: cylinders to check depend on r

      if (sqrt(pos(1)**2+pos(2)**2).gt.rintk+rmed) then

* Outside the Tyvek barrier

*  First check the inner tall cylinder
         norf=0
         call crstnksk(pos,dir,zptktk,rintk+rmed,norf,ic,in,hpos)
      
         if (ic.ne.0) then

* The photon hits the inner cylinder

            if (ic.eq.1) then

* Crosses inner barrel            
               
               tyvknorm(1) = pos(1) 
               tyvknorm(2) = pos(2)
               tyvknorm(3) = 0.

            else if (ic.eq.2) then

* Crosses inner top -- should never happen!

               tyvknorm(1) = 0.
               tyvknorm(2) = 0.
               tyvknorm(3) = 1.

               write(6,*)'Error: crosses inner tall top from outside ',pos,dir
               ierr=1

            else if (ic.eq.3) then

* Crosses inner bottom

               tyvknorm(1) = 0.
               tyvknorm(2) = 0.
               tyvknorm(3) = -1.

               write(6,*)'Error: crosses inner tall bottom from outside ',pos,dir
               ierr=1
            endif

         else

* Photon did not hit inner cylinder...look for intersection with outer cyl

            call crstnksk(pos,dir,zptktk,rtktk,norf,ic,in,hpos)

            if (ic.eq.1) then
               
*     Crosses outer barrel            
               
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
               
               write(6,*)'Error in sgtyvk!  No intersection: pos,dir ',pos,dir
               ierr = 1
               return

            endif  ! End of ic case
         endif ! End of did not hit inner case

      else   ! Inside the barrier on the top


*  First check the inner cylinder (top of inner detector)
         norf=0
         call crstnksk(pos,dir,zpintk+zmed,rintk+rmed,norf,ic,in,hpos)
      
         if (ic.ne.0) then

* The photon hits the inner cylinder

            if (ic.eq.1) then

* Crosses inner barrel -- should never happen!
               
               tyvknorm(1) = pos(1) 
               tyvknorm(2) = pos(2)
               tyvknorm(3) = 0.

               write(6,*)
     &           'Error: crosses inner cylinder from inside barrier ',
     &                   sqrt(pos(1)**2+pos(2)**2),pos,dir
               ierr=1

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

* Photon did not hit inner cylinder...look for intersection with barrier

            call crstnksk(pos,dir,zptktk,rintk+rmed,norf,ic,in,hpos)

            if (ic.eq.1) then
               
*     Crosses outer barrel            
               
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
               
               write(6,*)'Error in sgtyvk!  No intersection: pos,dir ',pos,dir
               ierr = 1
               return

            endif  ! End of ic case
         endif ! End of did not hit inner case



      endif

* Normalize

      norm = sqrt (tyvknorm(1)**2 + tyvknorm(2)**2 + tyvknorm(3)**2)

      do i=1,3
         tyvknorm(i) = tyvknorm(i)/norm
      enddo
      ierr = 0

      return 
      end


