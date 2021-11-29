************************************************************************
*     ----------------------
      FUNCTION RFBSSG(COSTH,ipol)
*     ----------------------
*
*     (Purpose)
*        Get the ratio of reflection on black sheet
*
*     (Input)
*        COSTH  : input angle
*
*     (Output)
*        RFBSSG : ratio of reflection
*
*     (Creation Date and Author)
*       1996. 6.27 ; First version by Y.Koshio
*         ref. (1) reflection value from K.Inoue
*                  mesured by ???
*       1998. 7.28 : Ensure function return value is always defined on return
*                    D.Casper
*       2001. 7.24 : Modified tuned by Shibata
*                            
************************************************************************

#include "skhead.h"
#include "skdetsim_parm.h"
#include "geant321/gctrak.inc"
c M Fechner -- tuning parameters inside this file
#include "reflections.h"

      integer ipol
      REAL COSTH,RFBSSG,A
      REAL NBS,B

      real sgtuned_bst, factor
      external sgtuned_bst

      IF ( (REFLECTION_TUNE .eq. LE_TUNE) .or.
     &     (REFLECTION_TUNE .eq. UNI_TUNE) ) THEN ! LE_,UNI_TUNE
         IF(SK_GEOMETRY.eq.SK_I) then 
            NBS=1.6/1.34
            if(costh.lt.0.) return
            A=ACOS(COSTH)
C     added for avoiding "divided by zero".
            if(costh.gt.0.99)a=acos(0.99)

            b=asin(sin(a)/nbs)
            if(ipol.eq.0)rfbssg=0.5*((sin(a-b)/sin(a+b))**2
     &           +(tan(a-b)/tan(a+b))**2)      
            if(ipol.eq.1)rfbssg=(tan(a-b)/tan(a+b))**2      
            if(ipol.eq.2)rfbssg=(sin(a-b)/sin(a+b))**2      
            blamda=1237.8E-09 / GEKIN
            RETURN

         ELSE IF(SK_GEOMETRY.eq.SK_II) then
            NBS=2.01/1.34       ! lowe (value by Parker 05/04/15)
         ELSE IF(SK_GEOMETRY.eq.SK_III) then
            NBS=2.01/1.34       ! lowe (value by Parker 05/04/15)
         ELSE IF(SK_GEOMETRY.eq.SK_IV) then
            NBS=2.01/1.34       ! lowe (value by Parker 05/04/15)
         ELSE IF(SK_GEOMETRY.eq.SK_V ) then ! just copy from SK-IV J.Kameda 2019/10/29
            NBS=2.01/1.34       
         ELSE
            write(*,*) 'RFBSSG: Not Supported SK_GEOMETRY',SK_GEOMETRY
            stop
         ENDIF
      ELSE IF(REFLECTION_TUNE .EQ. AP_TUNE) THEN
         NBS=1.6/1.33           ! atmpd
      else
         write(*,*) 'RFBSSG: Not supported REFLECTION_TUNE.', reflection_tune
         stop
      END IF

      RFBSSG = 0.
      IF (abs(costh) .GT. 1.) RETURN
CCC      A=ACOS(COSTH)	! bug, this caculate in radian
      A=ACOS(COSTH)*180./3.141593

      IF(A.LT.0. .OR. A.GT.90.)RETURN	! compare in degree

      A=ACOS(COSTH)

      b=asin(sin(a)/nbs)
cc      if(ipol.eq.0)rfbssg=0.5*((sin(a-b)/sin(a+b))**2
cc     &     +(tan(a-b)/tan(a+b))**2)      
cc      if(ipol.eq.1)rfbssg=(tan(a-b)/tan(a+b))**2      
cc      if(ipol.eq.2)rfbssg=(sin(a-b)/sin(a+b))**2      

      if (ipol .eq. 1) then ! P-wave
         rfbssg = ((nbs*cos(a) - cos(b)) / (nbs*cos(a) + cos(b)))**2
      elseif (ipol .eq. 2) then ! S-wave
         rfbssg = ((cos(a) - nbs*cos(b)) / (cos(a) + nbs*cos(b)))**2
      else
         rfbssg = 0.5*((nbs*cos(a) - cos(b)) / (nbs*cos(a) + cos(b)))**2
     &           +0.5*((cos(a) - nbs*cos(b)) / (cos(a) + nbs*cos(b)))**2
      endif

c
      blamda=1237.8E-09 / GEKIN
c

      if ((reflection_tune .ne. ap_tune) .and.   ! correction
     &   (sk_geometry .ge. SK_III)) then
         factor = sgtuned_bst(costh,blamda)
         rfbssg = rfbssg*factor
      endif

      IF (REFLECTION_TUNE .EQ. AP_TUNE) THEN
         xfac=tunedat(blamda,tuned_bs,at,0,ntuning)
         rfbssg=rfbssg*(1.-xfac)+xfac ! atmpd
      END IF
      RETURN
      END

******************************************************************
******************************************************************
      real function sgtuned_bst(costh,lambda)
*
*     A correction factor to Black sheet reflectivity.
*     Depending on incident angle and wave length.
*    
*     This correction factors are obtained from
*     direct measurement of the reflected light using 
*     laser injector data.
*     
*     BSTREF = 0.5 is assumed. Please see dsinip.F .
*     Assuming the correction is 1.0 at 0 degree, because 
*     laser injector data for water transparency is indicating.
*     Assuming the correction is 0.5 at 90 degree.
*
*     2009/06/24 Jun Kameda
*

      implicit none

      integer i_wave,i_theta, i
      real costh, lambda
      real theta

      real wave(3), ang(5), weight(5,3)
      real r_wave,r_theta

      data wave/337.0, 400.0, 420.0/
      data ang/0.0, 30.0, 45.0, 60.0, 90.0/
      data weight/

c     &  2.00000, 1.76324, 1.48192, 1.18563, 1.00000, 
c     &  2.00000, 1.31726, 1.17261, 0.94415, 1.00000,
c     &  2.00000, 1.34504, 1.20195, 0.92380, 1.00000/

c     &  2.00000, 1.67508, 1.40782, 1.12634, 1.00000, 
c     &  2.00000, 1.31726, 1.17261, 0.94415, 1.00000,
c     &  2.00000, 1.34504, 1.20195, 0.92380, 1.00000/

c     &  2.00000, 1.55019, 1.32194, 1.12634, 1.00000, 
c     &  2.00000, 1.16981, 1.14212, 0.92697, 1.00000,
c     &  2.00000, 1.18382, 1.16093, 0.88408, 1.00000/

c     &  2.00000, 1.56569, 1.33516, 1.13760, 1.00000, 
c     &  2.00000, 1.16981, 1.14212, 0.92697, 1.00000,
c     &  2.00000, 1.18382, 1.17254, 0.88408, 1.00000/

* Tarek Akiri, July 2012
*     Old best tuning values affected by Rayleigh scattering bug
**     &  2.00000, 1.56569, 1.33516, 1.14897, 1.00000, 
**     &  2.00000, 1.16981, 1.14212, 0.92697, 1.00000,
**     &  2.00000, 1.19565, 1.18423, 0.89292, 1.00000/

*     New tuning values with bug corrected
     &     2.00000, 1.60904, 1.36482, 1.20353,  1.00000, 
     &     2.00000, 1.24411, 1.13971, 0.970104, 1.00000,
     &     2.00000, 1.24317, 1.14026, 0.951023, 1.00000/

      sgtuned_bst = 1.0

      theta = acos(costh)*180.0/3.141593  ! same as rfbssg.F

      if ((theta .lt. 0) .or. (theta .gt. 90.0)) return
      if (lambda .lt. 0) return

cc scan wave length

      if (lambda .le. wave(1)) then
         r_wave = 0.0
         i_wave = 1
         goto 100
      endif

      if (lambda .ge. wave(3)) then
         r_wave = 1.0
         i_wave = 2
         goto 100
      endif
         
       do i = 2, 3
          if (lambda .le. wave(i)) then
             r_wave = (lambda - wave(i-1))/(wave(i)-wave(i-1))
             i_wave = i -1
             goto 100
          endif
       enddo
       
       write(*,*) 'sgtuned_bst:: Oops'

 100   continue

cc angle

       do i = 2, 5
          if (theta .le. ang(i)) then
             r_theta = (theta - ang(i-1))/(ang(i)-ang(i-1))
             i_theta = i -1
             goto 200
          endif
       enddo

 200   continue

      sgtuned_bst =
     &     weight(i_theta  ,i_wave  )*(1.0 - r_theta)*(1.0 - r_wave)
     &   + weight(i_theta+1,i_wave  )*r_theta*(1.0 - r_wave)
     &   + weight(i_theta  ,i_wave+1)*(1-r_theta)*r_wave
     &   + weight(i_theta+1,i_wave+1)*r_theta * r_wave

      return
      end
