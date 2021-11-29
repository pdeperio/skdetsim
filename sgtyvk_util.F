************************************************************************
* This is just a file containing the subroutines needed                *   
* by sgtyvk.F and sgytvk_odsegm.F                                      *
*                                                                      *
*         Fanny Dufour, July 2005                                      *  
*     
*  2008.01 : changes to incorporate new reflectivity model             *
*            see sgtyvk.F,  R.Wendell (rvw)                            *
*  2020.05 : Add clamp to incident angle calculation to avoid possible *
*            NaN                                                       *
*                                                                      *
************************************************************************



      subroutine sglambert(v)
      
* Return a vector according to Lambert's Law reflection,
* assuming normal is (0,0,1)
* th returned is 0 in xy plane
       
      implicit none
      
      real twopi
      parameter (twopi= 6.2831853)
      real v(3)
      real ph,th
      real RDUMMY
      real cosph,sinph,costh,sinth
      
* Uniform phi
      
      CALL GRNDM(RDUMMY,1)
      ph=RDUMMY*twopi
      
*     theta according to cos(theta) distribution
      
      CALL GRNDM(RDUMMY,1)
      th = asin(2.*RDUMMY-1)
      
*     -pi/2 to pi/2
      
      cosph=cos(ph)
      sinph=sin(ph)
      costh=cos(th)
      sinth=sin(th)
      v(1)=sinth*cosph
      v(2)=sinth*sinph
      v(3)=costh
      
      return
      end
      
******************************************************
      
      subroutine sgisotropic(v)
      
      implicit none
      
      real ph,th
      real v(3)
      real sinth,cosph,sinph
      real rdummy
      real twopi
      parameter (twopi= 6.2831853)
      
*     Return vector from uniform isotropic distribution over 2pi sr
*     positive hemisphere for normal (0,0,1)
                                                                                
      CALL GRNDM(RDUMMY,1)
      ph=RDUMMY*twopi
      
      CALL GRNDM(RDUMMY,1)
      th = acos(RDUMMY)
       
      cosph=cos(ph)
      sinph=sin(ph)
      sinth=sin(th)
      
      v(1)=sinth*cosph
      v(2)=sinth*sinph
      v(3)=RDUMMY
      
      
      return
      end

************************************************************
      
      subroutine sggauss_refl(inc_dir, v)
      
* Return a vector with Gaussian theta, sigma = 15 deg
* assuming normal is (0,0,1)
                                                                                
      implicit none
      
      real pi
      parameter (pi=3.14159)
      real twopi
      parameter (twopi= 6.2831853)
      real inc_dir(3)
      real v(3)
      real theta_in, phi_in     ! Incident photon angles
      real sig_theta
!     data sig_theta /0.262/  ! 15.0 deg, the old model
      real ph,th
      real RDUMMY
      real cosph,sinph,costh,sinth
      real grndm
      external grndm
      
      ! array of gaussian width's measured by A. Chavarria, for incident
      ! angles in degrees
      real sig_array(9)
      data sig_array/31.3, 35.4, 36.0, 35.8, 29.8, 22.9,22.2, 26.6,34.5/
      integer angbin

      call sgget_angles(inc_dir,theta_in,phi_in)
      
* Probably should have this isotropic in phi around reflected vector?
                                                                        
      ph = phi_in
*  Select proper gaussian width from measured values (rvw)
      angbin = nint( ( 3.14159 - theta_in ) * 180./3.14159/10.) + 1
      if (angbin.le.0 ) angbin=1
      if (angbin.gt.9 ) angbin=9                  ! A.Chavarria didn't measure incidence > 80 degrees
      sig_theta = sig_array(angbin)*3.14159/180.  ! convert to radians                
             
                                    
*     Theta is Gaussian around specular reflected direction
      th = 10000.0
      do while( abs(th) .gt. 3.14159/2.0 )  
        call RNORMX(RDUMMY,1,GRNDM)
        th = pi -theta_in + rdummy*sig_theta
      enddo      
                                                                        
      cosph=cos(ph)
      sinph=sin(ph)
      costh=cos(th)
      sinth=sin(th)
      v(1)=sinth*cosph
      v(2)=sinth*sinph
      v(3)=costh
     
      return
      end
      
      
*************************************************************

      subroutine sgget_angles(dir,th,ph)
      
*     Return theta,phi for given input vector.  Convention
*     is that th=0 is for (0,0,1) pointing straight up
      
      real th,ph
      real dir(3)
      real pi
      parameter(pi=3.141596)
      
      if(dir(3).lt.-1.) dir(3) = -1.
      if(dir(3).gt.1.) dir(3) = 1.
      th = acos(dir(3))
      if (dir(1).gt.0) then
         if (dir(2).gt.0) then
            ph = atan(dir(2)/dir(1))
         else
            ph = atan(dir(2)/dir(1))+ 2.*pi
         endif
      else if (dir(1).lt.0) then
         ph = atan(dir(2)/dir(1))+pi
      else
         if (dir(1).eq.0) then
            ph = 0.
         else if (dir(2).gt.0) then
            ph = pi/2.
         else
            ph = pi*1.5
         endif
      endif
      
      return
      end
      
      
      
