      subroutine skdcy_muplus(p_e,e_eplus)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  A routine to decay positive muons.  
C  Polarized and Unpolarized versions are available.
C
C  Last modified November 12, 2009 by J. Kameda, J. Albert 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Do we want to use the polarized version?
C#define USE_POLARIZED_VERSION

      implicit none
#include"weakparm.h"
#ifdef  USE_POLARIZED_VERSION
#include"geant321/gctrak.inc"
      integer polsn
      real costh, phi, polzn, sinth, e_max
      real costh_muon, sinth_muon, cosph_muon, sinph_muon
#endif
      real prob, entry, e_eplus
      real maxent, pabs
      real p_e(3)
      integer i

#ifdef  USE_POLARIZED_VERSION
C     If we're going to be polarized, what is our polarization?
C
C     POLZN = (N(-)-N(+))/(N(-)+N(+))
C
C     note that this is positive for mu+
C
      polzn = 0.25   ! default 0.25


C     first determine the muon spin orientation
C     if the muon has a pre-determined spin,
C       this is where to change the code
C     muplus tends to have spin anti-aligned with momentum
      call grndm(prob,1)
      polsn = 1  ! polarization sign (helicity)
      if ((prob*2-1).le.(polzn)) then
         polsn = -1
      endif

C     Kinematic maximum possible energy for decay e
      e_max = (mu_mass**2 + e_mass**2)/2.0/mu_mass

C     MC loop for E & theta
 200     continue
         call grndm(prob,1)
         e_eplus = prob * (mu_mass/2 -e_mass) + e_mass
         call grndm(prob,1)
         costh = prob * 2 - 1

C        Probability neglecting electron mass
C        From Cottingham & Greenwald Standard Model
C         maxent = 32/mu_mass**4 * e_eplus*e_eplus*((0.75*mu_mass - e_eplus)
C     &      -costh * (0.25*mu_mass - e_eplus))

C        More detailed decay rate
C        From *NEED CITATION*
         pabs = sqrt(e_eplus**2 - e_mass**2)
         maxent = 
     &   192.0/6.0*1.0/mu_mass**4 * pabs * e_eplus *
     &   (1.5*e_max - e_eplus - 0.5*e_mass**2/e_eplus
     &    +1.0*costh*(
     &    -0.5*e_max + e_eplus -0.5*e_mass**2/mu_mass
     &           )
     &   )

         call grndm(prob,1)
C        max prob slightly higher than max possible value of maxent
         entry = prob * 0.045
         if (entry.gt.maxent) goto 200

C     Now we get more angle stuff
      call grndm(prob, 1)
      phi = prob * 3.14159
      sinth = 1.0 - costh**2
C     Just to be safe......
      if (sinth .le. 0) then
         sinth = 0
      else
         sinth = sqrt(sinth)
      endif
C     And now let's sort out the momenta...
      p_e(1) = pabs * sinth * cos(phi) * 0.001
      p_e(2) = pabs * sinth * sin(phi) * 0.001
      p_e(3) = pabs * costh * 0.001
      e_eplus = e_eplus * 0.001

C     Next, we rotate to be aligned with the muon spin
C     Muon spin 

       costh_muon  = polsn*vect(6)
       sinth_muon  = sqrt(1-costh_muon**2)    ! 0<theta<pi
       cosph_muon  = polsn*vect(4)/sinth_muon
       sinph_muon  = polsn*vect(5)/sinth_muon

       call gdrot(p_e, costh_muon, sinth_muon, cosph_muon, sinph_muon)

C      Warning
       write(*,*) '::::::::::::::POLARIZED DECAY::::::::::::::::::::'
       write(*,*) 'decay e z-component is: ', p_e(3)
   
#endif
#ifndef USE_POLARIZED_VERSION      
 123     continue

         Call grndm(prob,1)
         entry = prob*2
         Call grndm(prob,1)
         e_eplus = prob

         maxent = e_eplus**2*(6-e_eplus*4)
         if(entry.le.maxent)then
            pabs = mu_mass*e_eplus/2.
         else
            pabs = 0.
         end if

         if (pabs .eq. 0.0) goto 123
         e_eplus = sqrt(pabs**2+e_mass**2)

         CALL RNDIR(p_e)
         p_e(1) = p_e(1)*pabs*0.001
         p_e(2) = p_e(2)*pabs*0.001
         p_e(3) = p_e(3)*pabs*0.001

#endif

         return
         end

