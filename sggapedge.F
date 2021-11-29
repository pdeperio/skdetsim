***********************************************************************
*     -----------------
      SUBROUTINE SGGAPEDGE(NOPMT)
*     -----------------
*
*     (Purpose)
*     IN SK-II, Cherenkov photons which travel from the 
*     water gap(water between 
*     pmt and acrylic) to the defined edge of this gap(into the 
*     non-sensitive region) are treated by this routine.
*     Physically these photons continue on behind the PMT.  This routine
*     provides the possibility of these photons hitting the PMT.  This is an
*     improvement over earlier methods that merely stopped these photons.
*     This routine checks to see if the photon direction will intersect with
*     the PMT.  If it interects then the photon is moved to the intersection 
*     point, and sgpmt called.  Otherwise the photon is stopped.
*
*     NOTE::
*     This is not normal photon travel.  Scattering and absorption do 
*     not occur.  Photons cannot return from this region.  This is a rough
*     approximation in this small region.
*
*     SEE ALSO SGACREDGE, SGHOLEBST
*     (Input)
*       NOPMT::  PMT cable number
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       2002.08.21 ; First version by Chris Sterner
*
************************************************************************
      implicit none
#include "skparm.h"
#include "guskge.h"
#include "geant321/gctrak.inc"
#include "guxyzpm.h"
#include "geotnk.h"
#include "geant321/gctmed.inc"
#include "geant321/gconst.inc"
#include "skdetsim_test.h"

*
*  ZBLST in guskge.h is the size of the shift of black-sheet
*  from that of SK-I.
*  It is set according to the SKGEOMETRY
*  in sggeom_sk1, sggeom_sk2, sggeom_sk3
*

C     INPUT
      integer NOPMT
C     LOCAL
      real x1, y1, z1     !position of small pmt sphere in
                          !local coordiante system
      parameter(x1 = 0.0)
      parameter(y1 = 0.0)
      parameter(z1 = 0.0)
      real x0, y0, z0 !photon starting location
      real a, b, c !photon direction
      real length, xcos, xsin
      real p, q , r !photon intersection point
      real rx, ry, rz, rn, dist, xback
      real BLAMDA, travel
      real EFFNSG


*CKOV Photons stopped by istop=1 in gustep.F
*Picture of geometry:

C                            _____________________________
C                           /
C                          /         ____________________
C                         /         /    WATER GAP
C               ID       /  ACR    /         ___________
C                       /         /         /
C                      /         /         /   PMT
C                     /         /         /
C            ________/_________/_________/
C              bst|    acredge/  gapedge/
C                 |          /         /
C                 ----------/         /
C            flange bottom           /
C                                   /
C            - - - - - - - - - - - -
C               old bst level
C
C
C     Simple model above
C

*
*     Move coordinate system origin to pmt position
*                this reduces numerical errors
      x0=vect(1)-gxyzpm(1,NOPMT)
      y0=vect(2)-gxyzpm(2,NOPMT)
      z0=vect(3)-gxyzpm(3,NOPMT)
      a = vect(4)
      b = vect(5)
      c = vect(6)
      length = sqrt(a**2+b**2+c**2)
      if(length.lt.1.E-6 .and. abs(length-1.0) .gt. 1.E-6) then
         a = a/length
         b = b/length
         c = c/length
      endif

*     
*     Find intersection of photon direction and geant pmt sphere
*
c  (p,q,r) is
c  cross point of
c    x-x0   y-y0   z-z0
c    ---- = ---- = ----    a^2+b^2+c^2=1
c     a      b      c
c  and
c   (x-x1)^2 + (y-y1)^2 + (z-z1)^2 = r1^2
c
c  dist is distance from (x1,y1,z1) to the line
c
      rx=x1-x0
      ry=y1-y0
      rz=z1-z0
      rn=sqrt(rx**2+ry**2+rz**2)
      xcos=(a*rx+b*ry+c*rz)/rn
      if(xcos.gt.1.)xcos=1.
      if(xcos.lt.0.) then
         return !outgoing photon
      endif
      xsin=sqrt(1.-xcos**2)
      dist=rn*xsin
c
      p=a*rn*xcos+x0
      q=b*rn*xcos+y0
      r=c*rn*xcos+z0
c
      if(dist.gt.RPMT) then !no intersection
         return
      endif
      if(rn*xcos.gt.50.0) then !Photon travels too far, may be problem
         write(*,*) 'Error in sggapedge 1, please investigate'
         return
      endif
c
      xback=sqrt(RPMT**2-dist**2)-PREC !step back to just inside sphere
      p=p-xback*a
      q=q-xback*b
      r=r-xback*c
C     Update time of flight
      BLAMDA = 1237.8E-09 / GEKIN
      travel=sqrt((x0-p)**2+(y0-q)**2+(z0-r)**2)
      TOFG = TOFG +travel/CLIGHT*EFFNSG(BLAMDA)
c     p q r are (x,y,z) of intersection point if it exists


*
*     Move intersetion point back into Master Reference System
*
      p=p+gxyzpm(1,NOPMT)
      q=q+gxyzpm(2,NOPMT)
      r=r+gxyzpm(3,NOPMT)

*
*     Check the intersection point is within old black sheet
*               Do not allow photon to hit below old black sheet

      if(sqrt(p**2+q**2) .gt. RINTK ) then
         return !radius too large
      else if(r .gt. 0.5*HIINTK ) then
         return !above top
      else if(r .lt. -0.5*HIINTK ) then
         return !below bottom
      endif
      
    
*
*     Check the intersection point is behind new black sheet 
*
      if(sqrt(p**2+q**2) .lt. (RINTK-ZBLST) ) then
         if(r .lt. (0.5*HIINTK-ZBLST) .and. r .gt. (-0.5*HIINTK+ZBLST)) then
            return 
         endif
      endif

*
*     move photon to intersection point
*
      vect(1)=p
      vect(2)=q
      vect(3)=r

*
*     Photon should hit geant pmt
*     call sgpmt
*
      call sgpmt(NOPMT)


      RETURN
      END
