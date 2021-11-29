
***********************************************************************
*     -----------------
      SUBROUTINE SGACREDGE(NOPMT,LOCATION)
*     -----------------
*
*     (Purpose)
*     Cherenkov photons which travel from the acrylic medium
*     to the defined edge of this medium(into the non-sensitive region) are
*     treated by this routine.
*     Physically these photons enter the flange region of the acrylic
*     cover which is behind the black sheet, but this
*     routine only approximates this region.
*     Similar to sggapedge, this routine finds the intersection of the
*     photon direction and the pmt surface.  Distance traveled in the acrylic
*     is also determined to account for absorption in the acrylic.  The 
*     acrylic volume is approximated by the continuation of the sphere.  
*     Thus the flange of the acrylic does not extend out.  Reflection
*     and refraction of the water acrylic interface is not treated.  
*     Scattering and absorption in the water gap are not treated.  Reflection 
*     from the pmt is not treated.  The bottom of the acrylic flange (at
*     ZACREDGE) stops all incident photons.
*
*     NOTE::
*     This is not normal GEANT photon travel.  We could not easily and quickly
*     model the region using GEANT volumes.  This routine corrects problems of
*     collecting photons encountered by not putting this region in sggeom_sk2.
*     Photons cannot return from this 
*     region.  This routine is only to help collect some extra photons and
*     get a more correct pmt angular acceptance distribution.
*
*     SEE ALSO SGGAPEDGE, SGHOLEBST
*     (Input)
*     NOPMT              Cable number of PMT
*     LOCATION           0 (=photon from ID water)
*                        1 (=photon from Acrylic Shell)
*     (Output)
*       none
*
*     (Creation Date and Author)
*       2002.08.21 ; First version by Chris Sterner
*       2006.04.29 ; Add BST collar in SK-III PMTs by Jun Kameda
*
************************************************************************
      implicit none

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
C         __ ________/_________/_________/
C       bst   | hole/  acredge/  gapedge/
C                  /         /         /
C                 ----------/         /
C            flange bottom           /
C                                   /
C            - - - - - - - - - - - -
C               old bst level
C
C
C     Simple model above
C     hole is determined by sgholebst, treated by sgacredge
#include "skparm.h"
#include "guxyzpm.h"
#include "geant321/gctrak.inc"
#include "guskge.h"
#include "geotnk.h"
#include "skdetsim_test.h"
#include "geant321/gctmed.inc"
#include "geant321/gconst.inc"
#include "skpmt.h"

*
*  ZBLST in guskge.h is the size of the shift of black-sheet
*  from that of SK-I.
*  It is set according to the SKGEOMETRY
*  in sggeom_sk1, sggeom_sk2, sggeom_sk3
*

C     INPUT
      integer NOPMT, LOCATION

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
      real BLAMDA, travel, zpass
      real EFFNSG, zintac, stepac, suacr, distacr, rdummy
      real acrnsg, acrsg, dista
      real ap,aq,ar,pq,pp,pr,aback,xpass,ypass
      real outap,outaq,outar
      real xr,yr,zr,prod_vn
      real xp,yp,zp,bp,cp,cos_vr,sin_vr
      real dist_p,alpha,prod_cn
      real xc,yc,zc,prod_rn
      real a_pro, b_pro, c_pro
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
      if(xcos.lt.0) then
         return !outgoing photon
      endif
      xsin=sqrt(1.-xcos**2)
      dist=rn*xsin
c
      p=a*rn*xcos+x0
      q=b*rn*xcos+y0
      r=c*rn*xcos+z0
c
      if(dist.gt.RPMT) then !no intersection with pmt
         return
      endif
      if(rn*xcos.gt.100.0) then !Photon travels too far, may be problem
         write(*,*) 'Error in sggapedge 1, please investigate'
         return
      endif
c

      xback=sqrt(RPMT**2-dist**2)-PREC !step back to just inside sphere
      pp=p-xback*a
      pq=q-xback*b
      pr=r-xback*c
C
C     TEST absorption by acrylic
C
c     determine distance traveled in acylic
c     find crossing point of inner surface
      rx=x1-(RACR - ZACCY)*gdxyzpm(1,NOPMT)-x0
      ry=y1-(RACR - ZACCY)*gdxyzpm(2,NOPMT)-y0
      rz=z1-(RACR - ZACCY)*gdxyzpm(3,NOPMT)-z0
      rn=sqrt(rx**2+ry**2+rz**2)
      xcos=(a*rx+b*ry+c*rz)/rn
      if(xcos.gt.1.)xcos=1.
      xsin=sqrt(1.-xcos**2)
      dista=rn*xsin    
c     nearest point is pqr
      p=a*rn*xcos+x0
      q=b*rn*xcos+y0
      r=c*rn*xcos+z0
c     cross out of acr at ap,aq,ar
      aback=sqrt((RACR-WIDACR/2.0)**2-dista**2)
      ap=p-aback*a
      aq=q-aback*b
      ar=r-aback*c
c
C
c     distance is then
C
      IF(LOCATION.eq.1) then !photon from acrylic 
         distacr=sqrt((x0-ap)**2+(y0-aq)**2+(z0-ar)**2)
      else if(LOCATION.eq.0) then !photon from ID Water
C        find crossing point of outer surface
         aback=sqrt((RACR+WIDACR/2.0)**2-dista**2)
         outap=p-aback*a
         outaq=q-aback*b
         outar=r-aback*c
         distacr=sqrt((outap-ap)**2+(outaq-aq)**2+(outar-ar)**2)
      else
         write(*,*) 'ERROR sgacredge LOCATION parameter'
         stop
      endif

      CALL GRNDM(RDUMMY,1)
      ZINTAC = -LOG(RDUMMY)
      BLAMDA = 1237.8E-09 / GEKIN
      STEPAC = ACRSG(BLAMDA)
ccyt cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (stepac.lt.1.0e-4) stepac = 1.0e-4
ccyt cccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUACR  = ZINTAC * 100./STEPAC
      IF(SUACR.LT.distacr) THEN
         return                 !absorb in acrylic
      endif

c
c    Blocked by the blacksheet collar
c    SK-III PMTs only have the BST collar inside the FRP case
c

      if ( NSKPMT(nopmt) .eq. 3) then   
*   relative position of the photon from PMT position
      xr = x0 -x1
      yr = y0 -y1
      zr = z0 -z1

*  inner product of relative position vector and normal vector
      prod_rn = xr*gdxyzpm(1,NOPMT)+yr*gdxyzpm(2,NOPMT)+zr*gdxyzpm(3,NOPMT)

*  inner product of direction vector and normal vector
      prod_vn = a*gdxyzpm(1,NOPMT)+b*gdxyzpm(2,NOPMT)+c*gdxyzpm(3,NOPMT)

*  projected vector of the relative position to the BST plane

      xp = xr - prod_rn*gdxyzpm(1,NOPMT)
      yp = yr - prod_rn*gdxyzpm(2,NOPMT)
      zp = zr - prod_rn*gdxyzpm(3,NOPMT)

*  projected vector of the direction to the BST plane

      a_pro = a - prod_vn*gdxyzpm(1,NOPMT)
      b_pro = b - prod_vn*gdxyzpm(2,NOPMT)
      c_pro = c - prod_vn*gdxyzpm(3,NOPMT)

*  Check whether the photon hit the BST collar
      cos_vr = 
     & -(xp*a_pro+yp*b_pro+zp*c_pro)/sqrt(xp**2+yp**2+zp**2)/sqrt(a_pro**2+b_pro**2+c_pro**2)

      if (cos_vr .lt.0) return  ! outgoind photon
      if (cos_vr .gt. 1) cos_vr = 1
      sin_vr = sqrt(1-cos_vr**2)
      dist_p = sqrt(xp**2+yp**2+zp**2)*sin_vr

       if (dist_p .lt. RFRP) then 
         alpha = sqrt(xp**2+yp**2+zp**2)*cos_vr - sqrt(RFRP**2-dist_p**2)
         alpha = alpha/sqrt(a_pro**2+b_pro**2+c_pro**2)
          xc = x0 + a * alpha
          yc = y0 + b * alpha
          zc = z0 + c * alpha
         prod_cn = (xc-x1)*gdxyzpm(1,NOPMT)
     &   +(yc-y1)*gdxyzpm(2,NOPMT)+(zc-z1)*gdxyzpm(3,NOPMT)
         if (prod_cn .le. ZCOLLAR) then ! hit BST collar
            return        
         endif

       endif

      endif

C     Update time of flight
      travel=sqrt((x0-pp)**2+(y0-pq)**2+(z0-pr)**2)
      TOFG = TOFG +(travel-distacr)/CLIGHT*EFFNSG(BLAMDA)
     &     +distacr/CLIGHT*ACRNSG(BLAMDA)

c     pp pq pr are (x,y,z) of intersection point if it exists
*
*     Move intersetion point back into Master Reference System
*
      pp=pp+gxyzpm(1,NOPMT)
      pq=pq+gxyzpm(2,NOPMT)
      pr=pr+gxyzpm(3,NOPMT)

C
C     test if photon hits bottom of acrylic flange (about 4cm tall)
C
      xpass=ap+gxyzpm(1,NOPMT)
      ypass=aq+gxyzpm(2,NOPMT)
      zpass=ar+gxyzpm(3,NOPMT)
      if(sqrt(xpass**2+ypass**2) .gt. RINTK-ZACREDGE) then
c         write(41,*) ap,aq,ar,nopmt,nskpmt(nopmt)
         return ! hit edge of acrylic on wall of tank
      else if(abs(zpass) .gt. 0.5*HIINTK-ZACREDGE ) then
c         write(41,*) ap,aq,ar,nopmt,nskpmt(nopmt)
         return ! hit edge at top or bottom of tank
      endif
*
*     Check the intersection point is within old black sheet
*               Do not allow photon to hit below old black sheet
*               Do not allow photon to hit backside of pmt 

      if(sqrt(pp**2+pq**2) .gt. RINTK ) then
         return !radius too large
      else if(pr .gt. 0.5*HIINTK ) then
         return !above top
      else if(pr .lt. -0.5*HIINTK ) then
         return !below bottom
      endif
      
*
*     Check the intersection point is behind new black sheet 
*
      if(sqrt(pp**2+pq**2) .lt. (RINTK-ZBLST) ) then
         if(pr .lt. (0.5*HIINTK-ZBLST) .and. pr .gt. (-0.5*HIINTK+ZBLST)) then
            return 
         endif
      endif

*
*     move photon to intersection point
*
      vect(1)=pp
      vect(2)=pq
      vect(3)=pr

*
*     Photon should hit geant pmt
*     call sgpmt
*
      call sgpmt(NOPMT)
 
      RETURN
      END
