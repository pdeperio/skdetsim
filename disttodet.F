      subroutine disttodet( pos, dist )
*     ********************************************************
*     Calculates distance out of the rock
*     input: pos(3)
*     output: dist
*     
*     2012.5.14: T. Wongjirad. Accounted for distance to the dome.
*
*     ********************************************************
      implicit none

#include "geotnk.h"
#include "sggeom_dome.h"

      REAL POS(3)      
      REAL dist
      REAL r
      REAL zabs
      REAL dz
      REAL dr

c     Get Radius and ABS(Z)
      dist = -1
      r = sqrt(pos(1)*pos(1)+pos(2)*pos(2))
      zabs = sqrt( pos(3)*pos(3) )

      if (r.le.RTKTK .and. zabs.le.ZPTKTK) then
c     if inside the detector, return 0
         dist = 0
      else if (r.gt.RTKTK .and. zabs.le.ZPTKTK) then
c     if outside r, but within end cap planes
         dist = r - RTKTK
      else if (pos(3).le.ZMTKTK.and.r.le.RTKTK) then
c     below bottom cap plane, but inside cylinder
         dist = zabs - ZPTKTK
      else if (pos(3).le.ZMTKTK.and.r.ge.RTKTK) then
c     bottom corner
         dr = r - RTKTK
         dz = zabs - ZPTKTK
         dist = sqrt(dr*dr+dz*dz)
      else if (pos(3).ge.ZPTKTK) then
c     above the water tank top cap plane
         if ( pos(3).le.ZPTKTK+400.and.r.ge.RTKTK ) then
c     below the spherical part of dome, outside of tank radius
            dist = r-RTKTK
         else if (pos(3).le.ZPTKTK+400.and.r.le.RTKTK) then
c     inside the dome
            dist = 0.0
         else 
c     in the region where spherical part of dome is closest surface
c     calculate distance to dome's center
            dist = sqrt( r*r 
     &           +(pos(3)-(dome_cyl_center(3)+dome_sphe_center(3)))**2)
            if (dist.ge.dome_sphe_radius) then
c     outside of dome in rock
               dist = dist - dome_sphe_radius
            else
c     inside dome in air
               dist = 0.
            endif
         endif
       ! end of if above top cap plane
      else
c     should never get here
         write(*,*) "ERROR disttotdet.F: ",
     &        "entered area in rock not accounted for."
         write(*,*), 'r=',r,'z=',pos(3)
         stop
      endif

      if (dist.lt.0) then
         write(*,*) "ERROR DIST OUTSIDE OF DET IS ",
     &        "NEGATIVE. SHOULD NEVER GET HERE"
         dist = 0
         stop
      endif

      return

      end
