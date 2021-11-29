************************************************************************
*     ----------------------------------------------
      SUBROUTINE SGHOLEBST(NEARPMT, INGAP, DISTANCE)
*     ----------------------------------------------
*
*     (Purpose)
*      The purpose of routine is to check whether the photon is hitting the
*      gap between the acrylic case and black sheet. This routine should be
*      called only when GEANT tracking says photon is hitting the black sheet.
*
*      SKAM version should be more than or equal to "03b".
*
*
*     (Input)
*       none
*
*     (Output)
*       NEARPMT : Cable number of the nearest PMT
*       INGAP   : =1 : Photon hit position is in the gap between black
*                      sheet and acrylic case
*                 =0 : Photon is actually hitting the black sheet
*       DISTANCE : Distance between the hit position and the axis of
*                  PMT
*
*     (Creation Date and Author)
*       2003.03.09 ; First version by M. Nakahata
************************************************************************
      implicit none
#include "skparm.h"
#include "skhead.h"
#include "skdetsim_parm.h"
#include "skcnst.h"
#include "geant321/gctrak.inc"
#include "geotnk.h"
#include "guskge.h"
#include "geopmt.h"
#include "guxyzpm.h"
#include "skmisch.h"

*
*  ZBLST in guskge.h is the size of the shift of black-sheet
*  from that of SK-I.
*  It is set according to the SKGEOMETRY
*  in sggeom_sk1, sggeom_sk2, sggeom_sk3
*

C -- Radius of the PMT hole on blacksheet 

      Real HOLE_IN_BLACKSHEET
      PARAMETER ( HOLE_IN_BLACKSHEET = 28.5)

      Integer INGAP, NEARPMT
      Real DISTANCE
      Real tiny
      PARAMETER ( tiny = 0.01)
      Integer ix, iy, iphi, iz, ixx, iyy
      Real phi, dphi, offset, x, y, stepz

***   This subroutine is only for SK_II and SK_III

      IF(SK_GEOMETRY .eq. SK_I) go to 9999

** where is the black sheet top, bottom or barrel ?

      IF (VECT(3).gt.(ZPINTK-ZBLST-tiny)) THEN     ! hitting top
         IF (VECT(1).ge.0.) then
            ix = IFIX(VECT(1)/STEPX) + 1
         ELSE
            ix = IFIX(VECT(1)/STEPX)
         END IF
         IF (VECT(2).ge.0.) then
            iy = IFIX(VECT(2)/STEPY) + 1
         ELSE
            iy = IFIX(VECT(2)/STEPY)
         END IF
         ixx = ix + MXTXPM/2
         iyy = iy + MXTYPM/2
         IF(ixx.lt.1 .or. ixx.gt.MXTXPM .or.
     &      iyy.lt.1 .or. iyy.gt.MXTYPM)  go to 9999
         NEARPMT = ISEQT(ixx, iyy)
         IF (NEARPMT.lt.1 .or. NEARPMT.gt.MAXPM) go to 9999
         IF (IMIS(NEARPMT) .ne.0) go to 9999
         DISTANCE = sqrt ( (gxyzpm(1,NEARPMT)-VECT(1))**2 +
     &                     (gxyzpm(2,NEARPMT)-VECT(2))**2)
         IF (DISTANCE.le.HOLE_IN_BLACKSHEET) then
            INGAP = 1
         else
            INGAP = 0
         END IF
      ELSE IF (VECT(3).lt.(ZMINTK+ZBLST+tiny)) THEN  ! hitting bottom
         IF (VECT(1).ge.0.) then
            ix = IFIX(VECT(1)/STEPX) + 1
         ELSE
            ix = IFIX(VECT(1)/STEPX)
         END IF
         IF (VECT(2).ge.0.) then
            iy = IFIX(VECT(2)/STEPY) + 1
         ELSE
            iy = IFIX(VECT(2)/STEPY)
         END IF
         ixx = ix + MXTXPM/2
         iyy = iy + MXTYPM/2
         IF(ixx.lt.1 .or. ixx.gt.MXTXPM .or.
     &      iyy.lt.1 .or. iyy.gt.MXTYPM)  go to 9999
         NEARPMT = ISEQB(ixx, iyy)
         IF (NEARPMT.lt.1 .or. NEARPMT.gt.MAXPM) go to 9999
         IF (IMIS(NEARPMT).ne.0) go to 9999
         DISTANCE = sqrt ( (gxyzpm(1,NEARPMT)-VECT(1))**2 +
     &              (gxyzpm(2,NEARPMT)-VECT(2))**2)
         IF (DISTANCE.le.HOLE_IN_BLACKSHEET) then
            INGAP = 1
         else
            INGAP = 0
         END IF
      ELSE                     !  barrel case
         phi = ATAN2(VECT(2), VECT(1))
         IF (phi.lt.0.0) phi = phi + 2 * PI
         dphi = 2*PI/float(MXWTPM)
         iphi = ifix(phi/dphi+0.5)
         If (iphi.ge.MXWTPM) iphi = 0
         offset = (float(MXWZPM)+1.0)/2.0
         STEPZ=HIWAL/FLOAT(MXWZPM)
         iz = IFIX(VECT(3)/STEPZ+offset+0.5)
         NEARPMT = iphi * MXWZPM + iz
         IF (NEARPMT.lt.1 .or. NEARPMT.gt.(MXWZPM*MXWTPM))go to 9999
         IF (IMIS(NEARPMT).ne.0) go to 9999
         x = gxyzpm(1,NEARPMT) * (RINTK-ZBLST)/RINTK
         y = gxyzpm(2,NEARPMT) * (RINTK-ZBLST)/RINTK
         DISTANCE = sqrt ( (x-VECT(1))**2 + (y-VECT(2))**2 +
     &              (gxyzpm(3,NEARPMT)-VECT(3))**2)
         IF (DISTANCE.le.HOLE_IN_BLACKSHEET) then
            INGAP = 1
         else
            INGAP = 0
         END IF
      END IF

      RETURN

C -- no relevant answer (either SK-I or problematic photon position)
 9999 CONTINUE
      INGAP = 0
      NEARPMT = 0
      RETURN
      END
