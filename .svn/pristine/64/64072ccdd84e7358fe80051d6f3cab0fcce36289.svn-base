***********************************************************************
*     -----------------
      SUBROUTINE SGACRNOR(POS,NORMAL,ACRTRAN)
*     -----------------
*
*     (Purpose)
*     This routine determines the vector normal to the surface (inner 
*     or outer) of the
*     acrylic shperes.  Because of limitations of the Geant standard fuctions
*     GGPERP, GMTOD, and GDTOM this routine does not use those functions.
*     NOTE!! This routine is for the specific case/geometry of photon tracking
*     near the surface of the acrylic covers as defined in the SK-II sggeom.F
*     as of 09/2002.  Assumption, acrylic covers are spheres, thus simplifying
*     calculation.
*     
*     This is not a general routine.
*
*
*     (Input)
*       POS(3)::     position of the photon near the water/acrylic interface
*                     in the Master Reference System
*       ACRTRAN(3)::    Position of the acrylic sphere's 
*                        reference system in MARS
*                       Filled in GUSTEP.F when photon is in acrylic
*                         if IENTER=0 then then this should reference the 
*                         sphere being exited not current volume
*
*     (Output)
*       NORMAL(3)::  direction normal to the water/acrylic interface
*                     in the Master Reference System.  May point in or out.
*
*     (Creation Date and Author)
*       2002.09.20 ; First version by Chris Sterner
*
************************************************************************

      IMPLICIT NONE
      
*     INPUT
      REAL POS(3), ACRTRAN(3)
*     OUTPUT
      REAL NORMAL(3)
*     LOCAL
      REAL DXL(3), UL(3), R

*     Transform position, POS, into the acrylic local reference system, DXL
*       (based on GMTOD)
      DXL(1) = POS(1) - ACRTRAN(1) 
      DXL(2) = POS(2) - ACRTRAN(2)
      DXL(3) = POS(3) - ACRTRAN(3)

*     Determine the radial direction, UL, pointing towards POS/DXL
      R = SQRT(DXL(1)**2+DXL(2)**2+DXL(3)**2)
      IF(R .NE. 0.0) then
         NORMAL(1) = DXL(1) / R
         NORMAL(2) = DXL(2) / R
         NORMAL(3) = DXL(3) / R
      ELSE
         WRITE(6,*) 'SGACRNOR:: ERROR in normal calculation'
      ENDIF

      RETURN
      END


