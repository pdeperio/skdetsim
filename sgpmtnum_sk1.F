************************************************************************
*     -------------------------------------------------
      SUBROUTINE SGPMTNUM_SK1(NAMES,n5,n6,n7,n8,NOPMT,iflag)
*     -------------------------------------------------
*     
*     (Purpose)
*     *     This routine is used to determine which pmt site a GEANT 
*     volume corresponds to.
*     This routine determines the pmt cable number from the GEANT volume's 
*     position in the geometrical heirarchy, as defined by NAMES 
*     and NUMBER(*).  As the Acrylic cover, Water gap, and inner pmt at a 
*     location are all daughters of the same mother, this routine returns the
*     same cable number for all three volumes at a pmt site.  
*     iflag notes if this volume corresponds to an irregularly 
*     removed pmt, i.e. if this pmt is defined in GEANT but is not in SK-I.
*     
*     (Input)
*     NAMES ; name of geometry (5th array)
*     n5    ; number(5) in gcvol
*     n6    ; number(6) in gcvol
*     n7    ; number(7) in gcvol
*     n8    ; number(8) in gcvol
*     
*     (Output)
*     NOPMT ; cable number
*     iflag ; 1:pmt exist, 0:pmt missing
*     
*     (Creation Date and Author)
*     2002.10.23 ; First version by Y.Koshio
*     2006.04.15 ; Split sgpmtnum.F into 3 routines, i.e. into SK-I,II,III
*                     by G.Mitsuka
*     
************************************************************************

      implicit none
#include "skdetsim_parm.h"
#include "skparm.h"
#include "geopmt.h"
#include "skhead.h"


      character*4 names
      integer n5,n6,n7,n8,nopmt,iflag

      IF(NAMES.EQ."PMPW") THEN
         NOPMT = N6 + (N7-1)*51

      ELSE IF(NAMES.EQ."PMPT") THEN
         NOPMT = ISEQT(N7,N8)
         
      ELSE IF(NAMES.EQ."PMPB") THEN
         NOPMT = ISEQB(N7,N8)
      ENDIF
      
      IFLAG=1
      IF(NOPMT.EQ.0) THEN
         IFLAG=0
      ENDIF
      
      RETURN
      END
