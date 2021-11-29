*********************************************************************
*
*     -----------------
      SUBROUTINE SGGEOM
*     -----------------
*
*     (Purpose)
*       Define the detector geometry and tracking media.
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       2003.01.07 ; First version by Y.Koshio
*       2006.04.29 ; Add Initialization of the Black Sheet Position
*                    for SK-II&SK-III
*       2007.09.06 ; Fix Black Sheet position bug
*       2008.09.30 ; Add SK_GEOMETRY = SK_IV  by J.Kameda
************************************************************************

      implicit none

#include "skhead.h"
#include "skdetsim_parm.h"
#include "guskge.h"

      IF(SK_GEOMETRY.eq.SK_I) then !Good Old SK-I
         call sggeom_sk1
      ELSE IF(SK_GEOMETRY.eq.SK_II) then !SK-II
         call sggeom_sk2
         write(*,*) 'GEOMETRY is set to SK_II'
      ELSE IF(SK_GEOMETRY.eq.SK_III) then !SK-III
         call sggeom_sk3
         write(*,*) 'GEOMETRY is set to SK_III'
      ELSE IF(SK_GEOMETRY.eq.SK_IV) then !SK-IV
         call sggeom_sk3
         write(*,*) 'GEOMETRY is set to SK_IV (same as SK-III)'
      ELSE IF(SK_GEOMETRY.eq.SK_V) then !SK-V
         call sggeom_sk3
         write(*,*) 'GEOMETRY is set to SK_V (same as SK-III)'
      ELSE                      !SK-IV not made yet         
         write(*,*) 'Error SGGEOM wrong SK_GEOMETRY flag',SK_GEOMETRY
         stop
      ENDIF
      return
      end
