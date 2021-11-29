      subroutine sggeom_dome

      implicit none

#include "geotnk.h"
#include "gumed_u.h"
#include "sggeom_dome.h"

      REAL ppckov(2),absco(2),effic(2),rindex(2)

      INTEGER IVOLERR

      write(*,*) '========================'
      write(*,*) '  Defining SK Dome'
      write(*,*) '========================'

*     DOME AIR MEDIUM 
*     uses standard material number 15 for air.
*     first zeroes refer to mag. field parameters
*     last zeroes indicate no user parameters
      write(*,*) "Air medium"
      CALL GSTMED(IMD_DOME_AIR,'Dome Air$',15,0,0,0,180,
     &     DOME_AIR_STEMAX, DOME_AIR_DEEMAX, DOME_AIR_EPSIL,
     &     DOME_AIR_STMIN, 0, 0)
      
      write(*,*) "DOMECYL"
      CALL GSVOLU('DMCL','TUBE',IMDMST,
     &     dome_cyl_par, npar_dome_cyl, IVOLERR )
      write(*,*) "DOMECYL (ERR=",IVOLERR,")"

      write(*,*) "dome_cyl_halfz: ",dome_cyl_halfz
      write(*,*) "dome_sphe_radius: ",dome_sphe_radius
      CALL GSVOLU('DMSP','SPHE', IMD_DOME_AIR,
     &     dome_sphe_par, npar_dome_sphe, IVOLERR )
      write(*,*) "DOMESPHE (ERR=",IVOLERR,")"


*     Give the Air Cherenkov Properties
      ppckov(1) = 1.768E-09
      ppckov(2) = 4.126E-09
      absco(1)  = 6000.
      absco(2)  = 6000.
      effic(1)  = 0.05
      effic(2)  = 0.05
      rindex(1) = 1.344
      rindex(2) = 1.344
      call gsckov(IMD_DOME_AIR,2,ppckov,absco,effic,rindex)

*     POSITION VOLUMES
*     For Mother A (MAST), We place B (DOMECYL) in A. 
*     Then place the protruding C (DOMESPHE) into B
      write(*,*) "GPOS: DOMECYL"
      CALL GSPOS('DMCL',1,'MAST',
     &     dome_cyl_center(1),dome_cyl_center(2),dome_cyl_center(3),
     &     1,'MANY')
      write(*,*) "GPOS: DOMESPHE"
      CALL GSPOS('DMSP',1,'DMCL',
     &     dome_sphe_center(1),dome_sphe_center(2),dome_sphe_center(3),
     &     1,'MANY')
      write(*,*) "========END OF DOME DEF========"
      return
      end
