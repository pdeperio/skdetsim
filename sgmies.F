************************************************************************
*     ----------------------
      REAL FUNCTION SGMIES(COSTH)
*     ----------------------
*
*     (Purpose)
*        Get the angle dist. of MIE scattering
*
*     (Input)
*        COSTH  : input angle
*
*     (Output)
*        SGMIES : angle dist. (max = 1)
*
*     (Creation Date and Author)
*       1997. 2.14 ; defined by M.Shiozawa
*       2003. 5.16 ; modified by J.Yoo
*       2003 .6.28 ; C. Sterner : replace SK-I tune removed by Yoo
*       2006. 4.18 ; K Bays : added SK-III case using SK-II data
************************************************************************
      IMPLICIT NONE
#include "skdetsim_parm.h"
#include "skhead.h"

      INTEGER IBIN_SK
      INTEGER IBINLE_SKI,IBINAP_SKI

      PARAMETER(IBINLE_SKI = 6,IBINAP_SKI=3)
      PARAMETER(IBIN_SK = 3)

      REAL    costh

      REAL    COSTLST_SK(IBIN_SK),MIELST_SK(IBIN_SK)
      REAL    COSTLSTLE_SKI(IBINLE_SKI),MIELSTLE_SKI(IBINLE_SKI)
      REAL    COSTLSTAP_SKI(IBINAP_SKI),MIELSTAP_SKI(IBINAP_SKI)
      INTEGER i
      REAL    rcost

      DATA COSTLSTLE_SKI /-1.,0.,0.866,0.940,0.984,1./
      DATA MIELSTLE_SKI  /1.0,1.0,1.5,2.,50.,100./
CCS  This angle distribution means front high peak
      DATA COSTLSTAP_SKI /-1.,0.,1./
      DATA MIELSTAP_SKI  /0.,0.,100./

*
* The Mie scattering angle distributions for 
* LE_TUNE mode and AP_TUNE mode are the same from (2003/05/16) for sk-II
* But for the purpose of the separation of those two mode, keep the routine itself.
*

      DATA COSTLST_SK /-1.,0.,1./
      DATA MIELST_SK  /0.,0.,100./

      IF(SK_GEOMETRY.ge.SK_II) THEN
         DO I=2,IBIN_SK
            IF(COSTH.LE.COSTLST_SK(I)) THEN
               RCOST = (COSTH - COSTLST_SK(I-1)) /
     /              (COSTLST_SK(I)-COSTLST_SK(I-1))
               SGMIES = MIELST_SK(I-1)*(1.-RCOST) +
     +              MIELST_SK(I)*RCOST
               GOTO 101
            END IF
         END DO

      ELSE IF(SK_GEOMETRY.eq.SK_I) THEN
         IF ((WATER_TUNE .EQ. LE_TUNE) .or. (WATER_TUNE .eq. UNI_TUNE)
     &		.or.(WATER_TUNE.eq.10) .or. (WATER_TUNE.eq.11) ) THEN
            DO I=2,IBINLE_SKI
               IF(COSTH.LE.COSTLSTLE_SKI(I)) THEN
                  RCOST = (COSTH - COSTLSTLE_SKI(I-1)) / 
     /                 (COSTLSTLE_SKI(I)-COSTLSTLE_SKI(I-1))
                  SGMIES = MIELSTLE_SKI(I-1)*(1.-RCOST) + 
     +                 MIELSTLE_SKI(I)*RCOST
                  GOTO 101
               END IF
            END DO
         ELSE
            DO I=2,IBINAP_SKI
               IF(COSTH.LE.COSTLSTAP_SKI(I)) THEN
                  RCOST = (COSTH - COSTLSTAP_SKI(I-1)) /
     /                 (COSTLSTAP_SKI(I)-COSTLSTAP_SKI(I-1))
                  SGMIES = MIELSTAP_SKI(I-1)*(1.-RCOST) + 
     +                 MIELSTAP_SKI(I)*RCOST
                  GOTO 101
               END IF
 2000       END DO
         END IF
         
      ENDIF
      
 101  CONTINUE

      SGMIES = SGMIES/100.
         
        
      RETURN
      END
