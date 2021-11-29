************************************************************************
*     ---------------------
      SUBROUTINE FADCWRINIT
*     ---------------------
*     Initialize the common : 
*       RFCHPLS,nwrchid,nwrhit,nwrtime,nwrsize,nwrpulse
*     
*     2001.09.16 ; T.Ooyabu
*     2001.11.05 ; added commons for fadcwrite      T.Ooyabu
*     2002.08.31 ; set parameters for FADC bank     T.Ooyabu
*     2003.05.27 ; change the bank format to consistent with observed data
*                                                   M.Ishitsuka
************************************************************************
      IMPLICIT NONE

#include "fadcparm.h"
#include "fadcbnk.h"
#include "fadc1pe.h"      /* RFCHPLS */
#include "fadcwrbnk.h"    /* nwrthreshold,nwrpedestal,rwrcount2mv */

      Integer ifcall

      if (ifcall.eq.1) goto 100
      ifcall = 1
*--   Set parameters for FADC bank
      nthreshold = nwrthreshold
      npedestal  = nwrpedestal
      count2mv   = rwrcount2mv
      if (nthreshold.lt.npedestal) then
         write(6,*)'nthreshold should be larger than npedestal'
         stop
      end if
      write(6,'(A35)')'====== FLASH ADC RUN CONTROL ======'
      write(6,'(A35)')'                                   '
      write(6,'(A15,I5)')'nthreshold :',nthreshold
      write(6,'(A15,I5)')'npedestal  :',npedestal
      write(6,'(A15,F5.3)')'count2mv  :',count2mv
      write(6,'(A35)')'                                   '
      write(6,'(A35)')'==================================='
 100  continue

      Call VZERO(RFCHPLS,NFADC_RANGE*MAXFADC)

      RETURN
      END


