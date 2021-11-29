************************************************************************
*     -------------------------------------
      SUBROUTINE FADCDARK
*     -------------------------------------
*
*     (Purpose)
*        Dark noise of PMT's for FADC
*        Note: the dark noise for FADC is not identical with original.
*     (Creation Date and Author)
*        2003.06.10 First version     M.Ishitsuka
*        2004.04.30 Added to CVS      M.Ishitsuka
************************************************************************

#include "skparm.h"
#include "skhead.h"
#include "dsparm.h"
#include "dsparma.h"
#include "skmc.h"
#include "dsdata.h"
#include "dsoffset.h"
#include "skmisch.h"

#include "fadcparm.h"
#include "fadctabl.h"
#include "fadc1pe.h"

      real darkexp

      darkexp = DARKDS * NFADC_RANGE * 2 * 1.E-09

      DO 30 IC=1,MAXPM

*        skip not-existing PMTs
         if(imis(ic).ne.0) goto 30

         CALL GRNDM(RDUMMY,1)
         IF(RDUMMY.GE.darkexp) GOTO 30
         CALL GRNDM(RDUMMY,1)
         T = RDUMMY*NFADC_RANGE*2
         call fadcrecord(IC,T)

 30   CONTINUE


      RETURN
      END
