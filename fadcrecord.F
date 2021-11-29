************************************************************************
*     ------------------------------------
      SUBROUTINE FADCRECORD(NOPMT,TOFGSK)
*     ------------------------------------
*     2001.09.16 ; First version by Takato Ooyabu
*     2001.10.25 ; set Saturation by Takato Ooyabu
*     Variables
*       NFADC_CH : FADC channel of hit PMT
*       IBFADC   : Bin Number
************************************************************************
      IMPLICIT NONE

#include "skparm.h"

#include "fadcparm.h"
#include "fadctabl.h"
#include "fadc1pe.h"

C     Arguments
      INTEGER NOPMT
      REAL TOFGSK
C --  Static variables
      INTEGER NFADC_CH,IBFADC
      REAL frac
      INTEGER ii

c     HIT PMT # to FADC channel
      NFADC_CH=IFADC(NOPMT)
      
      frac = -1.
      Do while ( frac .LE. 0. )
         Call skrn1pe(frac)
      End do

C     Calculate RFCHPLS
      DO ii=1,50
         IBFADC = int((ii+TOFGSK)*NS2BIN) 
         If (IBFADC.GT.0 .AND. IBFADC.LE.NFADC_RANGE ) then
            RFCHPLS(IBFADC,NFADC_CH)=
     $           RFCHPLS(IBFADC,NFADC_CH) + frac*FADC_AMP_FAC*RPACKET(ii)
* --        if pulse height is over 1000.mV
            if ( RFCHPLS(IBFADC,NFADC_CH) .gt. 1000. ) then
               RFCHPLS(IBFADC,NFADC_CH)=1000.
            endif
* --
         Endif
*         write(*,*) RFCHPLS(IBFADC,NFADC_CH)
      END DO

      RETURN
      END
