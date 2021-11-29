*
      subroutine dstrgnum
*
*     Set dummy TRG count number of the event and 
*     relative timing of the event to the data block.
*
*     input  : none
*
*     output : ntrg       : dummy TRG count number (16bit)
*              trg_offset : time from start point of the data block
*                           to GEANT zero point.
*
*     2008/07/04  Jun Kameda for SK-4
*

      implicit none

#include "geant321/gcflag.inc"
#include "skparm.h"
#include "skmc.h"
#include "dsparm.h"
#include "dsdata.h"

      integer itrg_clock,ndummy
      real dummy
c
c --- Set relative time to the start time of the data block.
c --- Every data block will be get every 1024 count of the 
c --- 60MHz master clock = 17.066 usec.
c 
         call grndm(dummy,1)
c         itrg_clock = int(dummy*1024)
c         trg_offset = real(itrg_clock*tdc_unit_sk4)
         trg_offset = real(dummy*tdc_unit_sk4*32768.0)

c --- calculate block number 

         ndummy = int(ievent*100)
*         ntrg = iand(65535,ndummy)  ! dummy TRG count is 16bit.
         ntrg = iand(32767,ndummy)  ! dummy TRG count is 16bit.
         nblock_mc = 5

         return
         end

