      SUBROUTINE SKTRGINFO

      implicit none

#include "skparm.h"
#include "skhead.h"
#include "skdetsim_parm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include <zbs.h>

#include "geant321/gcflag.inc"

      INTEGER*4 NDATA0,NDATA1,IDATA1(1000),IDATA0(1000)
      INTEGER*4 ITMP(3)
      INTEGER*4 IERR, LBANK, lsegm
      INTEGER*4 ICOUNT, OSUB
      CHARACTER*8 CD,CT,DUM1,DUM2

      CALL KZBLOC('HEADER',LBANK)
      IF (LBANK.eq.0) then
         write(*,*) 'SKTRGINFO: No Header bank was found. Someting Strainge.'
c         stop
      endif


      CALL KZGET0('HEADER',0,NDATA0,IDATA0)
      CALL KZGET0('HEADER',1,NDATA1,IDATA1)
      CALL KZBDEL('HEADER')

 150  CALL KZBCR0('HEADER',IERR)
      IF (IERR.ne.0) THEN
         write(*,*) "Could not create HEADER BANK"
         stop
      ENDIF

      call KZREP0('HEADER',0,'I',NDATA0,IDATA0)
      
      if (SK_GEOMETRY .ge. SK_IV) then
       lsegm     = 19
       idata1(10) = idtgsk   ! idtgsk
       idata1(11) = ntrg+2   ! nevhwsk
       idata1(12) = 0        ! hidhsk
       idata1(13) = 0        ! hidlsk
       idata1(14) = it0sk    ! it0sk   --> dsswtrg
       idata1(15) = 76800    ! gatewsk , 9600+ 67200
       idata1(18) = numhwsk  ! numhwsk --> dsdigi_sk4
       idata1(19) = 0        ! ifevsk
      else
       lsegm     = NDATA1
      endif

      CALL KZREP0('HEADER',1,'I',lsegm,IDATA1)

      RETURN
 100  FORMAT(I2,A1,I2,A1,I2)

      END


