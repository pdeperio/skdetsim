      SUBROUTINE skmkruninf

      implicit none

#include "skhead.h"
#include "skruninf.h"
#include "geant321/gcflag.inc"
#include "mcevent.h"
#include "softtrg_tblF.h"

      INTEGER*4 NDATA,IDATA(1000),IDATA1(1000)
      INTEGER*4 ITMP(3)
      INTEGER*4 IERR,LBANK
      INTEGER*4 ICOUNT, OSUB
      CHARACTER*8 CD,CT,DUM1,DUM2
      integer isegm, i

      call vzero(idata,1000)
      call vzero(idata1,1000)

      CALL KZBLOC('RUNINF',LBANK)
      IF (LBANK.eq.0) GOTO 150

      isegm = 0
      CALL KZGET0('RUNINF',isegm,NDATA,IDATA)
      isegm = 1
      CALL KZGET0('RUNINF',isegm,NDATA,IDATA1)
      CALL KZBDEL('RUNINF')

 150  CALL KZBCR0('RUNINF',IERR)
      IF (IERR.ne.0) THEN
         write(*,*) "Could not create RUNINF BANK"
      ENDIF

      isegm = 0
      ndata = SK_GEOMETRY   
      idata(1) = 0
      idata(2) = 999999    ! = nrunsk

      CALL KZREP0('RUNINF',isegm,'I',ndata,IDATA)

      
      do i = 0,31

          if (swtrgtype(i) .eq. 0) then  ! LE trigger
             idata1(90+i*5+1) = 0  ! detector
             idata1(90+i*5+2) = 47 ! threshold
             idata1(90+i*5+3) = 0   ! offset 
             idata1(90+i*5+4) = 9600  ! pre_t0 
             idata1(90+i*5+5) = 67200 ! post_t0
          endif

          if (swtrgtype(i) .eq. 1) then  ! HE trigger
             idata1(90+i*5+1) = 0  ! detector
             idata1(90+i*5+2) = 50 ! threshold
             idata1(90+i*5+3) = 0   ! offset 
             idata1(90+i*5+4) = 9600  ! pre_t0 
             idata1(90+i*5+5) = 67200 ! post_t0
          endif


          if (swtrgtype(i) .eq. 2) then  ! SLE trigger
             idata1(90+i*5+1) = 0  ! detector
             idata1(90+i*5+2) = 34 ! threshold
             idata1(90+i*5+3) = 384  ! offset 
             idata1(90+i*5+4) = 960  ! pre_t0 
             idata1(90+i*5+5) = 1923 ! post_t0
          endif

          if (swtrgtype(i) .eq. 3) then  ! OD trigger
             idata1(90+i*5+1) = 0     ! detector
             idata1(90+i*5+2) = 22    ! threshold
             idata1(90+i*5+3) = 100   ! offset 
             idata1(90+i*5+4) = 9600  ! pre_t0 
             idata1(90+i*5+5) = 67200 ! post_t0
          endif

      enddo
               idata1(251) = 0 ! mask

      isegm = 1
      ndata = 251

      CALL KZREP0('RUNINF',isegm,'I',ndata,IDATA1)
      
      RETURN

      END


