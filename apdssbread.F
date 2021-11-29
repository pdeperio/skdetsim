************************************************************************
*     -------------------------------------
      SUBROUTINE APDSSBREAD(ILUN,IERR)
*     -------------------------------------
*
*     (Purpose)
*        Read data and make sub-event information
*
*     (Input)
*        ILUN ; Logical Unit #
*
*     (Output)
*        IERR ; error status
*               1 ; Read Error
*               2 ; End of File
*
************************************************************************
      implicit none
      Integer ILUN,IERR

#include "skparm.h"
#include "skhead.h"
#include "skbadc.h"
#include "apversion.h"
#include "apscndprt.h"

      Integer remevent
      Data remevent/0/
      Integer nevent_oya, nrun_oya
      REAL SBPOS(3,maxnsubnev),SBPMOM(3,maxnsubnev)
     &     ,SBTT(maxnsubnev)
      INTEGER SBIP(maxnsubnev), SBNVCT
      Integer nsub, nvct, i
      Integer nsevnt,isevnt

      IERR = 0

      If (remevent.gt.0) Goto 2000
*     read oya event
 2    Continue
      Call kzeclr
      Call kzread(ilun,IERR)
      Call headsk(1002)
      If (IERR.eq.0) Goto 1003
      Return
 1002 Continue
      IERR = 1
      Return
 1003 Continue
      call skverchk
      nevent_oya = nevsk
      nrun_oya = nrunsk
*     check sub-event bank
      Call kzmnum(nsevnt,isevnt)
      If (nsevnt.gt.1) then
         print *, 'WARNING: there is already sub-event bank. SKIP.'
         Goto 2
      Endif
*     check sub-event in CONVVECT bank
      Call apfill2nd(numsubev)
      If (numsubev.le.0) Goto 2
      remevent = numsubev
      Call kzeclr
      write(99,*) nrun_oya,nevent_oya,numsubev
      write(*,*) nrun_oya,nevent_oya,numsubev

 2000 Continue
*     make vector
      remevent = remevent - 1
      nsub = numsubev - remevent
*     set arguments
      SBNVCT = nevinsube(nsub)
      Do i=1,sbnvct
         nvct = nsublist(nsub,i)
         SBPOS(1,i) = vtx2nd(1,nvct)
         SBPOS(2,i) = vtx2nd(2,nvct)
         SBPOS(3,i) = vtx2nd(3,nvct)
         SBPMOM(1,i) = p2nd(1,nvct)
         SBPMOM(2,i) = p2nd(2,nvct)
         SBPMOM(3,i) = p2nd(3,nvct)
*         SBTT(i) = real(ifix(t2nd(nvct)-nevstartt(nsub)))*(1.E-09)
         SBTT(i) = real(ifix(t2nd(nvct)))*(1.E-09)
         SBIP(i) = iprt2nd(nvct)
      Enddo
      Do i=1,sbnvct
         print *, 'RING,PID,TDIF=',i,SBIP(i),SBTT(i)
      Enddo

      CALL APVCMPMK(SBNVCT,SBPOS,SBPMOM,SBTT,SBIP)

      Return
      END










