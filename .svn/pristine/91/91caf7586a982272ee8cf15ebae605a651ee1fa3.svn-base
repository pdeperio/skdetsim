*     
      subroutine dsdigi_anti_sk4
*     
*     purpose: Select hits in the timing window and
*     put them to TQAZ bank.
*     
*     Input  : SKTQAZ common block
*     output : TQAZ bank
*     
*     2008/07   First version by Jun Kameda 
*     
      implicit none
#include "skparm.h"
#include "skhead.h"
C     #include "sktq.h"
#include "sktq_raw.h"
#include "skdetsim_parm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "zbs.h"
#include "geant321/gcflag.inc"

      integer itaskz_tmp(30*MAXPMA),iqaskz_tmp(30*MAXPMA)
      integer icabaz_tmp(30*MAXPMA)
      real    taskz_tmp(30*MAXPMA),qaskz_tmp(30*MAXPMA)
      integer nhitaz_tmp
      integer ihtflz_tmp

      integer nblock(30*MAXPMA), nblock_current, nhit_current, nhit_block(100)
      integer ntrg_min, ntrg_max

c     TQAZ bank

      integer iii,jjj
      integer ndata,isegm, idata(MAXPMA*20), ierr
      integer num_block
      integer itmin, itmax, it_diff
      integer iblock, ihit
      integer it_tmp,it_upper

c     tqareal bank

      integer ldat,ihitz,nent
      integer ibuf_tqareal(7),aibufz(3,MAXPMA*30)
      real    rbuf_tqareal(7),arbufz(3,MAXPMA*30)

      equivalence(ibuf_tqareal(1),rbuf_tqareal(1))
      equivalence(aibufz(1,1),arbufz(1,1))

      integer od_start, od_end
      
ccc   temporal to test
      integer it1,flg, i
      real t_tmp

      call vzero(nhit_block,100)
      do i = 1, maxpma*30
         aibufz(1,i) = 0
         aibufz(2,i) = 0
         aibufz(3,i) = 0
         arbufz(1,i) = 0.
         arbufz(2,i) = 0.
         arbufz(3,i) = 0.

      enddo


      if (SK_GEOMETRY .ge. SK_IV) then

c     c   taskz is measred in GEANT time to avoid rounding error.
c     c   itaskz is measured in Super-K time
c     tmin = real(it0sk)*tdc_unit_sk4 -  5000.0
c     tmax = real(it0sk)*tdc_unit_sk4 + 35000.0

         itmin = it0sk-int(5000.0/tdc_unit_sk4) !  -5 micro seccnd
#ifdef NEUTRON
         itmax = it0sk+int(535000.0/tdc_unit_sk4) ! +35 micro second
#else
         itmax = it0sk+int(35000.0/tdc_unit_sk4) ! +35 micro second
#endif

         ntrg_min = int(itmin/2**15)
         ntrg_max = int(itmax/2**15) 

         nhitaz_tmp = 0

         do jjj = 1, nhitaz_raw
            
            nblock_current =  int(itaskz_raw(jjj)/2**15)

c     c  Select hit PMTs inside the data blocks 
c     c which include the 40micro sec. time window

            if (
     &           (nblock_current .ge. ntrg_min) .and. (nblock_current .le. ntrg_max)) then

               nhitaz_tmp = nhitaz_tmp + 1
               taskz_tmp(nhitaz_tmp) =  taskz_raw(jjj)
               qaskz_tmp(nhitaz_tmp) =  qaskz_raw(jjj)
               itaskz_tmp(nhitaz_tmp) = itaskz_raw(jjj)
               iqaskz_tmp(nhitaz_tmp) = iqaskz_raw(jjj)
               icabaz_tmp(nhitaz_tmp) = icabaz_raw(jjj)
               nblock(nhitaz_tmp)     = nblock_current

               nhit_block(nblock_current - ntrg_min +1)
     &              =  nhit_block(nblock_current - ntrg_min +1) + 1

               flg = 0
               if (
     $              (itaskz_raw(jjj) .ge. itmin) .and.
     $              (itaskz_raw(jjj) .le. itmax)
     $              ) then ! in-gate flag
c     iqaskz_tmp(nhitaz_tmp) = iqaskz_tmp(nhitaz_tmp) + 2**11
                  iqaskz_tmp(nhitaz_tmp) = ibset(iqaskz_tmp(nhitaz_tmp) ,11)
                  flg = 1
               endif
            endif

         enddo

         nhitaz_raw = nhitaz_tmp
         call vzero(itaskz_raw,30*MAXPMA)
         call vzero(iqaskz_raw,30*MAXPMA)
         call vzero(icabaz_raw,30*MAXPMA)
         call vzero( taskz_raw,30*MAXPMA)
         call vzero( qaskz_raw,30*MAXPMA)

         call ucopy(itaskz_tmp, itaskz_raw,nhitaz_raw)
         call ucopy(iqaskz_tmp, iqaskz_raw,nhitaz_raw)
         call ucopy(icabaz_tmp, icabaz_raw,nhitaz_raw)
         call ucopy(taskz_tmp,   taskz_raw,nhitaz_raw)
         call ucopy(qaskz_tmp,   qaskz_raw,nhitaz_raw)


c     c  create TQAZ bank

         call kzbcr0('TQAZ',IERR)

         if (IERR .ne. 0) then
            write(*,*) 'Could not create TQAZ bank'
            stop
         endif

         num_block =ntrg_max-ntrg_min+1
         numhwsk   = num_block

         do iblock = 1,num_block
            idata(int(iblock-1)*3+1) = nhit_block(iblock)  
            idata(int(iblock-1)*3+2) = 0  
            idata(int(iblock-1)*3+3) = ntrg_min + iblock - 1   
         enddo

         ndata = num_block * 3
         isegm = 0
         call kzrep0('TQAZ',isegm,'I',ndata,idata)

c     c   TQ info. for iblock-th HW trigger,
c     c   stored in iblock-th segment of TQAZ bank.

         do iblock = 1, num_block
            ndata = 0
            nhit_current = 0

            do ihit = 1, nhitaz_raw
               if (nblock(ihit) .eq. (ntrg_min + iblock - 1) ) then
                  nhit_current = nhit_current + 1
                  it_tmp   = iand(itaskz_raw(ihit),2**15-1) ! lower 15bits is TDC count
                  it_upper = ishft(it_tmp,16)

C
C- In order to use software trigger 20000 should be added well before.
C-  Y.Hayato ( 2009/05/06 )
C
C                  idata(2*nhit_current-1)   = icabaz(ihit) + 20000
C
                  idata(2*nhit_current-1)   = icabaz_raw(ihit)
                  idata(2*nhit_current)     = it_upper + iqaskz_raw(ihit)
               endif
            enddo

            isegm = (iblock-1)*2+1
            ndata = nhit_current*2
            if( ndata .eq. 0 )then
               ndata = 2
               idata(1) = 0
               idata(2) = 0 
            endif
            call kzrep0('TQAZ',isegm,'I',ndata,idata)

         enddo

c     c Create tqareal bank

         ibuf_tqareal(1) = nhitaz_raw
         rbuf_tqareal(2) = 2.465 ! PC2PE (taken from tqrawsk.F)
         ibuf_tqareal(3) = 2 ! Electronics version is QBEE
         ibuf_tqareal(7) = it0sk

         CALL KZBCR1('TQAREAL',7,'1I 1F 5I',nhitaz_raw,3,'1I 2F',IERR) 
         CALL KZREP1('TQAREAL',0,7,ibuf_tqareal)

         od_start = -1000 + it0sk ! Taken from set_timing_gate.F
         od_end   =  1496 + it0sk ! 

         do ihitz = 1, nhitaz_raw
            ihtflz_tmp = ishft(iqaskz_raw(ihitz),-11) ! select 15bit-11bit
            ihtflz_tmp = ishft(ihtflz_tmp,1) ! copy to 5bit-1bit
            it_diff = (ntrg+2)*2**15 - it0sk
            taskz_raw(ihitz) 
     $           = taskz_raw(ihitz) 
     $           + real(it_diff)*tdc_unit_sk4 + trg_offset +1000.0 

            t_tmp = real(itaskz_raw(ihitz)-it0sk)*tdc_unit_sk4 +1000.0

            if (  
     $           iand(ihtflz_tmp,2).eq.2 .and.
     $           itaskz_raw(ihitz).ge.od_start .and. 
     $           itaskz_raw(ihitz).le.od_end) then
               ihtflz_tmp = ihtflz_tmp + 1 ! 0bit
            endif
            

C
C- In order to use software trigger 20000 should be added well before.
C-  Y.Hayato ( 2009/05/06 )
C
C            aibufz(1,ihitz) = icabaz(ihitz) + 20000 + (ihtflz_tmp * 2**16)
C
            aibufz(1,ihitz) = icabaz_raw(ihitz) + (ihtflz_tmp * 2**16)
            arbufz(2,ihitz) = taskz_raw(ihitz)
            arbufz(3,ihitz) = qaskz_raw(ihitz)
         enddo

         if (nhitaz_raw.NE.0) THEN
            call kzldat('TQAREAL',ldat)
            nent = 3*nhitaz_raw
            call lbfcopy(aibufz,izbs(ldat+7+1),nent)
         end if

      endif ! end of SK-IV routines


      return
      end
