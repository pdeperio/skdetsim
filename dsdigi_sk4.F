*     
      subroutine dsdigi_sk4
*     
*     purpose: Select hits in the timing window and
*     put them to TQZ bank.
*     
*     Input  : SKTQZ common block
*     output : TQZ bank
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

      integer itiskz_tmp(30*MAXPM),iqiskz_tmp(30*MAXPM)
      integer icabiz_tmp(30*MAXPM)
      real    tiskz_tmp(30*MAXPM),qiskz_tmp(30*MAXPM)
      integer nqiskz_tmp

C-- only for compatibility -- ( Y.Hayato / 2009/05/06 )
      integer iab(MAXPM)

      integer nblock(30*MAXPM), nblock_current, nhit_current, nhit_block(100)
      integer ntrg_min, ntrg_max

      integer ihtiflz_tmp

c     TQZ bank

      integer iii,jjj
      integer ndata,isegm, idata(MAXPM*20), ierr
      integer num_block
      integer itmin, itmax, it_diff
      integer iblock, ihit
      integer it_tmp,it_upper

c     TQREAL bank

      integer ldat,ihitz,nent
      integer ibuf_tqreal(7),ibufz(3,MAXPM*30)
      real    rbuf_tqreal(7),rbufz(3,MAXPM*30)

      equivalence(ibuf_tqreal(1),rbuf_tqreal(1))
      equivalence(ibufz(1,1),rbufz(1,1))

      integer id_start, id_end
      
ccc   temporal to test
      integer it1,flg
      real t_tmp

      call vzero(nhit_block,100)

      if (SK_GEOMETRY .ge. SK_IV) then

c     c   Tiskz is measred in GEANT time to avoid rounding error.
c     c   ITISKZ is measured in Super-K time
c     tmin = real(it0sk)*tdc_unit_sk4 -  5000.0
c     tmax = real(it0sk)*tdc_unit_sk4 + 35000.0

         itmin = it0sk-int(5000.0/tdc_unit_sk4) !  -5 micro seccnd
#ifdef NEUTRON
         itmax = it0sk+int(535000.0/tdc_unit_sk4) ! +35 micro second
         print *,"NEUTRON setting ", itmax
#else
         itmax = it0sk+int(35000.0/tdc_unit_sk4) ! +35 micro second
#endif
         ntrg_min = int(itmin/2**15)
         ntrg_max = int(itmax/2**15) 

         nqiskz_tmp = 0

         do jjj = 1, nqisk_raw
            
            nblock_current =  int(itiskz_raw(jjj)/2**15)

c     c  Select hit PMTs inside the data blocks 
c     c which include the 40micro sec. time window

            if (
     &           (nblock_current .ge. ntrg_min) 
     $           .and. (nblock_current .le. ntrg_max)) then

               nqiskz_tmp = nqiskz_tmp + 1
               tiskz_tmp(nqiskz_tmp)  =  tbuf_raw(jjj)
               qiskz_tmp(nqiskz_tmp)  =  qbuf_raw(jjj)
               itiskz_tmp(nqiskz_tmp) = itiskz_raw(jjj)
               iqiskz_tmp(nqiskz_tmp) = iqiskz_raw(jjj)
               icabiz_tmp(nqiskz_tmp) 
     $              = iand(icabbf_raw(jjj),65535)
               nblock(nqiskz_tmp)     = nblock_current

               nhit_block(nblock_current - ntrg_min +1)
     &              =  nhit_block(nblock_current - ntrg_min +1) + 1

               flg = 0
               if ((itiskz_raw(jjj) .ge. itmin) 
     $              .and. (itiskz_raw(jjj) .le. itmax)) 
     $              then ! in-gate flag
c                 iqiskz_tmp(nqiskz_tmp) = iqiskz_tmp(nqiskz_tmp) + 2**11
                  iqiskz_tmp(nqiskz_tmp) = ibset(iqiskz_tmp(nqiskz_tmp) ,11)
                  flg = 1
               endif
            endif

         enddo

         nqisk_raw = nqiskz_tmp
         call vzero(itiskz_raw,30*MAXPM)
         call vzero(iqiskz_raw,30*MAXPM)
         call vzero(icabbf_raw,30*MAXPM)
         call vzero( tbuf_raw,30*MAXPM)
         call vzero( qbuf_raw,30*MAXPM)

         call ucopy(itiskz_tmp, itiskz_raw,nqisk_raw)
         call ucopy(iqiskz_tmp, iqiskz_raw,nqisk_raw)
         call ucopy(icabiz_tmp, icabbf_raw,nqisk_raw)
         call ucopy(tiskz_tmp,  tbuf_raw,  nqisk_raw)
         call ucopy(qiskz_tmp,  qbuf_raw,  nqisk_raw)

c     c for check

c     do iii = 1, nqiskz
c     write(79,*) ievent,it0sk,itiskz(iii)
c     enddo

c     c  create TQZ bank

         call kzbcr0('TQZ',IERR)

         if (IERR .ne. 0) then
            write(*,*) 'Could not create TQZ bank'
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
         call kzrep0('TQZ',isegm,'I',ndata,idata)

c     c   TQ info. for iblock-th HW trigger,
c     c   stored in iblock-th segment of TQZ bank.

         do iblock = 1, num_block
            ndata = 0
            nhit_current = 0

            do ihit = 1, nqisk_raw
               if (nblock(ihit) .eq. (ntrg_min + iblock - 1) ) then
                  nhit_current = nhit_current + 1
                  it_tmp   = iand(itiskz_raw(ihit),2**15-1) ! lower 15bits is TDC count
                  it_upper = ishft(it_tmp,16)

                  idata(2*nhit_current-1)   = icabbf_raw(ihit)
                  idata(2*nhit_current)     = it_upper + iqiskz_raw(ihit)
               endif
            enddo

            isegm = (iblock-1)*2+1
            ndata = nhit_current*2
            call kzrep0('TQZ',isegm,'I',ndata,idata)

         enddo

c     c Create TQREAL bank

         ibuf_tqreal(1) = nqisk_raw
         rbuf_tqreal(2) = 2.465 ! PC2PE (taken from tqrawsk.F)
         ibuf_tqreal(3) = 2 ! Electronics version is QBEE
         ibuf_tqreal(7) = it0sk

         CALL KZBCR1('TQREAL',7,'1I 1F 5I',NQISK_raw,3,'1I 2F',IERR) 
         CALL KZREP1('TQREAL',0,7,ibuf_tqreal)

         call vzero(iab,MAXPM)
         id_start = -1000 + it0sk ! Taken from set_timing_gate.F
         id_end   =  1496 + it0sk ! 

         do ihitz = 1, nqisk_raw
            ihtiflz_tmp = ishft(iqiskz_raw(ihitz),-11) ! select 15bit-11bit
            ihtiflz_tmp = ishft(ihtiflz_tmp,1)         ! copy to 5bit-1bit
c           tiskz(ihitz) = tiskz(ihitz) - real(it0sk)*tdc_unit_sk4 + 1000.0  
            it_diff = (ntrg+2)*2**15 - it0sk
            tbuf_raw(ihitz) 
     $           = tbuf_raw(ihitz) 
     $           + real(it_diff)*tdc_unit_sk4 
     $           + trg_offset +1000.0 

            t_tmp = real(itiskz_raw(ihitz)-it0sk)*tdc_unit_sk4 +1000.0

C-- Why do you set IAB??
C
            if (icabbf_raw(ihitz).le.MAXPM .and. 
     $           iand(ihtiflz_tmp,2).eq.2 .and.
     $           iab(icabbf_raw(ihitz)).eq.0.and.
     $           ITISKZ_raw(ihitz).ge.id_start .and. 
     $           ITISKZ_raw(ihitz).le.id_end) then
               ihtiflz_tmp = ihtiflz_tmp + 1 ! 0bit
               iab(iand(icabbf_raw(ihitz),65535)) = 1
            endif


            IBUFZ(1,ihitz)    = ICABBF_RAW(ihitz) + (ihtiflz_tmp * 2**16)
            icabbf_raw(ihitz) = IBUFZ(1,ihitz)
            RBUFZ(2,ihitz)    = TBUF_RAW(ihitz)
            RBUFZ(3,ihitz)    = QBUF_RAW(ihitz)
         enddo

         IF (NQISK_raw.NE.0) THEN
            CALL KZLDAT('TQREAL',LDAT)
            nent = 3*nqisk_raw
            CALL LBFCOPY(IBUFZ,IZBS(LDAT+7+1),nent)
         END IF

      endif ! end of SK-IV routines


      return
      end
