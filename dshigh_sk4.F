************************************************************************
*     -------------------------------------
      SUBROUTINE DSHIGH_SK4
*     -------------------------------------
*
*     (Purpose)
*        Calculate T, Q of each hit and fill common variables for 
*        TRACK_PHOTON = 2(old AP mode) or 4(saving memory mode).
*
*        This routine uses a simple simulation.
*         (1) Timing shift due to dark hit is not taken into account.
*         (2) Charge integration efficiency in QTC gate is not taken into account.
*        If you need precise simulation, Please use other photon tracking mode.
*
*     (Creation Date and Author)
*       2008. 7.22 : First version based on dshigh.F
*                    Jun Kameda
*       2009.09.02 : Fix bugs and remove unnecessary lines.
*
************************************************************************
      IMPLICIT NONE

#include "skparm.h"
#include "skhead.h"
C#include "sktq.h"
#include "sktq_raw.h"
#include "dsparm.h"
#include "dsparma.h"
#include "skmc.h"
#include "dsdata.h"
#include "dsoffset.h"
#include "skmisch.h"
#include "skdetsim_parm.h"

#include <zbs.h>

#ifdef FADC      
#include "fadcwrbnk.h"
#endif

      INTEGER*4 JBUF(3),KBUF(MAXRNG*2+3)
      REAL*4    SBUF(3)
      EQUIVALENCE(JBUF(1),SBUF(1))

      INTEGER*4 ICTQ(3,MXPHDS),IFLAG(MXGATE)
      REAL*4    RCTQ(3,MXPHDS)
      EQUIVALENCE(ICTQ,RCTQ)

      Integer offset, i
      Integer tbuf(20000)

C-- For KZGET

      INTEGER*4 ISEGM,ISSEGM,NDATA

C-- 
      INTEGER*4 IRING,ICNT,IPH,IC,NDARK

      INTEGER*4 IPMT,IQTMP,IFTHR

      REAL*4    T,Q,TMIN,TMAX,TSIGE,TSIGL,PE1
      REAL*4    EXPD
      REAL*4    PE_SATU

      REAL*4    DUMMY,TRND

      INTEGER*4 ITMP,JJJ

c      integer*4 ihtiflz_tmp

**   For SK-IV electronics

      integer iq_small, iq_medium, iq_large, flag_gain,ngate
      real pc2pe_mc, cnt2pc_s, cnt2pc_m, cnt2pc_l
      parameter(pc2pe_mc=2.465)
      parameter(cnt2pc_s=0.100, cnt2pc_m=0.732, cnt2pc_l=5.141)

      integer index_gate(MXGATE)
      integer itiskz_tmp(MXGATE),iqiskz_tmp(MXGATE)
      integer icabiz_tmp(MXGATE)
      real    tiskz_tmp(MXGATE),qiskz_tmp(MXGATE)


      If (SK_GEOMETRY.lt.SK_IV) then   
         write(*,*) 'DSHIGH_SK4 is used. This routine is only for SK-IV or later.
     &   Please Check SK_GEOMETRY. Stop'
         stop
      endif

      call vzero(itiskz_tmp,MXGATE)
      call vzero(iqiskz_tmp,MXGATE)
      call vzero(tiskz_tmp,MXGATE)
      call vzero(qiskz_tmp,MXGATE)
      call vzero(icabiz_tmp,MXGATE)

      COMOFST=0.

C --- Read data in MCRAWTQ' bank

      isegm = 0
      issegm = 0
      CALL KZGET2('MCRAWTQ',isegm,issegm,NDATA,KBUF)

      NRNGDS = KBUF(1)
      NPHADS = 0
      DO 1 IRING=1,NRNGDS
         IPARDS(IRING) = KBUF(2*IRING)
         NPHRDS(IRING) = KBUF(2*IRING+1)
         NPHADS = NPHADS+NPHRDS(IRING)
 1    CONTINUE

      COMOFST=OFFSET

*     TRG 48clk. info.
      nt48sk(1) = 0
      nt48sk(2) = 0
      nt48sk(3) = ifix((offset0-offset)/20.)

      ICNT = 0
      CALL VZERO(IFLAG,MAXPM)

      DO 10 IRING=1,NRNGDS
         DO 20 IPH=1,NPHRDS(IRING)

            CALL KZGET2('MCRAWTQ',IRING,IPH,NDATA,JBUF)

            ICNT = ICNT + 1
            IC = JBUF(1)
            T  = SBUF(2) 
            Q  = FLOAT(JBUF(3))
            ICTQ(1,ICNT) = IC
            RCTQ(2,ICNT) = T
            RCTQ(3,ICNT) = Q
            IFLAG(IC)    = ICNT

 20      CONTINUE
 10   CONTINUE   

cc Generate dark hits

      tmin = -trg_offset - twid_block*2
#ifdef NEUTRON
      tmax = twid_block
#else
      tmax = (twid_block - trg_offset) + twid_block*3
#endif
      EXPD = DARKDS*(TMAX-TMIN)*(1.E-09)

      DO 110 IPMT=1,MAXPM

       if(imis(IPMT).ne.0) goto 110 ! skip non-existing PMTs

      call gpoiss(expd,ndark,1)

      if (ndark .gt. 0) then

         do 100 i = 1,ndark

         ICNT = ICNT + 1
          
         IF(ICNT.GT.MXPHDS) THEN
            WRITE(6,*) 'ERROR:  IN DSHIGH_SK4'
            WRITE(6,*) 'Too much hits. Stop'
            STOP
         ENDIF

         call grndm(dummy,1)
         ICTQ(1,ICNT) = IPMT
         RCTQ(2,ICNT) = DUMMY*(TMAX-TMIN)+TMIN
         RCTQ(3,ICNT) = 1

 100     continue
      endif
 110  continue

cc  sort array

      IF(ICNT.GE.1) THEN
         CALL SORTI(RCTQ,3,ICNT,1)
      END IF

cc Consider resolution

      ngate = 0

      DO 30 IC=1,ICNT   ! loop over hits

*        skip not-existing PMTs
         if(imis(ICTQ(1,IC)).ne.0) goto 30

         IPMT = ICTQ(1,IC)
         T    = RCTQ(2,IC)
         Q    = RCTQ(3,IC)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-IV (Just copy from SK-III)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          IQTMP = int(Q +.1)
          Q=0.

          do 38 ITMP=1,IQTMP
             CALL SKRN1PE_SK3(IPMT,PE1)
             Q=Q+PE1
 38       continue

             CALL DSTHR_SK4(IPMT,Q,IFTHR)
             IF(IFTHR.EQ.1) GOTO 30
             call dssatu_sk4(Q,pe_satu)  ! saturation effect of SK-IV
             Q = PE_satu

             CALL SKTSIG_sk4_asy(Q,TSIGE,TSIGL) 
             t =t + 4.4*exp(-7*q**2)+0.7
             CALL rntwogauss_asy(t, TSIGE,TSIGL, trnd)
             T = TRND

cc     fill variables in SKTQ bank

         ngate = ngate + 1
         ICABIZ_tmp(ngate) = IPMT
         QISKZ_tmp(ngate)  = Q

cc   TISKZ is measured in GEANT time. ITISKZ is measred in Super-K time
cc   to avoid rounding error of TISK.

         TISKZ_tmp(ngate)  =  T
         itiskz_tmp(ngate) = int((t+trg_offset)/tdc_unit_sk4) + (ntrg+2)*2**15 

c  Q digitization.
c  A very simple approximation.
c  More precise calibration results should be put here.

         iq_small  = int(Q/pc2pe_mc/cnt2pc_s) + 961
         iq_medium = int(Q/pc2pe_mc/cnt2pc_m) + 961
         iq_large  = int(Q/pc2pe_mc/cnt2pc_l) + 961

c  Decide the ADC gain level

           flag_gain = 0
           IQISKZ_TMP(ngate) = iq_small 

        if (iq_small  .ge. 1350) then
           flag_gain = 1
           IQISKZ_TMP(ngate) = iq_medium
        endif

        if (iq_medium .ge. 1350) then
           flag_gain = 2
           IQISKZ_TMP(ngate) = iq_large 
        endif

        IQISKZ_TMP(ngate) = IQISKZ_TMP(ngate) + flag_gain*2**14 ! range flag
 
 30   CONTINUE  ! end of the loop over hits

cc  Sorting variables.

       nqisk_raw = ngate
       call sortzv(itiskz_tmp,index_gate,ngate,1,0,0)

       do jjj = 1,ngate
          itiskz_raw(jjj) = itiskz_tmp(index_gate(jjj))
          iqiskz_raw(jjj) = iqiskz_tmp(index_gate(jjj))
          icabbf_raw(jjj) = icabiz_tmp(index_gate(jjj))
          icabbf_raw(jjj) = iand(icabbf_raw(jjj),65535)
          icabbf_raw(jjj) 
     $         = icabbf_raw(jjj)
     $         +(iand(ishft(iqiskz_raw(jjj),-11),62)*65536)
           tbuf_raw(jjj) =  tiskz_tmp(index_gate(jjj))
           qbuf_raw(jjj) =  qiskz_tmp(index_gate(jjj))
       enddo

      RETURN
      END
