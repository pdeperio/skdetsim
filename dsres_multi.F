************************************************************************
*     -------------------------------------
      SUBROUTINE DSRES_MULTI
*     -------------------------------------
*
*     (Purpose)
*        Assemble the individual photons' information.
*        Consider ADC and TDC resolutions
*
*     (Input)
*        'MCTQ' bank
*
*     (Output)
*        'DATATMP' bank (will be deleted in "dstrig.F")
*        ITISKZ(i) (i = 1, to number of hit) for software trigger
*       
*     (Creation Date and Author)
*       2008.07.08 ; First version by Jun Kameda
*                    based on dsres.F
*                    Multi-hit is taken into account.
*                    SK-IV is supported.
*       2017.04.10 ; After-pulse simulation added by Y. Nakajima
*
************************************************************************

      implicit none

#include "skparm.h"
#include "skhead.h"
C#include "sktq.h"
#include "sktq_raw.h"
#include "skdetsim_parm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "geant321/gcflag.inc"
#include "dsbnkf.h"
#include "skpmt.h"

      integer nphoton
      integer ph_pmt_list(MXPHDS)
      real    ph_t_list(MXPHDS)
      integer i,j,k,kkk

      real    tlist(MXGATE) ,QLIST(MXGATE),charge_list(MXGATE)
      integer iclist(MXGATE),iqlist(MXGATE)
      integer ngate_raw,ngate

      integer IBUF(3)
      real RBUF(3)
      EQUIVALENCE (IBUF(1),RBUF(1))
      integer iq_small, iq_medium, iq_large, flag_gain
      real pc2pe_mc, cnt2pc_s, cnt2pc_m, cnt2pc_l
      parameter(pc2pe_mc=2.465)
      parameter(cnt2pc_s=0.100, cnt2pc_m=0.732, cnt2pc_l=5.141)

**
      integer index_gate(MXGATE)
      integer itiskz_tmp(MXGATE),iqiskz_tmp(MXGATE)
      integer icabiz_tmp(MXGATE)
      real    tiskz_tmp(MXGATE),qiskz_tmp(MXGATE)
      integer jjj

**

      integer iring,iph,ndata,ic,iq,ibef,ierr,ifthr,icab,
     &     nwmast,nsegm,lsegm,isseg,isegm
      real t,pe1,pe,pe_satu, thr,tsig,trnd,qpc,qbuf
      real tsigE, tsigL

      integer idata(1000)

cc
cc  added by yokozawa for Timing resolution tuned by t-TOF-T0 in Linac
cc
      real tunefac
      parameter(tunefac=0.945)

cc Variables for after-pulse generation

      real x_dummy, t_dummy
      real after_pulse_prob, after_pulse_offset, after_pulse_sigma
c     values are from Okajima-san's study of the after pulses
      parameter(after_pulse_prob = 0.0056, after_pulse_offset = 14000, after_pulse_sigma = 2100)


C --- Consider T & Q resolutions

      CALL VZERO(QLIST,MXGATE)	
      CALL VZERO(charge_list,MXGATE)	
      
cc
cc  make list
cc

      nphoton = 0
      ngate  = 0

      call vzero(ph_t_list,MXPHDS)
      call vzero(ph_pmt_list,MXPHDS)

      DO 10 IRING=1,NRNGDS+1
         DO 20 IPH=1,NPHRDS(IRING)
c            CALL KZGET2('MCTQ',IRING,IPH,NDATA,IBUF)
            CALL KZGET2('MCPHOTON',IRING,IPH,NDATA,IBUF)

            nphoton = nphoton + 1
            ph_t_list(nphoton) = RBUF(2)
            ph_pmt_list(nphoton) = IBUF(1)

c     YN: Here adding (very primitive) after-pulses simulation
            if (ADD_AFTER_PULSE .eq. 1) then
               call grndm(x_dummy,1)
               if (x_dummy .lt. after_pulse_prob) then

                  call rngausk2(after_pulse_offset,after_pulse_sigma,t_dummy)
c                  write(*,*) 'DSRES_MULTI: Generating after pulse!'
c                  write(*,*) 'DSRES_MULTI: ',RBUF(2), t_dummy

                  nphoton = nphoton + 1
                  ph_t_list(nphoton) = RBUF(2)+t_dummy
                  ph_pmt_list(nphoton) = IBUF(1)
               endif
            endif

 20      CONTINUE
 10   CONTINUE   
         
C --- delete MCPHOTON bank to save Memory size

      IF(NBNKDS(1).EQ.0) CALL KZBDEL('MCPHOTON')
      
C --- calculate number of photons in each QTC gate

      call dscnt(nphoton, ph_t_list, ph_pmt_list, ngate_raw,
     & iclist,tlist,iqlist,charge_list)

C --- Save data into a temporary bank

      If (sk_geometry.le.SK_III) then
         NWMAST = 0
         NSEGM = 1
         LSEGM = 3
         CALL KZBCR2('DATATMP',NWMAST,'*I',NSEGM,LSEGM,'1I 2F',IERR)
      endif
      ISSEG = 0


      DO 30 I=1,  ngate_raw

         If (sk_geometry.eq.SK_I) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            DO 35 IQ=1,IQLIST(I)
               CALL SKRN1PE(PE1)
               QLIST(I)=QLIST(I)+PE1
 35         CONTINUE
            PE = QLIST(I)
            

            if((DS_TUNE .eq. LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) then
               thr=0.32
               CALL DSTHR(PE,IFTHR,thr)
            else if(DS_TUNE.eq. AP_TUNE) then
               CALL DSTHR(PE,IFTHR,thr)
            else
               write(*,*) 'DSRES_MULTI: Not Supported DS_TUNE',DS_TUNE
               stop
            endif
            
            IF(IFTHR.EQ.1) GOTO 30
            if(
     & ((DS_TUNE .eq. LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) .and.
     & (DS_TUNE.eq.1.and.iqlist(i).eq.1)) pe=pe*0.935
            QLIST(I) = PE
            
            if((DS_TUNE .eq. LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) then
               call sktsig_le(qlist(i),tsig)
               tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
               CALL rntwogauss(tlist(i),tsig,trnd)
               TLIST(I) = TRND
               
            else if(DS_TUNE.eq.AP_TUNE) then
               qpc = qlist(i) * 2.0555
               call sktsig(qpc,tsig)
c     CALL SKTSIG(QLIST(I),TSIG)
               CALL RNGAUS(0.,TSIG,TRND)
               TLIST(I) = TLIST(I)+TRND
            else
               write(*,*) 'DSRES_MULTI: Not Supported DS_TUNE',DS_TUNE
               stop
            endif

         Else if (sk_geometry.eq.SK_II) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-II
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            DO 36 IQ=1,IQLIST(I)
               CALL SKRN1PE_SK2(PE1)
               QLIST(I)=QLIST(I)+PE1
 36         CONTINUE
            PE = QLIST(I)
            
            CALL DSTHR(PE,IFTHR,thr)
            
            IF(IFTHR.EQ.1) GOTO 30
            QLIST(I) = PE
            
            call sktsig_le(qlist(i),tsig)
            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss(tlist(i),tsig,trnd)
            TLIST(I) = TRND
               
         Else if (sk_geometry.eq.SK_III) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-III
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            icab = iclist(i)

            do IQ=1,IQLIST(I)
               CALL SKRN1PE_SK3(icab,PE1)
               QLIST(I)=QLIST(I)+PE1
            enddo
            PE = QLIST(I)
            
            CALL dsthr_sk3(icab,PE,IFTHR)
            
            IF(IFTHR.EQ.1) GOTO 30
            QLIST(I) = PE

            call sktsig_le_sk3(qlist(i),tsig)
            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss(tlist(i),tsig,trnd)
            TLIST(I) = TRND

         Else if (sk_geometry.eq.SK_IV) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-IV
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            icab = iclist(i)

c            do IQ=1,IQLIST(I)
c               CALL SKRN1PE_SK3(icab,PE1)
c               QLIST(I)=QLIST(I)+PE1
c            enddo
c            PE = QLIST(I)
            PE = charge_list(i)
            
c            CALL dsthr_sk3(icab,PE,IFTHR)
            CALL dsthr_sk4(icab,PE,IFTHR)
            
            IF(IFTHR.EQ.1) GOTO 30
            call dssatu_sk4(pe,pe_satu)  ! saturation effect of SK-IV
            QLIST(I) = PE_satu

            if(nrunsk4_for_gain.eq.0) then ! for ATMPD analysis (not meant gain)
               call sktsig_sk4_asy(qlist(i),tsigE,tsigL)
               tsigE = tsigE * tunefac
               tsigL = tsigL * tunefac
            else                ! for LOWE analysis (not meant gain)
               qbuf = PE_satu/gain_corr(icab)
               call sktsig_sk4_asy(qbuf,tsigE,tsigL)
            endif
c            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss_asy(tlist(i),tsigE,tsigL,trnd)
            TLIST(I) = TRND


         Else if (sk_geometry.eq.SK_V) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-V (just copy from SK-IV J.Kameda 2019/10/29)
C                                    (modified dsthr_sk5 S.Izumiyama 2019/12/27)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            icab = iclist(i)
            PE = charge_list(i)

            CALL dsthr_sk5(icab,PE,IFTHR)

            IF(IFTHR.EQ.1) GOTO 30
            call dssatu_sk4(pe,pe_satu)  ! saturation effect of SK-IV
            QLIST(I) = PE_satu

            if (nrunsk4_for_gain .ne. 0) then
               qbuf = PE_satu/gain_corr(icab)
            else
               qbuf = PE_satu
            endif

C     SK PMT
            if ( (nskpmt(icab) .eq. 2) .or. (nskpmt(icab) .eq. 3)) then
             Call sktsig_sk5_asy(qbuf,tsigE,tsigL)
            elseif (nskpmt(icab) .eq. 5 ) then
C     HK PMT
             Call sktsig_hk_asy(qbuf,tsigE,tsigL)
            else
             write(6,*) 'DSRES_MULTI: PMT category is strange. Stop.'
             stop            
            endif

            CALL rntwogauss_asy(tlist(i),tsigE,tsigL,trnd)
            TLIST(I) = TRND

         Else
            write(6,*) 'supported only sk_geometry 1,2,3,4,5 in dsres'
            stop
         Endif

         IBUF(1) = ICLIST(I)
         RBUF(2) = TLIST(I)
         RBUF(3) = QLIST(I)

         ISEGM = 1
         NDATA = 3
         ISSEG = ISSEG + 1

         If (sk_geometry.le.3)
     &        CALL KZREP2('DATATMP',NSEGM,ISSEG,NDATA,IBUF)

      if (SK_GEOMETRY .ge. SK_IV) then

cc     fill variables in SKTQ bank

         ngate = ngate + 1
         ICABIZ_tmp(ngate) = ICLIST(I)
         QISKZ_tmp(ngate)  = QLIST(I)
c         TISKZ_tmp(ngate)  =  TLIST(I) + real(ntrg+2)*twid_block +trg_offset
c         itiskz_tmp(ngate) = int(tiskz_tmp(ngate)/tdc_unit_sk4)

cc   TISKZ is measured in GEANT time. ITISKZ is measred in Super-K time
cc   to avoid rounding error of TISK.

         TISKZ_tmp(ngate)  =  TLIST(I) 
         itiskz_tmp(ngate) = int((tlist(i)+trg_offset)/tdc_unit_sk4) + (ntrg+2)*2**15 

c  Q digitization.
c  A very simple approximation.
c  More precise calibration results should be put here.

         iq_small  = int(QLIST(I)*pc2pe_mc/cnt2pc_s) + 961
         iq_medium = int(QLIST(I)*pc2pe_mc/cnt2pc_m) + 961
         iq_large  = int(QLIST(I)*pc2pe_mc/cnt2pc_l) + 961

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

c   Saturation of the ADC

c        if (iqiskz_tmp(ngate) .ge. 2**11-1) then
c           iqiskz_tmp(ngate) = 2**11-1
c            qiskz_tmp(ngate) = (2**11-1-961)*cnt2pc_l/pc2pe_mc
c        endif

        IQISKZ_TMP(ngate) = IQISKZ_TMP(ngate) + flag_gain*2**14 ! range flag
 
      endif  

 30   CONTINUE   

      if (SK_GEOMETRY .ge. SK_IV) then

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

      endif


 900  format(I7,I7,I7,I7,F16.5)
C ---
      RETURN
      END
