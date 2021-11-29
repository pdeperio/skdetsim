************************************************************************
*     -------------------------------------
      SUBROUTINE DSTHR(PE,IFLAG,thr_le)
*     -------------------------------------
*
*     (Purpose)
*        simulate threthold effect
*
*     (Creation Date and Author)
*       1996.??.?? ; First version by K.Ishihara
*       1998.03.24 ; Modified by Ishihara 
*                    probabrity get by (nomal 1pe)/(PMT'gain higher 1pe)
*       1998.07.30 ; modified by D.Casper
*                    check x value before calculating threshold to
*                    avoid floating overflow
*       2001.02.09 ; RNDM -> GRNDM
*       2003.01.07 ; Merge le/ap by Y.Koshio
*       2004.04.17 ; New one for SK-II by M.Shiozawa
*
************************************************************************

      implicit none

#include "skdetsim_parm.h"
#include "skhead.h"

      integer iflag
      real pe, thr_le, x, thr, rdummy
      real thr_err, ped_err, th, rnpe, err
      parameter(thr_err=0.03, ped_err=0.3)
      real pc2pe_sk2

      if (sk_geometry.eq.1) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C     SK-I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
         if((DS_TUNE.eq.LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) then
            call rngausk2(thr_le, THR_ERR, th)
            iflag = 0
            if (pe .lt. th) then
               pe = 0.0
               iflag = 1
            else
               call rngausk2(pe, PED_ERR, rnpe)
               pe = rnpe
            endif
            
         else if(DS_TUNE.eq.AP_TUNE) then
            iflag = 0
            x=pe+0.10
            
            if (x.lt.1.1) then
               thr = min(1.0,
     $              -0.06374+x*(3.748+x*(-63.23+x*(452.0+x*(-1449.0+x*(2513.0
     $              +x*(-2529.+x*(1472.0+x*(-452.2+x*(51.34+x*2.370))))))))))
            else
               thr = 1.0
            endif
C     if (thr.lt.rndm(DUM)) then
            CALL GRNDM(RDUMMY,1)
            if (thr.lt.RDUMMY) then
               pe = 0.0
               iflag = 1
               return
            else
               call rngaus(0.0, 0.03, err)
               pe = pe+err
            endif
         else
           write(*,*) 'DSTHR: Not Supported DS_TUNE',DS_TUNE
           stop
         endif
      else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C     SK-II, common for lowe and highe analyses
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
         iflag = 0
         th=0.45                !pC
         pc2pe_sk2=2.297
         if (pe.lt.th/pc2pe_sk2) then
            pe=0.
            iflag = 1
         else
            call rngaus(0.0, 0.06, err)
            pe = pe+err
         endif
      endif


      return
      end
