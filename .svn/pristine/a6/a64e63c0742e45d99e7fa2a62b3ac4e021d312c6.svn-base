************************************************************************
*     -------------------------------------
      SUBROUTINE dsthr_sk5(icab,PE, IFLAG)
*     -------------------------------------
*
*     (Purpose)
*        simulate QBee threthold effect for SK5 case
*        copy from dsthr.
*        For SK5, QBee resolution is also considered here and revised
*        from dsthr_sk4.f.
* 
*
*     (INPUT)
*        icab  : cable number
*        pe    : photo electron value
*            
*      (OUTPUT)
*        iflag : above threshold or not
*                0 means hit, and 1 means non-hit
*        pe    : photo electron value
*            
*     (Creation Date and Author)
*       2019.12.05 : First version by S.Izumiyama
*         Modified threshold function of SK4 in order to match SK5 data
*       2021.05.20 : J.Kameda
*         Modify to change the channel identifier to nskpmt(i).    
************************************************************************

      implicit none

#include "skparm.h"            
#include "skpmt.h"
      
      integer icab, iflag
      real pe,x,thr,pc,pc2pe
c     parameter (pc2pe=2.243) ! sk3 pC2pe(used for calibration)       
      parameter (pc2pe=2.463) ! sk5 global pC2pe       

      real thrSlope_sk2, thrSlope_sk3
      real thrShift_sk2, thrShift_sk3
      parameter (thrSlope_sk2=3.915)
      parameter (thrSlope_sk3=3.307)   
      parameter (thrShift_sk2=0.4771)
      parameter (thrShift_sk3=0.5657)   
 
CCC   convert pe -> pC temporaly(using sk5 pc2pe)
      pc = pe*pc2pe

CCC   sk5 threshold
CCC   using function from electrics calibration
CCC   f(x) = 0.5 * ( Erf(a*(x-b))+1 )
CCC   a,b= free parameter(tuned by analysis)

         if( nskpmt(icab) .eq. 2) then
            thr = 0.5*(Erf(thrSlope_sk2*(pc-thrShift_sk2))+1)
         else if(nskpmt(icab).eq.3) then
            thr = 0.5*(Erf(thrSlope_sk3*(pc-thrShift_sk3))+1)
         else 
            thr = 1
         endif
         
      call ranlux(x,1)

      iflag = 0
      if (x.ge.thr) then
         iflag = 1
      endif
       

CCC  resume unit p.e. for SK5(using sk5 pc2pe)
       pe = pc/pc2pe

      return
      end
