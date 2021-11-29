************************************************************************
*     -------------------------------------
      SUBROUTINE dsthr_sk4(icab,PE, IFLAG)
*     -------------------------------------
*
*     (Purpose)
*        simulate QBee threthold effect for SK4 case
*        copy from dsthr.
*        For SK4, QBee resolution is also considered here.     
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
*       2009.01.18 ; First version by T.Tanaka
*
************************************************************************

      implicit none

      integer icab, iflag
      real pe,x,thr,pc,pc2pe,rpc
      parameter (pc2pe=2.243) ! sk3 pC2pe(used for calibration)       

      real th_paramsk2,th_paramsk3
      parameter (th_paramsk2=1.74909)
      parameter (th_paramsk3=1.56925)   
 
      integer maxpm, lun
      parameter (maxpm=11146, lun=99)
      integer n0, n1, n2,isk23(maxpm),isave/0/
      save isave, isk23

      character*132 file_name
      data file_name/'/skam/const/pmtinf.dat'/


      if(isave.eq.0) then
         isave = 1
         call findconsts(file_name,'pmtinf','pmtinf.dat',0)
         open(LUN,file=file_name,status='old')
 3       continue
         read(LUN,*,END=4) n0, n1, n2
         isk23(n0) = n1
         goto 3
 4       continue
         close(LUN)
      endif



CCC   convert pe -> pC temporaly(using sk3 pc2pe)
      pc = pe*pc2pe

CCC   sk4 threshold
CCC   using function from electrics calibration
CCC   f(x) = 0.5 * ( Erf(17.09*(a*x-0.649))+1 )
CCC   a= free parameter(tuned by analysis)

         if(isk23(icab).eq.2) then
         thr = 0.5 * ( Erf(17.09*(th_paramsk2*pc-0.649))+1 )            
         else if(isk23(icab).eq.3) then
         thr = 0.5 * ( Erf(17.09*(th_paramsk3*pc-0.649))+1 )
         else
         thr = 1
         endif

      call ranlux(x,1)

      iflag = 0
      if (x.ge.thr) then
         rpc=0.
         iflag = 1
      endif
       
CCC  consider resolution effect near threshold
       if (iflag.eq.0) then !! for hit case
       call dsreso_qb(icab,pc,rpc)
       endif

CCC  resume unit p.e. for SK4(using sk4 pc2pe)
       pe = rpc/pc2pe

      return
      end
