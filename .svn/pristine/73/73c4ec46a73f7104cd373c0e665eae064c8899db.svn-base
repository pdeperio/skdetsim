************************************************************************
*     -------------------------------------
      SUBROUTINE dsthr_sk3(icab, PE, IFLAG)
*     -------------------------------------
*
*     (Purpose)
*        simulate ATM threthold effect for SK3 case
*        copy from dsthr.
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
*       2007.01.22 ; First version by Y.Koshio
*
************************************************************************

      implicit none

      integer icab, iflag
      real pe, x, a2(9), a3(9), thr
      data a2/0.12199E-02,-0.54256E-01,6.7221,-282.26,7063.5,-99031.,
     &     0.72999E+06,-0.24579E+07,0.30070E+07/
      data a3/0.53046E-03,0.82390E-01,-1.3782,-32.843,2517.9,-46947.,
     &     0.37387E+06,-0.12522E+07,0.14798E+07/

      integer maxpm, lun
      parameter (maxpm=11146, lun=99)
      integer n0, n1, n2, isk23(maxpm), isave/0/
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

      if(isk23(icab).eq.2) then
         if (pe.lt.0.275) then
            thr = min(1.0,
     $           a2(1)+pe*(a2(2)+pe*(a2(3)+pe*(a2(4)+pe*(a2(5)+pe*(a2(6)
     $        +pe*(a2(7)+pe*(a2(8)+pe*a2(9)))))))))
         else
            thr = 1.0
         endif

      else if(isk23(icab).eq.3) then
         if (pe.lt.0.285) then
            thr = min(1.0,
     $           a3(1)+pe*(a3(2)+pe*(a3(3)+pe*(a3(4)+pe*(a3(5)+pe*(a3(6)
     $           +pe*(a3(7)+pe*(a3(8)+pe*a3(9)))))))))
         else
            thr = 1.0
         endif

      else
         write(0,*) 'warning in dsthr_sk3 missing cable =',icab
         pe = 0.
         iflag = 1
         return
      endif

      call ranlux(x,1)

      iflag = 0
      if (x.ge.thr) then
         pe=0.
         iflag = 1
      endif

      return
      end
