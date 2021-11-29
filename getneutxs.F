***********************************************************************
*     -----------------
      SUBROUTINE GETNEUTXS(pin,tot,absb,inel,scx,dcx,elas)
*     -----------------
*
*     (Purpose)
*       Load NEUT pion interaction cross sections
*
*     (Input)
*       Input table file
*       pin: pion momentum (MeV/c)
*
*     (Output)
*       Cross sections in mb
*
*     (Creation Date and Author)
*       2010.07.09 ; First version   by  P.de Perio
*
************************************************************************    
     
      implicit none

#include "apcard.h"
      
      character first_line
      character*128 fname
      character*200 filename
      integer   int, i
      
      integer nbins
      parameter (nbins=76)
 
      real tpi,pin
      
      real tpi_tbl(nbins)
      real tot_tbl(nbins),absb_tbl(nbins),inel_tbl(nbins)
      real scx_tbl(nbins),dcx_tbl(nbins),elas_tbl(nbins)
  
      real tot,absb,inel,scx,dcx,elas
      real dummy

      real xint,frac
      
      logical first /.true./
      common /savexs/ tot_tbl,absb_tbl,inel_tbl,scx_tbl,dcx_tbl,elas_tbl
      
      if (IAPIFSIVAR.eq.1) then
         filename='pitbl_qe0.6_ab0.7_cx0.5.dat'
      else if (IAPIFSIVAR.eq.2) then
         filename='pitbl_qe0.6_ab0.7_cx1.6.dat'
      else if (IAPIFSIVAR.eq.3) then
         filename='pitbl_qe0.7_ab1.6_cx0.4.dat'
      else if (IAPIFSIVAR.eq.4) then
         filename='pitbl_qe0.7_ab1.6_cx1.6.dat'
      else if (IAPIFSIVAR.eq.5) then
         filename='pitbl_qe1.4_ab0.6_cx0.6.dat'
      else if (IAPIFSIVAR.eq.6) then
         filename='pitbl_qe1.3_ab0.7_cx1.6.dat'
      else if (IAPIFSIVAR.eq.7) then
         filename='pitbl_qe1.5_ab1.5_cx0.4.dat'
      else if (IAPIFSIVAR.eq.8) then
         filename='pitbl_qe1.6_ab1.6_cx1.6.dat'
      else
         filename='pitbl_qe1.0_ab1.1_cx1.0.dat'
      end if

      if (first) then
         call findconsts(fname,'NEUT Pion Hadronic xsec Table',filename,0)
         open (unit=10,file=fname,status='old')
         read (10,'(a)') first_line
         do i = 1, nbins
            read (10,*) tpi_tbl(i),tot_tbl(i),absb_tbl(i),inel_tbl(i),
     &           scx_tbl(i),dcx_tbl(i),elas_tbl(i)
         end do
         close (10)

C     Write to log file
         write(*,*)
         write(*,*) 'T(MeV) Total(mb)  Absorb(mb)    QE(mb)      SCX(mb)     DCX (mb)    Elastic(mb)'
         
         do i=1,nbins
            write(*,*) tpi_tbl(i),tot_tbl(i),absb_tbl(i),inel_tbl(i),
     &           scx_tbl(i),dcx_tbl(i),elas_tbl(i)
         enddo
         
         first = .false.
      end if

      tpi=sqrt(139.57**2+pin**2)-139.57

      if (tpi.gt.380.) then
         write (6,*) 'GETNEUTXS: Error NEUT cross sections not valid for 
     &                tpi (',tpi,') > 380'
         stop
      else
         xint = tpi/10. + 1.
         int  = ifix(xint)
         frac = (tpi - 10.*float(int-1))/10.

         tot = tot_tbl(int) + frac*( tot_tbl(int+1) - tot_tbl(int) )
         absb = absb_tbl(int) + frac*( absb_tbl(int+1) - absb_tbl(int) )
         inel = inel_tbl(int) + frac*( inel_tbl(int+1) - inel_tbl(int) )
         scx = scx_tbl(int) + frac*( scx_tbl(int+1) - scx_tbl(int) )
         dcx = dcx_tbl(int) + frac*( dcx_tbl(int+1) - dcx_tbl(int) )
         elas = elas_tbl(int) + frac*( elas_tbl(int+1) - elas_tbl(int) )
      end if
      
      return
      end
