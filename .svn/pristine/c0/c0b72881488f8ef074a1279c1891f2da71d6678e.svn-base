      subroutine seg_openf
      implicit none
#include "aptrkseg.h"

      integer lrec, istat, lnblnk, ierr
C
      if (seg_filebuff(1).ne.0) then
        call hlimit(-nwpawc)
        iquest(10) = 1000000
        lrec = 8192
        print *,'Opening hbook file:',seg_filename(1:lname),lnblnk(seg_filename),lname
        call hropen(lun_seg, 'segments',seg_filename(1:lname),'NQX',lrec, istat)

        call hbnt(seg_idh,'Track segments',' ')
        call hbset('BSIZE',lrec,ierr)
        call hbname(seg_idh,'segments',nseg,'nseg[0,12288]:i,'//
     &     'ntot:i, time(nseg):r, theta(nseg):r, phi(nseg):r, nph(nseg):i')
        call hbname(seg_idh,'tracks',npart,'npart[0,1000]:i,'//
     &     'nevt:i, ndk:i, subevt:i, ip(npart):i, p(5,npart):r, thp(npart):r, php(npart):r, cls(npart)[0,7]:i')
      endif

      return
      end
