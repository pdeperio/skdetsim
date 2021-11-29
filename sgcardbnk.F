************************************************************************
*     -----------------
      SUBROUTINE SGCARDBNK
*     -----------------
*
*     (Purpose)
*       store card file into ZBS as text
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30 ; Koshio
*
************************************************************************
      implicit none
      
#include "dsbnkf.h"
      integer nfflines
      character*80 chfftexts(200),chdum
      integer ifftexts(20*200)
      equivalence (chfftexts(1),ifftexts(1))
      save chfftexts
      integer ifirst /0/
      save ifirst
      integer ierr,nseg,i

C -- extract VECT card

      if(ifirst.eq.0)then
         ifirst = 1
         OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
 100     read(90,'(a)',end=900,err=900)chdum
         nfflines = nfflines + 1
         chfftexts(nfflines) = chdum
         goto 100
 900     continue
         CLOSE(UNIT=90)
         write(6,*)'sgcardbnk : nfflines =',nfflines
         do i=1,nfflines
            write(6,*)chfftexts(i)
         enddo
      endif

      call kznseg('APDETCAR',nseg)
      if(nseg.ge.0)call kzbdel('APDETCAR')
      call kzbcr1('APDETCAR',nfflines*20,'-I',0,0,' ',ierr)
      call kzrep1('APDETCAR',0,nfflines*20,ifftexts)
      
      RETURN
      END


