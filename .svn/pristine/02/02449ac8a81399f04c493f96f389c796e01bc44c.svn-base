C Output information about the event to be used to check if new
C code has changed the behavior of skdetsim.
C 
C Currently, this means dump the ADC/TDC values of all hit tubes
C to a text file, so we can use the power of "diff"
C
C If you've got something cleverer, add it here.
C
C The VALIDATE logical card controls the use of this feature.
C
C See validate.sh for a script that mostly automates things.

      subroutine do_validate
      implicit none

#include "skdetsim_parm.h"
#include "skparm.h"
#include "sktq.h"
#include "skhead.h"


      integer ierr
      integer i
      logical first /.true./

      integer ICHK(MAXPM,2)
      COMMON /SKRDLC/ICHK



      if (first) then        ! Initialize things
         open(file='validation.dat',unit=69,form='formatted',
     &        status='new')
         first = .false.
      endif

cc --- initialize inner common

      ierr = 0

      NQISK = 0
      QISMSK = 0.
      QIMXSK = 0.
      TIMNSK = 100000.
      TIMXSK = -100000.
      
      call lbfset0(ICHK  ,MAXPM*2)
      call lbfset0(QISK  ,MAXPM  )
      call lbfset0(TISK  ,MAXPM  )
      call lbfset0(IQABSK,MAXPM*2)
      call lbfset0(ITABSK,MAXPM*2)

cc --- initialize anti common      

      nhitaz = 0
      nqask = 0
      qasmsk = 0.
      qamxsk = 0.
      tamnsk = 100000.
      tamxsk = -100000.
      call lbfset0(QASK  ,MAXPMA )
      call lbfset0(TASK  ,MAXPMA )
      call lbfset0(IHTFLZ,MAXPMA*16)
      call lbfset0(ICABAZ,MAXPMA*16)
      call lbfset0(TASKZ ,MAXPMA*16)
      call lbfset0(QASKZ ,MAXPMA*16)
      
C Fill commons
      call headsk
c      call aheadsk
      call tqsk
      call cmtqsk
      call tqask

C Loop over TQ, TQA banks, writing out T & Q for each PMT
      write (69,100) nevsk

      if (SK_GEOMETRY.ge.SK_IV) then
         
         call tqrealsk
         call tqarealsk(ierr)
         
C     First do the TQ bank
         write (69,101) 'ID',nqiskz,qismsk
         do i=1,nqiskz
            write(69,102) icabiz(i),qiskz(i),tiskz(i)
         end do
         
C     Next do the TQA bank
         write (69,101) 'OD',nhitaz,qasmsk
         do i=1,nhitaz
            write(69,102) icabaz(i),qaskz(i),taskz(i)
         end do
         write (69,*)
         
      else
                  
C     First do the TQ bank
         write (69,101) 'ID',nqisk,qismsk
         do i=1,nqisk
            write(69,102) ihcab(i),qisk(ihcab(i)),tisk(ihcab(i))
         end do
         
C     Next do the TQA bank
         write (69,101) 'OD',nhitaz,qasmsk
         do i=1,nhitaz
            write(69,102) icabaz(i),qaskz(i),taskz(i)
         end do
         write (69,*)
         
      endif


 100  format(1x,'Event # ',i8)
 101  format(1x,a2,' hit list has ',I6,' hit tubes and ',f10.1,
     &       ' total pe')
 102  format(1x,i6,f8.3,f9.2)

      return
      end
