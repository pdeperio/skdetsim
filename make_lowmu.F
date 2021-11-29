      program skdetsim
CC  Main routine for detector simulation 
CC  for Super-Kamiokande.
CC  Merged from 02a/apdetsim.F and ledetsim_4.0/ledetsim.F
CC  by J.Kameda
CC

      implicit none

#include "dsbnkf.h"
#include "gcbank.h"
#include "skparm.h"
c#include "sk1pe.h"
#include "apcard.h"
#include "aptrkseg.h"
#include "skdetsim_parm.h"
#include "skdetsim_test.h"
#include "skhead.h"
#include "guxyzpm.h"

      Integer iargc
      Integer ierr1, ierr2
      Integer icskip
      character*132 fskip
      integer ii
C
C+  Declare GEANT user routines external, so the linker finds them in our
C+  library
C

      external gukine,gutrak,gutrev,gustep,guphad,guhadr,grun,gtelec,
     &         gtckov,ggckov,gtmuon,gpart,gcalor,gnsphr,guout,skrn1pe

      Data lin /10/
      Data lout/20/

#ifdef LINAC
#include "linpara.h"
#include "lininfo.h"
#include "lininfo2.h"
#include "skday.h"
#endif

#ifdef LINAC
      integer run_num, i, irun
      real wt,watrds,dark
      character*1 dummy
      character*80 argv

CPC      call skrunday
CPC      call skwt

      call getarg(3,argv)
      read(argv,*) run_num

      call linsum2(run_num)

c      call getarg(2,argv)
c      read(argv,*) pinit

***  Dark rate should be set each run..
      lin_dark = -1
      if(lin_badrun(run_num).eq.-1) then
         write(6,*) 'Dark rate cannot be set in LINAC run...'
         write(6,*) 'Check Run number and'
         if(run_num.lt.30000) then
            write(6,*) '/home/sklowe/linac/const/linac_sk2_runsum.dat'
         else
            write(6,*) '/home/sklowe/linac/const/linac_sk3_runsum.dat'
         endif
         stop
      else

*  try dark rate from badsel
         open(98,file='/skam/const/runsum',
     &        status='old',form='formatted')
         do i=1,7
            read(98,'(a)')dummy
         enddo

         irun=0
         do while(irun.lt.lin_badrun(run_num))
            read(98,'(i6,57x,f9.5)',end=3)irun,dark
         enddo

 3       if(irun.ne.lin_badrun(run_num)) then
            write(6,*) 'cannot find dark rate from the table...'
            write(6,*) 'Check Run number and'
            if(run_num.lt.30000) then
               write(6,*) '/home/sklowe/linac/const/linac_sk2_runsum.dat'
            else
               write(6,*) '/home/sklowe/linac/const/linac_sk3_runsum.dat'
            endif
            close(98)
            stop
         endif
         lin_dark = dark
      endif

      pinit = lin_mom(run_num)
      if (pinit.le.0) then
         write(*,*) 'Input momentum is unknown !!'
         write(*,*) 'I use default value of',lin_mode(run_num),
     &        'MeV mode.'
         if ( lin_mode(run_num) .eq. 15) then
            pinit = 16.294
         elseif ( lin_mode(run_num) .eq. 12) then
            pinit = 13.644
         elseif ( lin_mode(run_num) .eq. 10) then
            pinit = 10.982
         elseif ( lin_mode(run_num) .eq. 8) then
            pinit = 8.861
         elseif ( lin_mode(run_num) .eq. 6) then
            pinit = 6.989
         elseif ( lin_mode(run_num) .eq. 5) then
            pinit = 6.034
         elseif ( lin_mode(run_num) .eq. 4) then
            pinit = 5.078
         else
            write(*,*) lin_mode(run_num),'MeV mode is unknown!'
            stop
         endif
      endif

CPC      call lfwater(relapse(run_num),wt)
CPC      WATRDS = wt
      WATRDS = 8800.

      sorx = lin_x(run_num)
      sory = lin_y(run_num)
      sorz = lin_z(run_num)

      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) 'RUN#',run_num,' is simulated.'
      write(*,*) 'Position =',sorx,sory,sorz
      write(*,*) 'Input momentum =',pinit
      write(*,*) 'Water transparency is',WATRDS
      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

#endif      

      Open(89,file='./lowmuvector.dat',status='old'
     &     ,form='formatted')

C     Initialization of ZEBRA system, dynamic core allocation

      Call gzebra(KWBANK)

C     User code for initialization

C --- get card file name (1st argument of the program)

      If (iargc().ge.1) then
         Call getarg(1, fcard)
      Else
        fcard = 'supersim.card'
      Endif
      Write(6,*) 'The card file you read is ',fcard
C
C+ Initialize GEANT
C
      Call sginit

C   
C   Skip events
C
      if (iargc().ge.2) then
	call getarg(2,fskip)
	read(fskip,'(I10)') icskip
	write(6,*) 'User requested skipping ',icskip,' events on command line'
	iskipev = max(iskipev,icskip)
      endif
C
C+ Open input data file if requested
C
      If (iflgv.eq.0.or.iflgv.eq.-2)
     &     Call skopenf(lin,1,'Z',ierr1)
      if (iflgv.eq.2) call skopenf(lin,1,'F',ierr1)
      if (ierr1.ne.0) write(6,*)'APDETSIM - error opening input file'
      if (OUTPUT_HBOOK .eq. OUT_SEG) then
         call seg_openf
      endif
C
C+ Open output file
C
      if  ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON)
     &     .and.(OUTPUT_ZBS)) Call skopenf(lout,1,'Z',ierr2)

C
C+ Prepare hbook for standard diagnostic test histograms
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         call sktest_init
      endif

C
C+ Mode -2 means read sub-event text file
C
      If (iflgv.eq.-2) 
     &     Open(99,file='sube.dat',form='formatted',status='unknown')
C
C+ Generate data
C
      Call grun
C
C+ Close input and output files
C
      If (iflgv.eq.0.or.iflgv.eq.-2) Call skclosef(lin)
      if  ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.
     &     (OUTPUT_ZBS)) CALL skclosef(lout)

      If (iflgv.eq.-2) close(99)
      If (OUTPUT_HBOOK.eq.OUT_SEG) call seg_closef
C
C+ Close Standard HBOOK diagnostic file
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         call sktest_close
      endif

C
C+  Standard GEANT3 termination
C
      Call glast
C
      Close(89)

C
      End


******************************************************
*     -----------------
      SUBROUTINE SGGVUS
*     -----------------
*
*     (Purpose)
*        User routine for Generating Vector
*
******************************************************
#include "vcwork.h"
#include "geotnk.h"
#include "apcard.h"

      REAL pos(3),pmom(3),dir(3)
      Real pp
      Integer i
      Real dummy
      Integer repeat_max, nrepeat
      Parameter (repeat_max=8)
      Data nrepeat/0/
      Real poffset
      Parameter (poffset=30.)
      Integer icall
      Data icall/0/
      
      If (icall.eq.0) then
         print *, 'skip ',(iskipev-1)*125,' events.'
         Do i=1,(iskipev-1)*125
            Read(89,*) pos(1),pos(2),pos(3)
     &           ,pmom(1),pmom(2),pmom(3),pp
         Enddo
         icall = 1
      Endif

      nrepeat=nrepeat+1
      if (nrepeat.gt.repeat_max) then
         nrepeat=1
      endif
      If (nrepeat.eq.1) then
         Read(89,*) pos(1),pos(2),pos(3),dir(1),dir(2),dir(3),pp
      endif

      pos(3) = zpintk-1.
C      PMOM(1)=(DIR(1) * (PP+poffset))*(1.0+real(nrepeat-25)/200.)
C      PMOM(2)=(DIR(2) * (PP+poffset))*(1.0+real(nrepeat-25)/200.)
C      PMOM(3)=(DIR(3) * (PP+poffset))*(1.0+real(nrepeat-25)/200.)
      PMOM(1)=(DIR(1) * PP)
      PMOM(2)=(DIR(2) * PP)
      PMOM(3)=(DIR(3) * PP)
      IPART =13
      print *, 'pos,mom=',pos(1), pos(2), pos(3)
     &     ,pmom(1), pmom(2), pmom(3)


      CALL VCMAKE(POS,PMOM,IPART)

      RETURN
      END


