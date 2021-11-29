      program skdetsim
CC  Main routine for detector simulation 
CC  for Super-Kamiokande.
CC  Merged from 02a/apdetsim.F and ledetsim_4.0/ledetsim.F
CC  by J.Kameda
CC
cc  09-JAN-2009 modified by y.takeuchi
cc    o remove RFLIST
cc    o added skroot/zbs for OUTPUT_TYPE=1
cc

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
      character*200 argv, fname_in, fname_out, fname

C
C+  Declare GEANT user routines external, so the linker finds them in our
C+  library
C

      external gukine,gutrak,gutrev,gustep,guphad,guhadr,grun,gtelec,
     &         gtckov,ggckov,gtmuon,gpart,gcalor,gnsphr,guout,skrn1pe
#ifdef LASER
     &       , gtnino, lsrgen
#endif

      Data lin /10/
      Data lout/20/

#ifdef NICKEL
#include "nicpara.h"
      character*132 fpos
#endif

#ifdef LINAC
#include "linpara.h"
#include "lininfo.h"
#include "lininfo2.h"
#include "skday.h"
#include "dsparm.h"

      integer run_num, i, irun
      real wt,dark
      character*1 dummy
#endif

#ifdef DTG
#include "dtginfo.h"
      integer dt_run_num
#endif

#ifdef WTSCTT
#include "guwtr.h"
#include "rayfac.h"
#include "wtsctt.h"
#include "dsparm.h"
#endif

#ifdef LASER
#include "laser.h"
#endif


c*** init. variables
      icskip = 0


c*** check & read arguments
#ifdef LINAC
      if (iargc().lt.3) then
         print *, "LINAC: skdetsim card output_file run_number"
         stop
      endif
      call getarg(3,argv)
      read(argv,*) run_num
#else
#ifdef NICKEL
      if (iargc().lt.5) then
         print *, "NI   : skdetsim card output_file x y z"
         stop
      endif
      call getarg(3,fpos)
      read(fpos,*) sorx
      call getarg(4,fpos)
      read(fpos,*) sory
      call getarg(5,fpos)
      read(fpos,*) sorz
#else
#ifdef DTG
      if(iargc().lt.3)then
         write(*,*) ' '
         write(*,*)'skdetsim <card name> <output_file> <dtg_runnum>'
         write(*,*) ' '
         stop
      end if

      call getarg(3,argv)
      read(argv,*) dt_run_num

      call dtgsum
      dtg_runnum = dt_run_num

      write(*,*) 'dt_run_num:',dt_run_num ! S.Ito
      write(*,*) 'dtg_runnum:',dtg_runnum ! S.Ito

      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) 'Calling Pope/Fry Attenuation readfile'
      call sgpope_read   ! To get Pope/Fry Attenuation Coefficients
      write(*,*) 'Read successful!'
      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) "DTG Monte Carlo for SK-IV"
      write(*,*) ' '
      write(*,*) 'RUN#',dtg_runnum,' is simulated.'
      write(*,*) 'Position x =',dtg_x(dtg_runnum)
      write(*,*) '         y =',dtg_y(dtg_runnum)
      write(*,*) '         z =',dtg_z(dtg_runnum)
c      write(*,*) 'Water transparency is',WATRDS
      write(*,*) 'Bad channel run:',dtg_badrun(dtg_runnum)
c -- Open vector file with MC inputs

c      fname = '/home/guillian/work/dtgen/mc1/ofl/n16_input_eric.txt'
c      fname = '/home/sekiya/skrep/n16_input_eric.txt'
      fname = '/home/sito/skrep/n16_input_eric.txt'
      open(4, FILE = fname)
      write(6,*) 'Choice of n16 input:'
      write(6,*) '  --> ', fname
#else
#ifndef LASER
#ifndef WTSCTT
      if (iargc().lt.2) then
         print *, "usage: skdetsim card output_file",
     &        "[input_file] [iskip]"
         stop
      endif
      if (iargc().ge.3) call getarg(3, fname_in)
      if (iargc().ge.4) then
         call getarg(4, fskip)
         read(fskip,*) icskip
      endif
#endif
#endif
#endif
#endif
#endif


c*** read common arguments
C --- get card file name (1st argument of the program)
      call getarg(1, fcard)
      Write(6,*) 'The card file you read is ',fcard
      call getarg(2, fname_out)

      

#ifdef LINAC
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
         else if(run_num.lt.40000) then
            write(6,*) '/home/sklowe/linac/const/linac_sk3_runsum.dat'
         else
            write(6,*) '/home/sklowe/linac/const/linac_sk4_runsum.dat'
         endif
         stop
      else

*  try dark rate from badsel
ccyt         open(98,file='/skam/const/runsum',status='old',form='formatted')
         call findconsts(fname,'runsum','runsum',2)
         open(98,file=fname,status='old',form='formatted')
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
            else if(run_num.lt.40000) then
               write(6,*) '/home/sklowe/linac/const/linac_sk3_runsum.dat'
            else
               write(6,*) '/home/sklowe/linac/const/linac_sk4_runsum.dat'
            endif
            close(98)
            stop
         endif
         lin_dark = dark
      endif
      write(6,*) '*************************************'
      write(6,*) 'Dark rate is set to ',lin_dark,' kHz'
      write(6,*) '   from Run #',irun,' ',lin_badrun(run_num)
      write(6,*) '*************************************'

      pinit = lin_mom(run_num)
      if (pinit.le.0) then
         write(*,*) 'Input momentum is unknown !!'
         write(*,*) 'I use default value of',
     &        lin_mode(run_num),'MeV mode.'
         if ( lin_mode(run_num) .eq. 18) then
            pinit = 18.938
         elseif ( lin_mode(run_num) .eq. 15) then
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
         elseif ( lin_mode(run_num) .eq. 3) then
            pinit = 4.700
         else
            write(*,*) lin_mode(run_num),'MeV mode is unknown!'
            stop
         endif
      endif

      call skrunday
      if(run_num.ge.60000) then ! for SK-IV
         call skwt_gain_corr
      else
         call skwt
      endif

      call lfwater(relapse(run_num),wt)
      WATRDS = wt
c      WATRDS = 13954.8
c      WATRDS = 12414.166          ! Aug. 2010 LINAC fes

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

c      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c      write(*,*) 'Calling Pope/Fry Attenuation readfile (lowe only)'
c      call readfile   ! To get Pope/Fry Attenuation Coefficients
c      write(*,*) 'Read successful!'
c      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

C     Initialization of ZEBRA system, dynamic core allocation

      Call gzebra(KWBANK)

C     User code for initialization

C
C+ Initialize GEANT
C
      Call sginit

C   
C   Skip events
C
      if (icskip.gt.0) then
	write(6,*) 'User requested skipping ',icskip,' events on command line'
	iskipev = max(iskipev,icskip)
      endif

#ifdef LASER
      if (iargc().ge.5) then
         call getarg(3,argv)
         read(argv,*)nphoto
         call getarg(4,argv)
         read(argv,*)blamda
         call getarg(5,argv)
         read(argv,*)sigbm
         sigbm = sigbm*3.14159265/180.
         write(*,*)"***********************************************"
         print *,'check nphoto,blamda,sigbm',nphoto,blamda,sigbm
         write(*,*)'***********************************************'
      else
         write(0,*) 'Error! A few arguments for Laser Simulation!'
         write(0,*) 'Usage: skdetsim', 
     &        '[card] [fskip] [nphoto] [blamda] [sigbm]'
         stop
      endif
#endif
 
#ifdef WTSCTT
CMM.. Set water parameters from comandline arguments

#ifdef LASER
      if (iargc().ge.8) then
         call getarg(6,argv)
         read(argv,*)yooabs
         call getarg(7,argv)
         read(argv,*)yooray
         call getarg(8,argv)
         read(argv,*)yoomie
#else
      if (iargc().ge.5) then
         call getarg(3,argv)
         read(argv,*)yooabs
         call getarg(4,argv)
         read(argv,*)yooray
         call getarg(5,argv)
         read(argv,*)yoomie
#endif
      else
         yooabs = -1.0
         yooray = -1.0
         yoomie = -1.0
      endif
      IF (WATER_TUNE.EQ.AP_TUNE) THEN
         print*, '*****************************************************'
         print*, 'Default parameters of abs, sym, asy : ', AABL,ARAS,AMIS
         print*, 'Factors of abs, sym, asy : ', yooabs, yooray, yoomie
         IF(yooabs.ge.0) AABL = AABL*yooabs
         IF(yooray.ge.0) ARAS = ARAS*yooray
         IF(yoomie.ge.0) AMIS = AMIS*yoomie
         print*, 'Set Parameters of abs, sym, asy : ', AABL,ARAS,AMIS
         print*, '*****************************************************'
      ELSE
         ablfac = 1.0
         IF(IMWATR.EQ.0) THEN
            print*, '**************************************************'
            print*, 'Default parameters of abs, abl, sym, asy : ',
     &           absfac, ablfac, rayfac, miefac
            print*, 'Factors of abs, sym, asy : ', yooabs, yooray, yoomie
            IF(yooabs.ge.0) then
               absfac=absfac*yooabs
               ablfac=yooabs
            ENDIF
            IF(yooray.ge.0) rayfac=rayfac*yooray
            IF(yoomie.ge.0) miefac=miefac*yoomie
         ENDIF
         print*, 'Set Parameters of abs, abl, sym, asy : ', absfac,
     &        ablfac, rayfac, miefac
         print*, '****************************************************'
      ENDIF
#endif

#ifdef LINAC
      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) 'RUN#',run_num,' is simulated.'
      write(*,*) 'Position =',sorx,sory,sorz
      write(*,*) 'Input momentum =',pinit
      write(*,*) 'Water transparency is',WATRDS
      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#endif

#ifdef DSSWTRIGSTOREINFO
      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) 'Recording the software event trigger.'
      write(*,*) 'Bank name: DSTRGOFF'
      write(*,*) 'Segments: 11'
      write(*,*) ' seg 0 (MASTER SEGMENT): 5I'
      write(*,*) '   (I) length of master: 5'
      write(*,*) '   (I) length of seg1 to seg10: 4'
      write(*,*) '   (I) IT0_OFFSET of event: random clock position'
      write(*,*) '   (I) geant0 in clock counts'
      write(*,*) '   (I) number of software triggers'
      write(*,*) ' seg 1-10: Stores info from soft trigger: 2I1F1I'
      write(*,*) '   note: seg1 is the primary trigger'
      write(*,*) '   (I) sw trig id:(-1=NONE),(Bit flag:1=LE,2=HE,3=OD)'
      write(*,*) '   (I) it0sk of trigger'
      write(*,*) '   (F) geant_t0 relative to software trigger'
      write(*,*) '   (I) trigger candidate number'
      write(*,*) 'Bank info filled in dsswtrg.F'
#endif


C
C+ Open input data file if requested
C
      If (iflgv.eq.0.or.iflgv.eq.-2) then
c*** ZBS file
c         if (fname_in(1:5) .eq. '/disk') then
c            call set_rflist(lin,fname_in,'DISK',' ','RED',' ',' ',
c     &           'recl=5670 status=old',' ',' ')
c         else
            call set_rflist(lin,fname_in,'LOCAL',' ','RED',' ',' ',
     &           'recl=5670 status=old',' ',' ')
c         endif
         Call skopenf(lin,1,'Z',ierr1)
      endif

      if (iflgv.eq.2) then
         call set_rflist(lin, fname_in, 'LOCAL',' ','RED',' ',' ',
     $        'recl=5670 form=formatted',' ',' ')
         call skopenf(lin,1,'F',ierr1)
       endif

      if (ierr1.ne.0) write(6,*)'APDETSIM - error opening input file'
      if (OUTPUT_HBOOK .eq. OUT_SEG) then
         call seg_openf
      endif
C
C+ Open output file
C
      if  ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.(OUTPUT_ZBS)) then

         if (OUTPUT_TYPE.eq.2) then ! atmpd mode
c***  ZBS file
c            if (fname_out(1:5) .eq. '/disk') then
c               call set_rflist(lout, fname_out,'DISK',' ','RPL',' ',' ',
c     $              'recl=5670 status=new',' ',' ')
c            else
               call set_rflist(lout, fname_out,'LOCAL',' ','RPL',' ',' ',
     $              'recl=5670 status=new',' ',' ')
c            endif
            call skopenf(lout,1,'Z',ierr2)
         else ! le mode 
            if (SK_FILE_FORMAT .eq. 1) then
               print *,'----- skroot: ', fname_out(1:lnblnk(fname_out))
               call skroot_open_write(lout,fname_out(1:lnblnk(fname_out)))
               call skoptn('31,30')
            else
               print *,'----- zbs: ', fname_out(1:lnblnk(fname_out))
c               if (fname_out(1:5) .eq. '/disk') then
c                  call set_rflist(lout, fname_out, 'DISK',' ','RPL',' ',' ',
c     $                 'recl=5670 status=new',' ',' ')
c               else
                  call set_rflist(lout, fname_out, 'LOCAL',' ','RPL',' ',' ',
     $                 'recl=5670 status=new',' ',' ')
c               endif
               call skopenf(lout,1,'Z',ierr2)
            endif

         endif
      endif

C
C+ Prepare hbook for standard diagnostic test histograms
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         call sktest_init
      endif
C
C+ Prepare fiTQun scattering table root tree
C
      if (OUTPUT_SCATTABLE .eq. 1) then
         call initscattable
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
     &     (OUTPUT_ZBS)) then
         if (OUTPUT_TYPE.eq.2) then ! atmpd mode
            CALL skclosef(lout)
         else ! le mode 
            if (SK_FILE_FORMAT.eq.0) then ! ZBS format
               CALL skclosef(lout)
            else
               call skroot_close(lout)
               call skroot_end()
            endif
         endif
      endif

            

      If (iflgv.eq.-2) close(99)
      If (OUTPUT_HBOOK.eq.OUT_SEG) call seg_closef
C
C+ Close Standard HBOOK diagnostic file
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         call sktest_close
      endif

C
C+ Close scattering table root file
C
      if (OUTPUT_SCATTABLE .eq. 1) then
         call writescattable
      endif

C
C+  Standard GEANT3 termination
C
      Call glast
C
      End

