************************************************************************
*     -----------------
      SUBROUTINE APRDCD
*     -----------------
*
*     (Purpose)
*       Read VECT card and setup common blocks.
*       Called by : SGINIT
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*
************************************************************************
*****************************************
* include file for ap-specific control
*
*
*
*
#include "dsbnkf.h"
#include "geant321/gcflag.inc"

#include "guflag.h"
#include "apcard.h"
#include "aptrkseg.h"
#include "skdetsim_parm.h"

#ifdef SITE_HAS_NEUT
#include "neutparams.h"
#endif

*****************************************
*"apcard.h"
*      INTEGER*4 IAPHADR,IAPION,RAPITHR
*      INTEGER*4 IAPIMOD,IAPIFSIVAR
*      INTEGER*4 IFLUKA,IGHEIS,ICALOR
*      PARAMETER (IFLUKA=1)
*      PARAMETER (IGHEIS=2)
*      PARAMETER (ICALOR=3)
*
*      COMMON /APCARD/IAPHADR,IAPION,RAPITHR,IAPIMOD,IAPIFSIVAR
*
*      INTEGER*4 ISKIPEV
*      COMMON /APCCRD/ISKIPEV
*****************************************

      CHARACTER*40 TMPFIL, COMMAND
      integer lnblnk
      byte iname(132)
      equivalence (iname, seg_filebuff)

C -- initialize parameter

      IAPHADR=ICALOR
      IAPION =1
      RAPITHR=500.
      IAPIMOD=0
      IAPIFSIVAR=0
      IAPPHONUC     = 1
      PNSCALE = 1.

      IAPMORETRUTH  = 0

      ISKIPEV=0
C -- extract VECT card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('AP', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=95, FILE=TMPFIL,STATUS='OLD')
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/ln -sf '//TMPFIL//' fort.95')
#endif
      CALL FFINIT(0)
      CALL FFSET('LINP', 95) 
      CALL FFKEY('HADR' ,IAPHADR, 1,'INTEGER')
      CALL FFKEY('PION' ,IAPION , 1,'INTEGER')
      CALL FFKEY('PITHR',RAPITHR, 1,'REAL')
      CALL FFKEY('PIMOD' ,IAPIMOD , 1,'INTEGER')
      CALL FFKEY('PIFSIVAR' ,IAPIFSIVAR, 1,'INTEGER')
      CALL FFKEY('PHONUC'   ,IAPPHONUC     , 1,'INTEGER')
      CALL FFKEY('PNSCALE'  ,PNSCALE , 1,'REAL')
      CALL FFKEY('MORETRUTH',IAPMORETRUTH  , 1,'INTEGER')
      CALL FFKEY('ISKIP',ISKIPEV, 1,'INTEGER')
      CALL FFKEY('SEGF',seg_filebuff,33,'INTEGER')
      CALL FFGO
      CLOSE(UNIT=95)
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/rm fort.95')
#endif
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)


      write(*,*) ''
      write(*,*) "**************************************************"
      write(*,*) "*****        AP-specific card options        *****"
      write(*,*) "**************************************************"

      if (ISKIPEV.lt.0) ISKIPEV=0
      write(*,*) "SKIP ",ISKIPEV, "Events before generating events."
         
      if (IAPHADR.eq.IFLUKA) THEN
         write(*,*) "You selected FLUKA"
      ENDIF
      if (IAPHADR.eq.IGHEIS) THEN
         write(*,*) "You selected GHEISHA"
      ENDIF
      if (IAPHADR.eq.ICALOR) THEN
         write(*,*) "You selected GCALOR"
      ENDIF

      if (IAPION.eq.1) THEN
         write(*,*) 
     $   "Momentum(PION)<",RAPITHR,"(MeV) is traced by ORIG."
         
         if (IAPIMOD.ge.1) THEN

C           Make sure cross section table is in working directory
C           (This needs to be updated to use "findfile")
            Call getneutxs(0,dum,dum,dum,dum,dum,dum)

            if (IAPIMOD.eq.2) THEN
#ifdef SITE_HAS_NEUT
               write(*,*) ''
               write(*,*) "**************************************************"
               write(*,*) "*****      NEUT-specific card options        *****"
               write(*,*) "**************************************************"

               call NECARD

               NUSIM=0
                             
#else
               write(*,*) 
     $       "Error: Trying to use NEUT-NUCEFF but NEUT_ROOT is not set"
               stop
#endif
            end if
            
         ENDIF
         
         write(*,*) 
     $   "Momentum(PION)>",RAPITHR,"(MeV) is traced by GEANT"

      ELSE
         write(*,*) 
     $   "PION is traced by GEANT"
      ENDIF

      if (IAPMORETRUTH.NE.0) THEN
         write (*,*) "Some tracking modified for more particle truth information: IAPMORETRUTH = ", IAPMORETRUTH
      endif


      if (OUTPUT_HBOOK.eq.OUT_SEG) then
         if(TRACK_PHOTON.ne.NOT_TRACK_PHOTON) then
            write(*,*) 'Cherenkov photons and track segments will be stored'
         else 
            write(*,*) 'Only track segments (and no Cherenkov photons) will be generated'
         endif
      endif

      if (OUTPUT_HBOOK.eq.OUT_SEG) then
         if (seg_filebuff(1).ne.0) then
            lname = 1
            do while (iname(lname).ne. 0)
               lname = lname + 1
            enddo
            lname = lname - 1
            call uhtoc(seg_filebuff,4,seg_filename,lname)
            write (*,*) 'Track segment Ntuple will be written on: ',seg_filename(1:lname)
         endif
      endif

      RETURN
      END


