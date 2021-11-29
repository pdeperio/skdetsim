************************************************************************
*     -----------------
      SUBROUTINE CNTLCARD
*     -----------------
*
*     (Purpose)
*       Read basic control variables from card and setup common blocks.
*       Called by : SGINIT
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*      2003.01.10 J.Kameda
*
************************************************************************
      implicit none

#include "dsbnkf.h"
#include "geant321/gcflag.inc"
#include "guflag.h"
#include "skdetsim_parm.h"
#include "skhead.h"

      CHARACTER*40 TMPFIL, COMMAND
      integer iret
C -- extract SKCNTL card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')      
      CALL EXCARD('SKCNTL', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=92, FILE=TMPFIL,STATUS='OLD')
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/ln -sf '//TMPFIL//' fort.92')
#endif

C-- initializiation
      SK_GEOMETRY = 0
      TRACK_PHOTON = 0
      DS_TUNE = 0
      OUTPUT_ZBS = .TRUE.
      OUTPUT_HBOOK = 0
      OUTPUT_SCATTABLE = 0
      OUTPUT_TEXT = .FALSE.
      VALIDATE = .FALSE.
      WATER_TUNE = 0
      REFLECTION_TUNE = 0
      sk_file_format = 0
C-- by Thomas Wester
      BADCHRUN = 0

      FLAG_TRACK_SW = 0

      IFMESS = 0

#ifdef WithG4
C     on/off switch for Geant4-simulation for neutron
      USEG4 = 0
C     on/off switch for Gd-simulation
C     GD_WTR = 0 --> pure water simulation
C     GD_WTR = 1 --> calculate Gd cross-section and simulate Gd-reaction when necessary
      GD_WTR = 0 
      GD_MODEL = 0 
      GD_CONC=0.2178

C     Initial seed for GEANT4
      G4_INITI_RAND = 6435459
#endif

      CALL FFINIT(0)
      CALL FFSET('SIZE', 20)
      CALL FFSET('LINP', 92)
      CALL FFKEY('SKGEOMETRY',SK_GEOMETRY, 1,'INTEGER')
      CALL FFKEY('TRACKPHOTON',TRACK_PHOTON,1,'INTEGER')
      CALL FFKEY('DSTUNE',DS_TUNE,1,'INTEGER')
      CALL FFKEY('OUTPUTZBS',OUTPUT_ZBS,1,'LOGICAL')
      CALL FFKEY('OUTPUTTYPE',OUTPUT_TYPE,1,'INTEGER')
      CALL FFKEY('OUTPUTHBOOK',OUTPUT_HBOOK,1,'INTEGER')
      CALL FFKEY('OUTPUTSCATTABLE',OUTPUT_SCATTABLE,1,'INTEGER')
      CALL FFKEY('OUTPUTTEXT',OUTPUT_TEXT,1,'LOGICAL')
      CALL FFKEY('VALIDATE',VALIDATE,1,'LOGICAL')
      CALL FFKEY('TRACKAUTOSW',FLAG_TRACK_SW,1,'INTEGER')
      call FFKEY('FILEFORMAT',SK_FILE_FORMAT,1,'INTEGER')
C --  by Thomas Wester
      call FFKEY('BADCHRUN',BADCHRUN,1,'INTEGER')

      CALL FFKEY('WATERTUNE',WATER_TUNE,1,'INTEGER')
      CALL FFKEY('REFLECTIONTUNE',REFLECTION_TUNE, 1,'INTEGER')

      CALL FFKEY('MESSAGEOFF',IFMESS, 1,'INTEGER')

#ifdef WithG4
      CALL FFKEY('USEG4',USEG4, 1,'INTEGER')
      CALL FFKEY('GDWTR',GD_WTR, 1,'INTEGER')
      CALL FFKEY('G4INITIRAND',G4_INITI_RAND, 1,'INTEGER')
      CALL FFKEY('GDMODEL',GD_MODEL, 1, 'INTEGER')
      CALL FFKEY('GDCONC',GD_CONC, 1, 'REAL')
#endif

      CALL FFGO
      CLOSE(UNIT=92)

      if (FLAG_TRACK_SW .eq. 1) then
         TRACK_PHOTON_ORG = TRACK_PHOTON
      write(*,*) '***********************************************'
      write(*,*) '* Tracking method will automatically switched *'
      write(*,*) '* to 4 for over 1TeV neutrinos.               *'
      write(*,*) '***********************************************'
      endif


#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/rm fort.92')
#endif
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)
#ifdef UPMU
      write(6,*) "************************************************"
      write(6,*) "***  This is skdetsim for UPMU              ****"
      write(6,*) "************************************************"

#endif

      write(6,*) "**************************************************"
      write(6,*) "***  skdetsim fundamental control parameters  ****"
      write(6,*) "**************************************************"

      IF( SK_GEOMETRY.EQ.SK_I.OR.SK_GEOMETRY.EQ.SK_II.OR.SK_GEOMETRY.EQ.SK_III.OR.SK_GEOMETRY.EQ.SK_IV.OR.SK_GEOMETRY.EQ.SK_V) THEN
        WRITE(6,610) sk_geometry
 610    FORMAT(' Detector Geometry version (1:SK I 2:SK II 3: SK III 4: SK IV 5: SK V): ',I4)
      ELSE
        WRITE(0,*)'cntlcard.F : Error, not supporetd sk geometry ', sk_geometry
        STOP 1
      ENDIF

      IF(TRACK_PHOTON.EQ.NOT_TRACK_PHOTON	  .OR.
     |   TRACK_PHOTON.EQ.TRACK_ALL_PHOTON	  .OR.
     |   TRACK_PHOTON.EQ.TRACK_QE_PHOTON	  .OR.
     |   TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_PHOTON .OR.
     |   TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_SAVE_MEMORY ) THEN
        WRITE(6,620) track_photon
C 620    FORMAT(' Photon tracking (1:all 2:QE multiplied 3: no tracking): ',I4)
 620    FORMAT(' Photon tracking (-1: no tracking 1:all 2:QE multiplied 3', 
     & ': QE/Abs. multiplied 4:QE/Abs. SAVE_MEMORY):',I4)
      ELSE
        WRITE(0,*) "cntlcard.F : Error, not suported TRACKPHOTON modes , ", track_photon
        STOP 1
      ENDIF

      IF(WATER_TUNE.EQ.0) THEN
        WATER_TUNE =  TRACK_PHOTON
        IF(TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_SAVE_MEMORY) WATER_TUNE = UNI_TUNE 
        WRITE(6,*) 'No WATERTUNE in card. Set the default WATER_TUNE.'
      ENDIF
C     WATER_TUNE == 10, is Gd effect due to scattering
C     WATER_TUNE == 11, is Gd effect due to absorption
      IF(WATER_TUNE.EQ.LE_TUNE.OR.WATER_TUNE.EQ.AP_TUNE.OR.WATER_TUNE.EQ.UNI_TUNE
     & .OR.WATER_TUNE.EQ.10.OR.WATER_TUNE.EQ.11) THEN  
        WRITE(6,630) water_tune
C 630    FORMAT(' Water tuning (1: low E  2: ATMPD) : ',I4)
 630    FORMAT(' Water tuning (1: low E  2: ATMPD  3: UNI mode) : ',I4)
      ELSE
        WRITE(0,*) "cntlcard.F : Error, not suported WATERTUNE modes , ", water_tune
        STOP 1
      ENDIF

      IF(REFLECTION_TUNE.EQ.0) THEN
        REFLECTION_TUNE =  TRACK_PHOTON
        IF(TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_SAVE_MEMORY) REFLECTION_TUNE = UNI_TUNE 
        WRITE(6,*) 'No REFLECTIONTUNE in card. Set the default REFLECTION_TUNE.'
      ENDIF
      IF(REFLECTION_TUNE.EQ.LE_TUNE.OR.REFLECTION_TUNE.EQ.AP_TUNE.OR.REFLECTION_TUNE.EQ.UNI_TUNE) THEN
        WRITE(6,640) reflection_tune
C 640    FORMAT(' Reflection tuning (1: low E  2: ATMPD) : ',I4)
 640    FORMAT(' Reflection tuning (1: low E  2: ATMPD  3: UNI mode) : ',I4)
      ELSE
        WRITE(0,*) "cntlcard.F : Error, not suported REFLECTIONTUNE modes , ", reflection_tune
        STOP 1
      ENDIF

      IF(DS_TUNE.EQ.LE_TUNE.OR.DS_TUNE.EQ.AP_TUNE.OR.DS_TUNE.EQ.UNI_TUNE) THEN
        WRITE(6,650) ds_tune
C 650    FORMAT(' DS tuning (1: low E  2: ATMPD)', I4)
 650    FORMAT(' DS tuning (1: low E  2: ATMPD  3: UNI_TUNE)', I4)
      ELSE
        WRITE(0,*) "cntlcard.F : Error, not suported DSTUNE modes , ", ds_tune
        STOP 1
      ENDIF

      IF(TRACK_PHOTON .EQ. NOT_TRACK_PHOTON ) THEN

        WRITE(6,*) "******************************************************"
        WRITE(6,*) "******************************************************"
        WRITE(6,*) "******************************************************"
        WRITE(6,*) "***                                                ***"
        WRITE(6,*) "*** You select NOT TRACK PHOTON mode.              ***"
        WRITE(6,*) "*** IF you want to TRACK PHOTONS,                  ***"
        WRITE(6,*) "*** please CHANGE TRACK PHOTON mode!!!             ***"
        WRITE(6,*) "***                                                ***"
        WRITE(6,*) "******************************************************"
        WRITE(6,*) "******************************************************"
        WRITE(6,*) "******************************************************"

      ELSEIF (SK_GEOMETRY .EQ. SK_I .OR. SK_GEOMETRY .EQ. SK_II) THEN

        IF (.NOT.(
CC LE official mode
     |      (TRACK_PHOTON    .EQ. TRACK_ALL_PHOTON .AND.
     &       WATER_TUNE      .EQ. LE_TUNE          .AND.
     &       REFLECTION_TUNE .EQ. LE_TUNE          .AND.
     &       DS_TUNE         .EQ. LE_TUNE               ) .OR.
CC AP official mode
     |      (TRACK_PHOTON    .EQ. TRACK_QE_PHOTON  .AND.
     &       WATER_TUNE      .EQ. AP_TUNE          .AND.
     &       REFLECTION_TUNE .EQ. AP_TUNE          .AND.
     &       DS_TUNE         .EQ. AP_TUNE               )     )) THEN
          WRITE(6,*) "cntlcard.F : Warning, unofficial modes!!"
          WRITE(6,*) "If you want use Unofficial Modes, Please modify cntlcard.F. Stop."
          stop
        ENDIF

      ELSEIF (SK_GEOMETRY .EQ. SK_III) THEN

        IF (.NOT.(
CC AP official mode
     |      (TRACK_PHOTON     .EQ. TRACK_QE_PHOTON .AND.
     &       WATER_TUNE       .EQ. AP_TUNE         .AND.
     &       REFLECTION_TUNE  .EQ. AP_TUNE         .AND.
     &       DS_TUNE          .EQ. AP_TUNE                       ) .OR.
CC LE official mode
     |      (TRACK_PHOTON     .EQ. TRACK_QE_OVER_ABS_PHOTON .AND.
     &       (WATER_TUNE.EQ.UNI_TUNE.OR.WATER_TUNE.EQ.10.OR.WATER_TUNE.EQ.11).AND.
     &       REFLECTION_TUNE  .EQ. UNI_TUNE                 .AND.
     &       DS_TUNE          .EQ. UNI_TUNE                      )     )) THEN
          WRITE(6,*) "cntlcard.F : Warning, unofficial modes!!"
          WRITE(6,*) "If you want use Unofficial Modes, Please modify cntlcard.F. Stop."
          WRITE(6,*) "TRACK, TQE_OverABS_photon" 
          WRITE(6,*) TRACK_PHOTON, TRACK_QE_OVER_ABS_PHOTON
          WRITE(6,*) "WATER_TUNE, REFLECTION_TUNE, DS_TUNE, UNI_TUNE"
          WRITE(6,*) WATER_TUNE, REFLECTION_TUNE, DS_TUNE, UNI_TUNE
          stop
        ENDIF

      ELSE IF (SK_GEOMETRY .EQ. SK_IV .OR. SK_GEOMETRY .EQ. SK_V ) THEN
CC forbidden mode
CC        IF ( ) THEN
CC          WRITE(6,*) "cntlcard.F : Error. Forbidden mode!!"
CC          STOP 1

        IF (.NOT.(
CC AP official mode
     |      (TRACK_PHOTON     .EQ. TRACK_QE_OVER_ABS_PHOTON .AND.
     &       (WATER_TUNE.EQ.UNI_TUNE.OR.WATER_TUNE.EQ.10.OR.WATER_TUNE.EQ.11)  .AND.
     &       ((REFLECTION_TUNE  .EQ. AP_TUNE).or.(REFLECTION_TUNE .eq. UNI_TUNE)) .AND.
     &       ((DS_TUNE          .EQ. AP_TUNE).or.(DS_TUNE         .EQ. UNI_TUNE))
     &                  ) .OR.
CC AP official mode (saving memory version)
     |      ((TRACK_PHOTON     .EQ. TRACK_QE_OVER_ABS_SAVE_MEMORY) .and.
     &       (WATER_TUNE.EQ.UNI_TUNE.OR.WATER_TUNE.EQ.10.OR.WATER_TUNE.EQ.11).and.
     &       ((REFLECTION_TUNE  .EQ. AP_TUNE).or.(REFLECTION_TUNE .eq. UNI_TUNE)) .AND.
     &       ((DS_TUNE .EQ. AP_TUNE).or.(DS_TUNE .EQ. UNI_TUNE))
     &                                                           ) .OR.
CC LE official mode ??
     |      (TRACK_PHOTON     .EQ. TRACK_QE_OVER_ABS_PHOTON .AND.
     &       (WATER_TUNE.EQ.UNI_TUNE.OR.WATER_TUNE.EQ.10.OR.WATER_TUNE.EQ.11) .AND.
     &       REFLECTION_TUNE  .EQ. UNI_TUNE                 .AND.
     &       ((DS_TUNE .EQ. LE_TUNE).or.(DS_TUNE .EQ. UNI_TUNE))
     &                                                           ))) THEN
          WRITE(6,*) "cntlcard.F : Warning, unofficial modes!!"
          WRITE(6,*) "If you want use Unofficial Modes, Please modify cntlcard.F. Stop."
          stop
        ENDIF
      ENDIF

      if (output_zbs) then
         write(6,*) "ZBS output   : ON "
         if(output_type.eq.1) then
            write(6,*) 'For Low-E'
         else if(output_type.eq.2) then
            write(6,*) 'For ATMPD'
         else
            write(0,*) 'cntlcard.F : Error, OUTPUTTYPE is wrong, ',output_type,'. Please set 1 or 2'
            stop 1
         endif
      else
         write(6,*) "ZBS output   : OFF "
      endif

      if (output_hbook.eq.OUT_SEG) then
         write(6,*) "HBOOK output : ON for SEG"
      else if(output_hbook.eq.OUT_TEST) then
         write(6,*) "HBOOK output : ON for Diagnostic Test "
      else if(output_hbook.eq.0) then
         write(6,*) "HBOOK output : OFF "
      else
         write(0,*) 'cntlcard.F : Error, OUTPUTHBOOK is wrong, ',output_hbook,'. Please set 0/1/2'
         stop 1
      endif

      if (output_scattable.eq.1) then
         write(6,*) "SCATTABLE output : ON "
      else if(output_scattable.eq.0) then
         write(6,*) "SCATTABLE output : OFF "
      else
         write(0,*) 'cntlcard.F : Error, OUTPUTSCATTABLE is wrong, ',output_scattable,'. Please set 0/1'
         stop 1
      endif

      if (output_text) then
       write(6,*) "TEXT output  : ON "
      else
       write(6,*) "TEXT output  : OFF "
      endif

      if (validate) then
       write(6,*) "Validation output  : ON "
      else
       write(6,*) "Validation output  : OFF "
      endif

cc check SK_FILE_FORMAT

      if ((SK_FILE_FORMAT .ne. 0) .and. (SK_FILE_FORMAT .ne. 1)) then
         write(6,*) 'Error: FILEFORMAT is not correctly set in card.', SK_FILE_FORMAT
         stop
      endif

      RETURN
      END
