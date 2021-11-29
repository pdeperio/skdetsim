************************************************************************
*     ------------------
*     INCLUDE 'dsparm.h'
*     ------------------
*
*       SGGATE ; gate width for high energy simulation (nsec)
*       DARKDS ; dark noise rate of PMT's (Hz)
*       TRESDS ; timing resolution (nsec)
*       QRESDS ; pulse height resolution (fraction of sigma) 
*       TWINDS ; hitsum pulse window (nsec)
*       TRIGDS ; trigger threshold (hit)
*       GATEDS ; dynamic range (nsec)
*       GLTRDS ; the timing of global trigger (nsec)
*       SIGWDS ; charge integration time (nsec)
*       DEADDS ; dead time for reject reflection (nsec)
*       OFFSET0; Default geant T0(nsec)
*
*     For TRACK_PHOTON 3 and SK1-3
*       TRGDS  ; TRIGGER TIMING in GEANT TIMING(nsec)
*       DTRGOFF; Difference between TRGDS & OFFSET(nsec)
*        T(real) = T(geant)+ OFFSET0 + DTRGOFF - TRGDS
*
*     For QBEE
*       SIGWDS_QBEE ; charge integration time (nsec)
*       DEADDS_QBEE ; dead time for reject reflection (nsec)
*       tdc_unit_sk4 ; one count of TDC in SK4 (nsec)
*       twid_block   ; time width of the one data block (nsec)
*
*       WATRDS ; water transparency
*
*     (Creation Date and Author)
*       1992.05.18 ; First version by K.S.Hirata
*         95.11.19 :      modified by Y.Koshio
*         98.11.21 :      modified by Y.Koshio (add water transp.)
*         09.05.28 :      modified by Y.Koshio (add top-bottom asym.)
*         09.08.06 :      added IWATERP for tuning studies by M. Shiozawa
*
************************************************************************

c      integer maxpm
c      parameter (maxpm=11146)

      INTEGER DOPMTDARKRT
      character*200 PMTNOISEFILE
      integer pmtnoiserun
   
      REAL  SGGATE, DARKDS, TRESDS, QRESDS, TWINDS, TRIGDS, GATEDS,
     &   GLTRDS, SIGWDS, DEADDS, OFFSET0, TRGDS, DTRGOFF,
     &   SIGWDS_QBEE, DEADDS_QBEE, DARKRT, GAIN_CORR

      integer ADD_AFTER_PULSE

      COMMON /DSPARM/ SGGATE, DARKDS, TRESDS, QRESDS, TWINDS, 
     &                TRIGDS, GATEDS, GLTRDS, SIGWDS, DEADDS,
     &                OFFSET0, TRGDS, DTRGOFF,
     &                SIGWDS_QBEE, DEADDS_QBEE, 
     &                DOPMTDARKRT, PMTNOISEFILE, DARKRT(11146),
     &                GAIN_CORR(11146), ADD_AFTER_PULSE

      INTEGER IMWATR
      REAL WATRDS
      COMMON /DSWATRA/ IMWATR, WATRDS

      INTEGER IWATERP
      COMMON /DSWTRTUNE/ IWATERP

      REAL tdc_unit_sk4, twid_block
      parameter(tdc_unit_sk4=0.520833)
      parameter(twid_block = 1024*TDC_unit_sk4*32)      

      INTEGER IMTBA
      REAL TBADS
      COMMON /DSTBA/ IMTBA, TBADS

      INTEGER nrunsk4_for_gain
      COMMON /DSGAIN/ nrunsk4_for_gain

