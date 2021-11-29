CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       
C     Include file to control SK detector simulation
C     ==============================================
C
C  --------------------------------------------------------------- 
C     Arguments            Purpose
C  --------------------------------------------------------------- 
C     TRACK_PHOTON         Control tracking photons
C     (INTEGER*4)          = NOT_TRACK_PHOTON (=-1) Not track photon at all 
C                          = TRACK_ALL_PHOTON (= 1) Track all photons
C                          = TRACK_QE_PHOTON  (= 2) Apply Q.E. and only
C                                             accepted photons are tracked
C                          = TRACK_QE_OVER_ABS_PHOTON  (= 3) Apply QE/ABS
C                                  and only accepted photons are tracked
C     WATER_TUNE           Water Parameter tuning
C     (INTEGER*4)          = LE_TUNE (= 1)  old LE group tuning         
C                          = AP_TUNE (= 2)  old ATMPD group tuning
C                               for only TRACK_QE_PHOTON
C                          = UNI_TUNE (= 3) tuning
C                               for TRACK_QE_OVER_ABS_PHOTON
C			   = GD_SCATTERING (= 10) Gd Effect due to scattering in SK-IV		
C				Works same as UNI_TUNE otherwise
C                          = GD_ABSORPTION (= 11) Gd Effect due to absorption in SK-IV
C                               Works same as UNI_TUNE otherwise
C
C
C     REFLECTION_TUNE      Reflection tuning
C     (INTEGER*4)          = LE_TUNE (= 1)  old LE group tuning         
C                          = AP_TUNE (= 2)  old ATMPD group tuning
C                               for TRACK_QE_PHOTON
C                          = UNI_TUNE (= 3) tuning
C                               for TRACK_QE_OVER_ABS_PHOTON
C     DS_TUNE              Tuning in DS routines
C     (INTEGER*4)          = LE_TUNE (= 1)  keep all p.e. and only a gate
C                          = AP_TUNE (= 2)  throw away p.e. out of gate and only a gate
C                          = LE_MULTI_TUNE (= 3) keep all p.e. and multi-gate
C                          = AP_MULTI_TUNE (= 4) throw away p.e. out of gate and multi-gate
C     OUTPUT_ZBS           Save ZBS data
C     OUTPUT_TYPE          For ATMPD group needs some additional bank
C     (INTEGER*4)          = LE group    (= 1) not add those banks
C                          = ATMPD group (= 2) add those banks
C     OUTPUT_HBOOK         Save HBOOK histograms
C     (INTEGER*4)          = No HBOOK ouput (= 0)
C                          = OUT_SEG (= 1)  Output for fastsim segments
C                          = OUT_TEST (= 2) Output standard diagnostic tests
C     OUTPUT_TEXT          Save TEXT records
C
C     OUTPUT_SCATTABLE     Output fiTQun scattering table root tree
C
C     VALIDATE             = Produce output for crosschecks
C  ---------------------------------------------------------------
C
      INTEGER*4 TRACK_PHOTON, WATER_TUNE, REFLECTION_TUNE, DS_TUNE 
      INTEGER*4 TRACK_PHOTON_ORG
      INTEGER*4 NOT_TRACK_PHOTON, TRACK_ALL_PHOTON, TRACK_QE_PHOTON
     &     ,TRACK_QE_OVER_ABS_PHOTON, TRACK_QE_OVER_ABS_SAVE_MEMORY
      INTEGER*4 LE_TUNE, AP_TUNE, LE_MULTI_TUNE, AP_MULTI_TUNE, UNI_TUNE
      INTEGER*4 FLAG_TRACK_SW, GD_ABSORPTION, GD_SCATTERING
      INTEGER*4 GD_MIX, GD_EGADS, GD_GADZOOKS
      REAL    GD_SCAT_FRAC
      LOGICAL OUTPUT_ZBS, OUTPUT_TEXT, VALIDATE
      INTEGER*4 OUTPUT_HBOOK, OUT_SEG, OUT_TEST, OUTPUT_TYPE
      INTEGER*4 OUTPUT_SCATTABLE
      PARAMETER (NOT_TRACK_PHOTON=-1)
      PARAMETER (TRACK_ALL_PHOTON=1)
      PARAMETER (TRACK_QE_PHOTON=2)
      PARAMETER (TRACK_QE_OVER_ABS_PHOTON=3)
      PARAMETER (TRACK_QE_OVER_ABS_SAVE_MEMORY=4)
      PARAMETER (LE_TUNE=1)
      PARAMETER (AP_TUNE=2)
      PARAMETER (UNI_TUNE=3)
      PARAMETER (GD_SCATTERING=10)
      PARAMETER (GD_ABSORPTION=11)
      PARAMETER (GD_MIX=13)
      PARAMETER (GD_EGADS=23)
      PARAMETER (GD_GADZOOKS=42)

      PARAMETER (LE_MULTI_TUNE=3)
      PARAMETER (AP_MULTI_TUNE=4)
      PARAMETER (OUT_SEG=1)
      PARAMETER (OUT_TEST=2)
      INTEGER*4 SK_I, SK_II, SK_III, SK_IV, SK_V, SK_VI
      PARAMETER (SK_I   = 1)
      PARAMETER (SK_II  = 2)
      PARAMETER (SK_III = 3)
      PARAMETER (SK_IV  = 4)
      PARAMETER (SK_V   = 5)  
      PARAMETER (SK_VI  = 6)  

c added by t.wester for specifiying badch from run number    
      INTEGER BADCHRUN

c added by M.Harada(Original is by Y.Ashida, R.Akutsu, and Ka Ming)
c for implementing the Gd-G4 simulation.
C     on/off switch for GEANT4 simulation
      INTEGER*4 USEG4
C     on/off switch for Gd-simulation
      INTEGER*4 GD_WTR, G4_INITI_RAND, GD_MODEL
      REAL GD_CONC

			INTEGER*4 OD_PHOTON, COUNT_ALL_OD_PHOTON 
      PARAMETER ( COUNT_ALL_OD_PHOTON  = 1 )
      COMMON/SKDETSIM_PARM/TRACK_PHOTON, WATER_TUNE, 
     &   REFLECTION_TUNE, DS_TUNE,
     &   OUTPUT_ZBS, OUTPUT_HBOOK, OUTPUT_TEXT, OUTPUT_TYPE, 
     &   VALIDATE, OD_PHOTON, TRACK_PHOTON_ORG, FLAG_TRACK_SW,
     &   OUTPUT_SCATTABLE, BADCHRUN, USEG4, GD_WTR, GD_CONC, GD_SCAT_FRAC, G4_INITI_RAND, GD_MODEL
