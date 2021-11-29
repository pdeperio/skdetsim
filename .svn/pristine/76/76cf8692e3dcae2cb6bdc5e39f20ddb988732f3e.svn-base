************************************************************************
*     ------------------
*     INCLUDE 'dsparma.h'
*     ------------------
*
*       FULWIDA ; full gate width for anti simulation (nsec)
*       DSANTT0 ; Antidetector T0 location (nsec)
*                 If -16 <-> +16 then DSANTT0 should be 16
*                 If -10 <-> +6  then DSANTT0 should be 10
*       DARKDSA ; dark noise rate of PMT's (Hz)
*       DSTYVKR ; (obsolete)Relative ratio for reflectivity of tyvek
*       APMTCOR ; PMT's correction efficiency (Max 1.)
*       APMTEFR ; Ratio of efficiency for old/new tubes (assume <1)
*       APMTEFT,B,W ; Relative efficiencies for top, bottom, wall tubes
*          (taking into account unsimulated absorption)
*       AWLSCOR ; Wave Length Shifter's correction efficiency(Max 1.)
*       ADSCTHR ; Anti discriminator threshold(unit:pe)
*       ADSTHRF ; Anti PMT threshold factor for old and new pmts(for thresh in pe)
*       ADRANGE ; (obsolete)Anti detector charge dynamic range, new tubes
*       ADSAPO  ; Anti detector saturation parameters, old tubes
*       ADSAPN  ; Anti detector saturation parameters, new tubes
*       IAGEQAL ; Quality of anti detector's geometry (Quick:1/Exact:2)
*       DSTYVKRTOP,BOT,WAL; Relative ratio of Tyvek reflection, for TOP, BOT, WAL
*       DSTYVKM ; Tyvek reflectivity model
*       DSSEGRTOP,BOT,WAL; Relative ratio of segmentation Tyvek reflection TOP, BOT, WAL(not used)
*       PROBTT,B   ; probability to have transmitted photon by TYVEK on
*                         top and bottom
*       GAINFACT ; Gaussian portion of response for old od-tubes
*       GAIN_PLATEAU ; Plateau portion of response for old od-tubes  
*       ODDARKRUN ; Use dark rates from specified file for OD
*       ODINDSPE ; If set use individual PMT SPE, default is no set (so use PMT type SPE)  
*       ODPEDSMEAR ; Pedestal smearing for OD for SK5, default = 0 in SK5 (already included 
*       							in SPE), in SK4 this is hardcoded as 0.12 in dsant_sk4.F
*       ODIMBSMEAR ; Extra charge smearing for IMB (also referred as "old" in the code...)
*      									  PMTs in SK5, default is 0 in SK5,
*       									in SK4 this is hardcoded as 0.2 in dsanti_sk4.F
*       PROBODSST ; effective absorption on stainless steel in top od region
*       PROBODSSB ; effective absorption on stainless steel in bottom od region
************************************************************************

      REAL*4 FULWIDA,DSANTT0,DARKDSA(3),DSTYVKR,APMTCOR,APMTEFR,
     $       APMTEFT,APMTEFB,APMTEFW,
     $       AWLSCOR,ADSCTHR,ADSTHRF,ADRANGE,
     $       ADSAPO(3), ADSAPN(3),
     $       DSTYVKRTOP(2),DSTYVKRBOT(2),DSTYVKRWAL(2),DSTYVKM,
     $       DSSEGRTOP,DSSEGRBOT,DSSEGRWAL,
     $       PROBTT, PROBTB, GAINFACT, GAIN_PLATEAU,
     $       ODPEDSMEAR, ODIMBSMEAR, PROBODSST, PROBODSSB
      INTEGER*4 IAGEQAL, ODDARKRUN
      LOGICAL ODINDSPE

      COMMON /DSPARMA/FULWIDA,DSANTT0,DARKDSA,DSTYVKR,APMTCOR,
     $           APMTEFR,APMTEFT,APMTEFB,APMTEFW,
     $           AWLSCOR,ADSCTHR,ADSTHRF,ADRANGE,
     $           ADSAPO,ADSAPN, IAGEQAL,
     $           DSTYVKRTOP,DSTYVKRBOT,DSTYVKRWAL,DSTYVKM,
     $           DSSEGRTOP,DSSEGRBOT,DSSEGRWAL,
     $           PROBTT, PROBTB, GAINFACT, GAIN_PLATEAU,
     $           ODDARKRUN, ODINDSPE, ODPEDSMEAR, ODIMBSMEAR,
     $           PROBODSST, PROBODSSB
