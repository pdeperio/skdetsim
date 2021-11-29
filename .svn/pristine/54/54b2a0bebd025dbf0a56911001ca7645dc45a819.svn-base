************************************************************************
*     --------------------
*     INCLUDE 'fadcwrbnk.h'
*     --------------------
*     (Purpose)
*       Set Run Control Parameters
*
*     (Creation Date and Author)
*       2001.10.29 ; T.Ooyabu
*       2004.04.30 ; Modified and put on CVS by M.Ishitsuka
*
************************************************************************
*     FADC Bank Control
      logical IFADCI
      common /fadcbnkctrl/IFADCI

*     FADC Event Flag
*     flag_fadc_sube = 0   --> The event has no sub-event
*                    = 1   --> The event has a sub-event
      Integer flag_fadc_sube
      common /fadcevent/ flag_fadc_sube

*     FADC parameters
      Integer nwrthreshold
      PARAMETER (nwrthreshold=1)
      Integer nwrpedestal 
      PARAMETER (nwrpedestal = 0)
      Real    rwrcount2mv
      PARAMETER (rwrcount2mv = 0.255)

C----- FADC global information
      Integer npedestal 
      Real    count2mv
      Integer nthreshold
      common /fadchdbnk/ npedestal
     &                 , count2mv, nthreshold

