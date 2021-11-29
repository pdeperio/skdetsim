*****************************************
* include file for ap-specific control
*
*
*
*
*****************************************

      INTEGER*4 IAPHADR,IAPION
      INTEGER*4 IAPIMOD,IAPIFSIVAR
      INTEGER*4 IFLUKA,IGHEIS,ICALOR

*     HAT 06/2012, photonuc on/off
      INTEGER*4 IAPPHONUC
*     HAT 06/2012  photonuc scale factor
      REAL*4    pnscale

      INTEGER*4 IAPMORETRUTH  

      REAL*4    RAPITHR
      PARAMETER (IFLUKA=1)
      PARAMETER (IGHEIS=2)
      PARAMETER (ICALOR=3)

      COMMON /APCARD/IAPHADR,IAPION,RAPITHR,IAPIMOD,IAPIFSIVAR,
     +               IAPPHONUC, pnscale

      INTEGER*4 ISKIPEV
      COMMON /APCCRD/ISKIPEV

