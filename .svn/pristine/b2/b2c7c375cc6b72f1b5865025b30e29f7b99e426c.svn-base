************************************************************************
*     ----------------------
*      FUNCTION RFPMSG(COSTH,ipol)
      SUBROUTINE RFPMSG(COSTH,ipol,prrefl,prabsl)
*     ----------------------
*
*     (Purpose)
*        Get the absorptance and reflectance of PMT photocathode
*
*     (Input)
*        COSTH  : input angle
*
*     (Output)
*CCCC        RFPMSG  : reflectance
*        prrefl  : reflectance
*        prabsl  : absorptance (used after since SK-III)
*
*     (Creation Date and Author)
*       1996. 2.28 ; First version by Y.Koshio
*         ref. (1) the reflection exp. by Y.Suzuki & K.Okumura 
*              (2) the calculation of air to water by S.Hatakeyama
*       2001. 7.24 ; Modified "tuned" by Shibata
*       2008. 3.xx : modified by Byeongsu Yang for SK-III LE mode
*         function -> subroutine
*         variable  rfpmsg -> prrefl
*         output : rfpmsg -> prrefl,prabsl
*       2008. 4.24 : LE_TUNE can be used at AP_TUNE. Byeongsu Yang
*       July 2012  : Changes following Rayleigh scattering bug by T. Akiri
*       2020. 5.20 : Upated for SK-V with new tuned values by J.Kameda
************************************************************************
      IMPLICIT NONE

#include "skdetsim_parm.h"
#include "skhead.h"
#include "geant321/gctrak.inc"
cM Fechner
#include "reflections.h"

      real rs(200),rp(200)
* Changes by Tarek Akiri in Dec 2012
* Tuning with SK4 has laser with 445 nm
* Need to increase the arrays 83 --> 108 
*      real rs_le(83,200),rp_le(83,200),as_le(83,200),ap_le(83,200),abst0(83)
      real rs_le(108,200),rp_le(108,200),as_le(108,200),ap_le(108,200),abst0(108)
      integer ical
      COMMON /RFABPM/ rs_le,rp_le,as_le,ap_le,abst0,ical
      real deg_dummy, wav_dummy,ts_dummy,tp_dummy


      integer i,j

      real ang,costh,blamda,tunedat
      real prrefl,prabsl,prabsl0
      integer iang,ipol,iwav

      save rs,rp
      real refamp
      
      character*128 fname

C      data tuned/0.467,0.533,0.933,1.133,1.067/  ! default
C       data tuned/0.388,0.533,0.60,1.133,1.067/   ! final tuning 1
c March 30,2007  : SK-III ATMPD reflection tuning
c commented out by Maximilien Fechner -- see reflections.h and /doc/
c         data tuned/0.467,0.570,0.650,1.133,1.067/  
c      data at/365.,400.,420.,460.,500./

      data refamp/0.8/
      data rs/
     &   0.4519,   0.3427,   0.3748,   0.3818,   0.3774,
     &   0.3653,   0.3762,   0.3761,   0.3724,   0.3716,
     &   0.3797,   0.3698,   0.3737,   0.3745,   0.3768,
     &   0.3757,   0.3748,   0.3730,   0.3756,   0.3783,
     &   0.3742,   0.3782,   0.3759,   0.3768,   0.3749,
     &   0.3783,   0.3747,   0.3775,   0.3810,   0.3773,
     &   0.3780,   0.3667,   0.3794,   0.3781,   0.3785,
     &   0.3806,   0.3816,   0.3829,   0.3821,   0.3822,
     &   0.3789,   0.3860,   0.3828,   0.3798,   0.3879,
     &   0.3846,   0.3853,   0.3841,   0.3853,   0.3897,
     &   0.3837,   0.3881,   0.3909,   0.3867,   0.3895,
     &   0.3964,   0.3888,   0.3951,   0.3866,   0.3944,
     &   0.3953,   0.3974,   0.3844,   0.3991,   0.3922,
     &   0.4002,   0.3888,   0.4180,   0.3682,   0.4712,
     &   0.2237,   0.2568,   0.4516,   0.3836,   0.4362,
     &   0.3959,   0.4026,   0.4041,   0.4193,   0.4016,
     &   0.4138,   0.3991,   0.4158,   0.4086,   0.4241,
     &   0.4140,   0.4131,   0.4124,   0.4143,   0.4196,
     &   0.4267,   0.4131,   0.4210,   0.4230,   0.4295,
     &   0.4200,   0.4329,   0.4270,   0.4287,   0.4262,
     &   0.4273,   0.4387,   0.4276,   0.4340,   0.4416,
     &   0.4368,   0.4368,   0.4385,   0.4399,   0.4416,
     &   0.4432,   0.4446,   0.4453,   0.4474,   0.4539,
     &   0.4540,   0.4457,   0.4605,   0.4477,   0.4620,
     &   0.4608,   0.4590,   0.4562,   0.4719,   0.4555,
     &   0.4673,   0.4693,   0.4796,   0.4701,   0.4753,
     &   0.4761,   0.4788,   0.4812,   0.4779,   0.4918,
     &   0.4816,   0.4827,   0.4917,   0.4950,   0.5009,
     &   0.4683,   0.4981,   0.5531,   0.5015,   0.4951,
     &   0.7366,   0.4347,   0.5608,   0.5036,   0.5570,
     &   0.5174,   0.5230,   0.5125,   0.5319,   0.5412,
     &   0.5475,   0.5462,   0.5447,   0.5506,   0.5428,
     &   0.5651,   0.5632,   0.5631,   0.5636,   0.5713,
     &   0.5879,   0.5794,   0.5937,   0.6099,   0.5985,
     &   0.6062,   0.6094,   0.6187,   0.6171,   0.6381,
     &   0.6326,   0.7124,   0.6373,   0.6779,   0.6621,
     &   0.6752,   0.6940,   0.7025,   0.7123,   0.7436,
     &   0.7439,   0.7355,   0.7633,   0.7413,   0.7814,
     &   0.7976,   0.8012,   0.8294,   0.8503,   0.8709,
     &   0.8892,   0.9091,   0.9298,   0.9531,   0.9958/
      data rp/
     &   0.4519,   0.3426,   0.3747,   0.3816,   0.3772,
     &   0.3650,   0.3756,   0.3755,   0.3715,   0.3706,
     &   0.3783,   0.3683,   0.3718,   0.3724,   0.3742,
     &   0.3728,   0.3716,   0.3693,   0.3715,   0.3737,
     &   0.3691,   0.3726,   0.3699,   0.3701,   0.3677,
     &   0.3704,   0.3664,   0.3684,   0.3710,   0.3668,
     &   0.3671,   0.3570,   0.3665,   0.3648,   0.3644,
     &   0.3654,   0.3655,   0.3659,   0.3643,   0.3634,
     &   0.3598,   0.3649,   0.3614,   0.3580,   0.3636,
     &   0.3598,   0.3593,   0.3575,   0.3574,   0.3599,
     &   0.3538,   0.3564,   0.3575,   0.3530,   0.3541,
     &   0.3583,   0.3512,   0.3546,   0.3473,   0.3518,
     &   0.3509,   0.3513,   0.3411,   0.3500,   0.3441,
     &   0.3480,   0.3394,   0.3564,   0.3242,   0.3864,
     &   0.2334,   0.2503,   0.3718,   0.3264,   0.3572,
     &   0.3324,   0.3339,   0.3341,   0.3416,   0.3298,
     &   0.3353,   0.3263,   0.3334,   0.3279,   0.3344,
     &   0.3280,   0.3262,   0.3241,   0.3241,   0.3247,
     &   0.3265,   0.3190,   0.3214,   0.3202,   0.3215,
     &   0.3168,   0.3194,   0.3162,   0.3153,   0.3133,
     &   0.3119,   0.3143,   0.3095,   0.3105,   0.3105,
     &   0.3082,   0.3068,   0.3060,   0.3051,   0.3043,
     &   0.3033,   0.3023,   0.3014,   0.3007,   0.3003,
     &   0.2991,   0.2972,   0.2977,   0.2957,   0.2958,
     &   0.2947,   0.2938,   0.2930,   0.2925,   0.2918,
     &   0.2911,   0.2904,   0.2895,   0.2893,   0.2887,
     &   0.2883,   0.2877,   0.2873,   0.2878,   0.2855,
     &   0.2874,   0.2868,   0.2861,   0.2869,   0.2847,
     &   0.2924,   0.2892,   0.2728,   0.2867,   0.3047,
     &   0.2054,   0.3295,   0.2712,   0.2991,   0.2793,
     &   0.2967,   0.2912,   0.3042,   0.2981,   0.2969,
     &   0.2954,   0.3051,   0.3060,   0.3052,   0.3165,
     &   0.3076,   0.3179,   0.3199,   0.3224,   0.3330,
     &   0.3283,   0.3406,   0.3335,   0.3448,   0.3541,
     &   0.3587,   0.3648,   0.3707,   0.3764,   0.3900,
     &   0.4058,   0.3954,   0.3924,   0.4260,   0.4390,
     &   0.4467,   0.4625,   0.4755,   0.4880,   0.4918,
     &   0.5203,   0.5373,   0.5534,   0.5870,   0.6014,
     &   0.6145,   0.6544,   0.6835,   0.6990,   0.7398,
     &   0.7783,   0.8205,   0.8774,   0.9056,   0.9585/


CC   2008. 4.24 : at SK_III AP_TUNE, can be used. Byeongsu Yang
* Changes by Tarek Akiri in Dec 2012
* Revert PMT reflection calculation to previous method
* Use of new reflection file
        
      IF (SK_GEOMETRY.ge.SK_III .and. REFLECTION_TUNE .NE. AP_TUNE  ) THEN

        if(ical.eq.0)then
            ical=1

          if (SK_GEOMETRY .eq. SK_V) then
            call findconsts(fname,'reflection table',
     $            'reflect.feb21.dat',0)
          else
            call findconsts(fname,'reflection table',
     $           'reflect.dec12.dat',0)
          endif

          open(unit=97,file=fname,status='old',
     $           form='formatted')
         
*            do i=1,83
            do i=1,108
             do j=1,200
               read(97,*)wav_dummy,deg_dummy,rs_le(i,j),rp_le(i,j),
     $               ts_dummy,tp_dummy,as_le(i,j),ap_le(i,j)
               if ( j .eq. 1 ) abst0(i) = 0.5*(ap_le(i,1)+as_le(i,1))
             enddo
            enddo
            close(97)
         endif
 
c     ang in degree
         ANG=ACOS(COSTH)*180./3.141593
c     iang=int(ang*200./90.)
         iang=int(ang*200./90.) + 1
         if(iang.gt.200) then
           iang = 200
         else if(iang.lt.1) then
           iang = 1
         endif

         blamda=1237.8E-09 / GEKIN

         iwav=int(blamda-336.0)
*         if(iwav.gt.83) then
*           iwav = 83
         if(iwav.gt.108) then
           iwav = 108
         else if(iwav.lt.1) then
           iwav = 1
         endif
 
         if(ipol.eq.1) then
           prrefl=rp_le(iwav,iang)
           prabsl=ap_le(iwav,iang)
         else if(ipol.eq.2) then
           prrefl=rs_le(iwav,iang)
           prabsl=as_le(iwav,iang)
         else if(ipol.eq.0) then
           prrefl=0.5*(rp_le(iwav,iang)+rs_le(iwav,iang))
           prabsl=0.5*(ap_le(iwav,iang)+as_le(iwav,iang))
         endif


      ELSE IF(REFLECTION_TUNE .EQ. LE_TUNE) THEN ! for SK-I/II LE tune
         if(ical.eq.0)then
            ical=1
            call findconsts(fname,'reflection table',
     $           'reflect.dec98.dat',0)
c            open(97,file='/home/sklowe/const/reflect.dec98.dat',
            open(97,file=fname,
     &           form='formatted',status='old')
             i=1
             do j=1,200
               read(97,*)deg_dummy,rs_le(i,j),rp_le(i,j),
     &               ts_dummy,tp_dummy,as_le(i,j),ap_le(i,j)
             enddo
            close(97)
         endif
         
c     ang in degree
         ANG=ACOS(COSTH)*180./3.141593
c     iang=int(ang*200./90.)
         iang=int(ang*200./90.) + 1
         if(iang.gt.200) then
           iang = 200
         else if(iang.lt.1) then
           iang = 1
         endif

         if(ipol.eq.1) then
           prrefl=rp_le(1,iang)
C           prabsl=ap_le(1,iang)
         else if(ipol.eq.2) then
           prrefl=rs_le(1,iang)
C           prabsl=as_le(1,iang)
         else if(ipol.eq.0) then
           prrefl=0.5*(rp_le(1,iang)+rs_le(1,iang))
C           prabsl=0.5*(ap_le(1,iang)+as_le(1,iang))
         endif

c
C         blamda=1237.8E-09 / GEKIN
c
         
      ELSE IF(REFLECTION_TUNE .EQ. AP_TUNE) THEN
c     ang in degree
         ANG=ACOS(COSTH)*180./3.141593
         iang=int(ang*200./90.)
         if ( iang.le.0 ) then
            iang=1
         else if ( iang.ge.200 ) then
            iang=200
         endif

C         refp=rp(iang)
C         refs=rs(iang)
C         if(ipol.eq.1)prrefl=refp
C         if(ipol.eq.2)prrefl=refs
C         if(ipol.eq.0)prrefl=0.5*(refp+refs)

C         prabsl0=0.5*(ap(1)+as(1))
         if(ipol.eq.1) then
           prrefl=rp(iang)
C           prabsl=ap(iang)
         else if(ipol.eq.2) then
           prrefl=rs(iang)
C           prabsl=ap(iang)
         else if(ipol.eq.0) then
           prrefl=0.5*(rp(iang)+rs(iang))
C           prabsl=ap(iang)
         endif

     
         blamda=1237.8E-09 / GEKIN
     
CMM   Add refamp
c         prrefl=refamp*prrefl*tunedat(blamda,tuned,at,0)
c M Fechner
         prrefl=refamp*prrefl*tunedat(blamda,tuned_pm,at,0,ntuning)
      ELSE
         write(*,*) 'RFPMSG: Combination of Reflection tune and SK geometry is incorrect. Stop.'
         stop
      ENDIF

      RETURN
      END

************************************************************************
*     ----------------------
      FUNCTION PRABSL0(blamda)
*     ----------------------
*
*     (Purpose)
*        Get the absorptance of photocathode at normal incident, 0 degree
*
*     (Input)
*        blamda : wavelength
*
*     (Output)
*        prabsl0 : average absorbance around 0 degree(0<= <0.45deg)
*                  (used after since SK-III) (normal incident)
*
*     (Creation Date and Author)
*       2008. 4.23 : Byeongsu Yang
*         output : prabsl0
************************************************************************
      IMPLICIT NONE
#include "skdetsim_parm.h"
#include "skhead.h"
      integer i,j

*      real rs_le(83,200),rp_le(83,200),as_le(83,200),ap_le(83,200),abst0(83)
      real rs_le(108,200),rp_le(108,200),as_le(108,200),ap_le(108,200),abst0(108)
      integer ical
      COMMON /RFABPM/ rs_le,rp_le,as_le,ap_le,abst0,ical
      real deg_dummy, wav_dummy,ts_dummy,tp_dummy

      real blamda, prabsl0
      integer iwav

      character*128 fname

      if(ical.eq.0)then
         ical=1

         if (SK_GEOMETRY .eq. SK_V) then
         call findconsts(fname,'reflection table',
     $           'reflect.feb21.dat',0)
         else
         call findconsts(fname,'reflection table',
     $           'reflect.dec12.dat',0)
         endif         

         open(unit=97,file=fname,status='old',
     $           form='formatted')

*     do i=1,83
         do i=1,108
          do j=1,200
            read(97,*)wav_dummy,deg_dummy,rs_le(i,j),rp_le(i,j),
     $            ts_dummy,tp_dummy,as_le(i,j),ap_le(i,j)
            if ( j .eq. 1 ) abst0(i) = 0.5*(ap_le(i,1)+as_le(i,1))
          enddo
         enddo
         close(97)
      endif
 
      iwav=int(blamda-336.0)
*      if(iwav.gt.83) then
*        iwav = 83
      if(iwav.gt.108) then
        iwav = 108
      else if(iwav.lt.1) then
        iwav = 1
      endif

      prabsl0 = abst0(iwav)
 
      RETURN
      END
