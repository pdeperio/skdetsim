************************************************************************
*     -------------------------------------
      SUBROUTINE DS
*     -------------------------------------
*
*     (Purpose)
*       Detector Simulator (in geant)
*
*     (Creation Date and Author)
*       1995.01.27 ; First version by Y.Koshio
*       2003.01.07 ; Modified for SK-II by Y.Koshio
*       2008.10.29 ; Modified for SK-IV by J.Kameda
*                    New tracking methods are supported.
************************************************************************

      implicit none

#include "skhead.h"
#include "skdetsim_parm.h"

#include "guflag.h"

      IF(TRACK_PHOTON.ne.NOT_TRACK_PHOTON) THEN

         if (TRACK_PHOTON.eq.TRACK_ALL_PHOTON) THEN

            if (SK_GEOMETRY .eq. SK_IV .OR. SK_GEOMETRY .eq. SK_V) then

               call DSTRGNUM   
               call DSNOIS_MULTI
               call DSRES_MULTI
               if(ifanti.ne.0) call dsanti_sk4
               call DSSWTRG
                call dsdigi_sk4
               if(ifanti.ne.0) call dsdigi_anti_sk4

            elseif ((SK_GEOMETRY .eq. SK_I) .or. (SK_GEOMETRY .eq. SK_II) 
     &         .or. (SK_GEOMETRY .eq. SK_III) ) then

               CALL DSNOIS
               CALL DSRES
               CALL DSTRIG
               CALL DSANTI

            else
               write(6,*) 'ERROR:DS: SK_GEOMETRY is irregular. Stop.'
               Stop
            endif

         elseif (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_PHOTON) then

            if (SK_GEOMETRY .eq. SK_IV .OR. SK_GEOMETRY .eq. SK_V) then            
               call DSTRGNUM   
               call DSNOIS_MULTI
               call DSRES_MULTI
               if(IFANTI.ne.0) call dsanti_sk4
               call DSCPTQZTORAW
               call DSSWTRG
                call dsdigi_sk4
               if(IFANTI.ne.0) call dsdigi_anti_sk4

            elseif ((SK_GEOMETRY .eq. SK_I) .or. (SK_GEOMETRY .eq. SK_II) 
     &         .or. (SK_GEOMETRY .eq. SK_III) ) then

               CALL DSNOIS_multi
               CALL DSRES_multi
               CALL DSTRIG_multi
               CALL DSANTI

            else
               write(6,*) 'ERROR:DS: SK_GEOMETRY is irregular. Stop.'
               Stop
            endif

         elseif (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY ) THEN

            if (SK_GEOMETRY .eq. SK_IV .OR. SK_GEOMETRY .eq. SK_V) then
               call DSTRGNUM   
               call DSHIGH_SK4
               if(IFANTI.ne.0) call dsanti_sk4
               call DSCPTQZTORAW
               call DSSWTRG
                call dsdigi_sk4
               if(IFANTI.ne.0) call dsdigi_anti_sk4

            elseif ((SK_GEOMETRY .eq. SK_I) .or. (SK_GEOMETRY .eq. SK_II) 
     &         .or. (SK_GEOMETRY .eq. SK_III) ) then

               CALL DSHIGH
               CALL DSANTI

            else
               write(6,*) 'ERROR:DS: SK_GEOMETRY is irregular. Stop.'
               Stop
            endif

         elseif (TRACK_PHOTON.eq.TRACK_QE_PHOTON) THEN

            if (SK_GEOMETRY .eq. SK_IV .OR. SK_GEOMETRY .eq. SK_V) then
               call DSTRGNUM
               call DSHIGH_SK4
               if(IFANTI.ne.0) call dsanti_sk4
               call DSSWTRG
                call dsdigi_sk4
               if(IFANTI.ne.0) call dsdigi_anti_sk4

            elseif ((SK_GEOMETRY .eq. SK_I) .or. (SK_GEOMETRY .eq. SK_II) 
     &         .or. (SK_GEOMETRY .eq. SK_III) ) then

               CALL DSHIGH
               CALL DSANTI

            else
               write(6,*) 'ERROR:DS: SK_GEOMETRY is irregular. Stop.'
               Stop
            endif

         else

            write(*,*) 'DS: TRACK_PHOTON is irregal !!!'

         endif


         CALL DSPRSV
         CALL DSBDEL

      endif

      RETURN
      END
