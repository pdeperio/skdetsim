************************************************************************
*     ------------------
*     INCLUDE 'dsdata.h'
*     ------------------
*
*       MXPHDS ; max. # of photons
*       MXPMDS ; max. # of PMT's
*       MXRGDS ; max. # of rings
*       MXGATE ; max. # of QBEE gates in a event.(only for SK4)
* 
*       NPHADS ; # of photons
*       NPHRDS ; # of photons for individual rings
*       NRNGDS ; # of rings
*       NTRG   ; TRG count
*       TRG_OFFSET ; offset from begining of the data
*                     block to GEANT time
*
*     (Creation Date and Author)
*       1992.05.20 ; First version by K.S.Hirata
*       1995.02.03 ; Modified by Koshio
*
************************************************************************

      INTEGER  MXPHDS, MXRGDS, NRNGDS, IPARDS, NPHADS, NPHRDS, MXGATE
      Integer  ntrg, nblock_mc
      real     trg_offset

#ifdef HIGH
      PARAMETER(MXPHDS=15000000, MXRGDS=MAXRNG+1)
#else
#ifdef MIDDLE
      PARAMETER(MXPHDS=5000000, MXRGDS=MAXRNG+1)
#else
      PARAMETER(MXPHDS=5000000, MXRGDS=MAXRNG+1)
#endif
#endif

#ifdef NEUTRON
      PARAMETER(MXGATE=MAXPM*30)
#else
      PARAMETER(MXGATE=MAXPM*10)
#endif
      COMMON /DSDATA/ NRNGDS, IPARDS(MXRGDS), NPHADS, NPHRDS(MXRGDS),
     &   ntrg, trg_offset, nblock_mc
