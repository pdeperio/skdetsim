************************************************************************
*     ----------------
*     INCLUDE 'sktq_raw.h'
*     ----------------
*
* ------- for QB -----------
*       NQISKZ	  ; Number of ALL ID hits (not just near sw trigger)
*       IHTIFLZ(I); Hit flags for ALL ID hits (I = 1, NQISKZ)
*                     11-6   (# of TRG EVENT COUNTER - 1) * 64 (0-63)
*                      5-4   charge range (0:Small, 1:Medium, 2:Large)
*                      3-2   trig ID (0: Narrow,   1: Wide
*                                     2: Pedestal, 3: Not used)
*                      1bit  In gate (1=in gate, 0=not in gate)
*                      0bit  In 1.3usec (1=in, 0=out)
*
*       ICABIZ(I) ; Cable number for ALL ID hits (I = 1, NQISKZ)
*       ITISKZ(I) ; raw T count for ALL ID hits (I = 1, NQISKZ)
*       IQISKZ(I) ; raw Q count for ALL ID hits (I = 1, NQISKZ)
*                    15-14   charge range (0:Small, 1:Medium, 2:Large)
*                    13-12   trig ID (0: Narrow,   1: Wide
*                                     2: Pedestal, 3: Not used)
*                    11bit   In gate (1=in gate, 0=not in gate)
*                    10-0bit ADC count
*
*       TISKZ(I)  ; T (ns) for ALL ID hits (I = 1, NQISKZ)
*       QISKZ(I)  ; Q (pe) for ALL ID hits (I = 1, NQISKZ)
*
*       ---
*       ITASKZ(I) ; raw T count for ALL OD hits (I = 1, NHITAZ)
*       IQASKZ(I) ; raw Q count for ALL OD hits (I = 1, NHITAZ)
*       for SK-IV OD, NHITAZ, IHTFLZ, ICABAZ, TASKZ, QASKZ are also used.
*
*
*       common/RAWTQINFO/   information before bad channel mask
*
*
*     (Creation Date and Author)
*       1992.08.27 ; First version by K.S.Hirata
*       1996.03.30 ; Added /SKTQAZ/  -- Flanagan, J.
*	1996.08.10 ; Added comments about ibada exclusion -- Flanagan, J.
*       2007.07.03 ; modified by y.koshio/y.takeuchi
*          o added QB 
*          o extend 16*MAXPMA -> 30*MAXPMA in SKTQAZ
*          o added ITASKZ & IQASKZ in SKTQAZ
*       2007.09.18 ; modified by y.takeuchi
*          o added RAWTQINFO
************************************************************************

      INTEGER qb_od_offset

      parameter (qb_od_offset = 20000)   ! cable # offset for OD

c*** to obtain raw TQREAL info. 
      INTEGER nqisk_raw, ICABBF_RAW(30*MAXPM)
      REAL pc2pe_raw, TBUF_RAW(30*MAXPM), QBUF_RAW(30*MAXPM)
      integer NHITAZ_raw, ICABAZ_raw(30 * MAXPMA)
      real TASKZ_raw(30 * MAXPMA), QASKZ_raw(30 * MAXPMA)
      integer ITISKZ_raw(30*MAXPM),IQISKZ_raw(30*MAXPM)
      integer ITASKZ_raw(30*MAXPMA),IQASKZ_raw(30*MAXPMA)
      common/RAWTQINFO/nqisk_raw,pc2pe_raw,ICABBF_RAW,TBUF_RAW,QBUF_RAW,
     &     nhitaz_raw, icabaz_raw, taskz_raw, qaskz_raw,
     &     ITISKZ_raw, IQISKZ_raw, ITASKZ_raw, IQASKZ_raw

c*** for 1.3 micro sec timing gate
      integer SKGATE_START_COUNT, SKGATE_END_COUNT
      parameter (SKGATE_START_COUNT = -1000)  ! start of 1.3usec gate
      parameter (SKGATE_END_COUNT   = +1496)  ! end of 1.3usec gate
 
#ifdef ICHI_KILO
      REAL PC2PE
      PARAMETER (PC2PE = 2.28333)
#endif
