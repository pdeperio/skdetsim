C ------------------------------
C SKDETSIM_TEST.H
C
C common for standard skdetsim diagnostic histograms
C
C

      Integer*4 TEST_events, TEST_gen_phot, TEST_coll_phot,TEST_Colla_phot
      REAL*4 TEST_Hitmap,TEST_pmtmap,TEST_Hitmapa, TEST_Nhitmap, TEST_Nhit
      integer*4 coll_phot_max, gen_phot_max, TEST_sum_coll, TEST_sum_colla
      COMMON/SK_TEST/ TEST_events, TEST_gen_phot, TEST_coll_phot
     & ,TEST_Colla_phot
      COMMON/SK_TEST2/ TEST_Hitmap(MAXPM),coll_phot_max, gen_phot_max     
     &       ,TEST_pmtmap(MAXPM),TEST_Hitmapa(MAXPMA), TEST_Nhitmap(MAXPM)
     &       , TEST_sum_coll, TEST_sum_colla, TEST_Nhit
      character*256 TEST_fname
      parameter (TEST_fname='stdtest.hbk')
      real TEST_gx, TEST_gy, TEST_gz
      integer  TEST_NOPMT, TEST_exist
      COMMON/SK_TEST_PMT/TEST_NOPMT,TEST_gx,TEST_gy
     &       , TEST_gz, TEST_exist
