      INTEGER    IMTVAC, IMTWTR, NMWTR, IMTACR, NMACR

      REAL       AWTR, ZWTR, WMWTR, AWTR1, AWTR2, ZWTR1, ZWTR2, DNWTR,
     &           WMWTR1, WMWTR2, 
     &           AACR, ZACR, WMACR, AACR1, AACR2, AACR3, ZACR1, ZACR2, ZACR3, 
     &           DNACR,  WMACR1, WMACR2, WMACR3

**  vacuum

      PARAMETER (IMTVAC=16)

**  water

      DIMENSION AWTR(2),ZWTR(2),WMWTR(2)
      PARAMETER ( IMTWTR = 17 )
      PARAMETER ( AWTR1 = 1.010, AWTR2 = 15.999 )
      PARAMETER ( ZWTR1 = 1.   , ZWTR2 = 8.     )
      PARAMETER ( DNWTR = 1.00 )
      PARAMETER ( NMWTR = -2 )
      PARAMETER ( WMWTR1= 2., WMWTR2= 1. )
*
* Parameters for Standard Rock
*    Add by M.Yoshida
*
       INTEGER    IMTRCK, NMRCK
       REAL       ARCK, ZRCK, WMRCK, ARCK1, ARCK2, ZRCK1, ZRCK2, DNRCK,
     &            WMRCK1, WMRCK2
       DIMENSION ARCK(2),ZRCK(2),WMRCK(2)
 
       PARAMETER ( IMTRCK = 99 )
       PARAMETER ( ARCK1 = 15.999, ARCK2 = 28.0855 )
       PARAMETER ( ZRCK1 = 8.   , ZRCK2 = 14.     )
       PARAMETER ( DNRCK = 2.7 )
       PARAMETER ( NMRCK = -2 )
       PARAMETER ( WMRCK1= 2., WMRCK2= 1. )
 

**  acryl

      DIMENSION AACR(3),ZACR(3),WMACR(3)
      PARAMETER ( IMTACR = 18 )
      PARAMETER ( AACR1 = 1.010, AACR2 = 15.999, AACR3 = 12.00 )
      PARAMETER ( ZACR1 = 1.   , ZACR2 = 8.    , ZACR3 = 6.    )
      PARAMETER ( DNACR = 1.19 )
      PARAMETER ( NMACR = -3 )
      PARAMETER ( WMACR1= 8., WMACR2= 2., WMACR3 = 5. )
