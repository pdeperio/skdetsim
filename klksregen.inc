c -*- FORTRAN -*-
      
      REAL KLKSLIFE_H20,KLKSLIF,KLKSREGL
      LOGICAL KLKSENABLE 

c     KL->KS regeneration is treated as a "decay" with a mean lifetime
c     given by this parameter in seconds. (see README_KLKSregen for more detail)
c      PARAMETER (KLKSLIFE_H20 = 3.E-7) !~5% prob (ad hoc)
      PARAMETER (KLKSLIFE_H20 = 120.E-7) ! ~0.1% prob (from reference Nuclear Instruments and Methods in Physics Research A 350 (1994) 144-149)


c     set true to allow kl->ks regeneration
      PARAMETER (KLKSENABLE = .FALSE.)


      COMMON/KLKSREG/KLKSLIF,KLKSREGL
