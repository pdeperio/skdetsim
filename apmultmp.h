C-Secondary PARTICLE STORAGE COMMON
      INTEGER*4 MAX2NDP
      PARAMETER (MAX2NDP=MAXRNG*40)
      real vpi0,ppi0,tpi0
      integer npi0, ippi0, itrcon, iprtpi0, imeca
      common/multmp/npi0,ippi0(MAX2NDP),vpi0(3,MAX2NDP)
     &     ,ppi0(3,MAX2NDP),itrcon(MAX2NDP),iprtpi0(MAX2NDP)
     &     ,tpi0(MAX2NDP),imeca(MAX2NDP)

      integer flag_warning_convvect
      common /warning_convvect/flag_warning_convvect
