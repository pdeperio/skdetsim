#ifndef SECMAXRNG
C      INTEGER*4  secmaxrng
      INTEGER  secmaxrng
      PARAMETER (secmaxrng=4000)
#define SECMAXRNG
#endif
      
C      integer*4 n2prt,itrk(secmaxrng),
      integer n2prt,itrk(secmaxrng),
     &          iprt2nd(secmaxrng),
     &          iprtorg(secmaxrng),lmec(secmaxrng)
C      real*4    vtx2nd(3,secmaxrng),p2nd(3,secmaxrng),
      real    vtx2nd(3,secmaxrng),p2nd(3,secmaxrng),
     &          t2nd(secmaxrng)

      common/secndprt/n2prt,itrk,vtx2nd
     &     ,p2nd,iprt2nd,t2nd
     &     ,iprtorg,lmec


C      INTEGER*4  maxnsubev
      INTEGER  maxnsubev
      PARAMETER (maxnsubev=10)
C      INTEGER*4  maxnsubnev
      INTEGER  maxnsubnev
      PARAMETER (maxnsubnev=100)
      Integer numsubev
      Integer nsublist(maxnsubev,maxnsubnev), nevinsube(maxnsubev)
      Real nevstartt(maxnsubev)
      Common/subevent/numsubev, nsublist, nevinsube, nevstartt

      Integer tingate
      Data tingate/950./

