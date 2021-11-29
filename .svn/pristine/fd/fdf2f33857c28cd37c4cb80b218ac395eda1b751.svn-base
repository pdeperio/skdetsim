C--
C     secondary particle informations
C     
C     nscndprt : Number                        of secondary particles
C     itrkscnd : Parent track number (GEANT)   of the secondary particle
C     istakscnd: Parent stack track number (GEANT) of the secondary particle
C     vtxscnd  : Vertex                        of the generated point
C     pscnd    : Momentum                      of the secondary particle
C     iprtscnd : Particle code                 of the secondary particle
C     tscnd    : Generated time                of the secondary particle
C     iprntprt : Parent particle code          of the secondary particle
C     lmecscnd : Interaction code  that produced the secondary particle
C     iprnttrk : Parent track number (MCVECT)  of the secondary particle
C     iorgprt  : Parent track PID code(MCVECT) of the secondary particle
C     iprntidx : Array index (Fortran) of this secondary particle\'s parent
C                   Negative: Index corresponding to MCVECT (NEUT) stack
C                   Positive: Index corresponding to this stack
C                   0: Unmatched (E.g. muon or hadron > 10 GeV/c)
C     nchilds  : Number of daughter particles (NGKINE in GEANT)
C     ichildidx: Array index (Fortran) in this stack of this secondary particle\'s first child
C                   0: No matched children
C     pprnt    : Momentum                      of the parent particle at interaction
C     pprntinit: Initial momentum              of the parent particle at birth
C     vtxprnt  : Vertex                        of the parent particle at birth
C     iflgscnd : User interaction code for better classification than lmecscnd
C                Currently for pions only, following convention which is 
C                hybrid of NEUT VCWORK and FSIHIST:
C                     (*100 for SKDETSIM non-NEFFECT interactions e.g. elastic SGPIEL;
C                          +0 Free Hydrogen, +1 Oxygen)
C                  -1 : ESCAPE
C                   0 : INITIAL (or unmatched parent vertex if I>1)
C                   1 :
C                   2 : 
C                   3 : ABSORPTION
C                   4 : CHARGE EXCHANGE
C                   5 : 
C                   6 : 
C                   7 : HADRON PRODUCTION 
C                   8 : QUASI-ELASTIC SCATTER
C                   9 : FORWARD (ELASTIC-LIKE) SCATTER
C
C--
C     Linking input particle stack (MCVECT) to this stack
C     
C     nchildsvc  : Number of daughter particles (NGKINE in GEANT)
C     ichildidxvc: Array index (Fortran) in this stack of this MCVECT particle\'s first child
C                   0: No matched children
C
C     1998/03/16 Y.Hayato
C     2012/01/19 P.de Perio - Add user interaction flag and momentum of parent 
C     2013/02/03 P.de Perio - Add more detailed parent information and direct parent
c                             and child indices
C--

#ifndef SECMAXRNG
      INTEGER*4  secmaxrng
      PARAMETER (secmaxrng=4000)
#define SECMAXRNG
#endif
      
      integer*4 nscndprt,itrkscnd(secmaxrng),istakscnd(secmaxrng),
     &          iprtscnd(secmaxrng),iprntprt(secmaxrng),
     &          lmecscnd(secmaxrng),iprnttrk(secmaxrng),
     &          iorgprt(secmaxrng),iflgscnd(secmaxrng),
     &          iflgscnd_temp, iprntid_temp, ichildid_temp,
     &          iprntidx(secmaxrng), nchilds(secmaxrng), ichildidx(secmaxrng)

      real*4    vtxscnd(3,secmaxrng),pscnd(3,secmaxrng),
     &          tscnd(secmaxrng),pprnt(3,secmaxrng),
     &          pprntinit(3,secmaxrng),vtxprnt(3,secmaxrng)

      common/secndprt/nscndprt,itrkscnd,istakscnd,vtxscnd
     &     ,pscnd,iprtscnd,tscnd,iprntprt,lmecscnd
     &     ,iprnttrk,iorgprt,pprnt,pprntinit, vtxprnt
     &     ,iflgscnd,iflgscnd_temp, iprntidx
     &     ,nchilds, ichildidx

      integer flag_warning_scndprt
      common /warning_scndprt/flag_warning_scndprt

C     Corresponds to size in vcwork.h (MAXVC)
      integer*4 npvcscnd, nchildsvc(100), ichildidxvc(100)

      common/secndprtvc/npvcscnd, nchildsvc, ichildidxvc
