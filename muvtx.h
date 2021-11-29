C A bank of muon vertices used by the upmu analysis
C Meaning:
C Word    Vertex
C ----------------------------------
C  1      Inner Entrance Point
C  2      Inner Exit Point
C  3      Stopping Point
C  4      OD Entrance Point
C  5      OD Exit Point
C
C Meaning of idmu variable:
C   0001 digit - has entered OD already?
C   0010 digit - has entered ID already?
C   0100 digit - has exited ID already?
C   1000 digit - has exited OD already?
C
C Stored in MUVTX bank
C Only last muon in the initial set of particles saved in gukine gets
C this treatment.  Yes, this might not be the muon we care about.
C But, it is most of the time and this is simple.  first_id_mu contains the
C ITRA number of that muon, that is what is used to tag it.
 
      integer nvtx
      parameter(nvtx=5)
      integer idmu
      integer first_id_mu
      integer last_medium
      real mupos(nvtx,3),mudir(nvtx,3),mumom(nvtx)
      common /muvtx/ mupos,mudir,mumom
      common /muvtx2/ idmu,first_id_mu,last_medium
