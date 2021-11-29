      integer nwpawc, ihpaw, nseg, nphseg, mxseg, iquest, nphtot
      integer igloo, mxpart, ipseg, npart, nevt, c_part
      real p_part, th_part, ph_part
      integer ip_part
      integer  lun_seg, seg_idh
      integer seg_filebuff(33), lname
      character*132 seg_filename
      parameter (nwpawc = 1000000)
      parameter (igloo = 2**5)  ! Binning level for igloo pixelization
      parameter (mxseg = 12288)
      parameter (mxpart = 1000)
C      parameter (seg_noseg = 0)
C      parameter (seg_both = 1)
C      parameter (seg_noph = 2)
CCS03 The functionality of flgseg has been replaced by two different flags
C     please now use OUTPUT_HBOOK =1 (SEG) and TRACK_PHOTON
      parameter (seg_idh = 1)
      parameter (lun_seg = 30)
      real phseg, thseg, pseg, tseg
      integer ndk, subevt
      common /pawc/ ihpaw(nwpawc)
      common /quest/ iquest(100)
      common /prtseg/ npart, nevt, ndk, subevt, ip_part(mxpart), p_part(5,mxpart), th_part(mxpart), ph_part(mxpart),
     &                c_part(mxpart)
      common /trkseg/ nseg, nphtot, tseg(mxseg),
     &     thseg(mxseg), phseg(mxseg), nphseg(mxseg)
      common /segmod/ seg_filename, seg_filebuff, lname
