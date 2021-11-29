c
c    Count # of photo-electrons in gates
c    Multi-hit is taken into account.
c
    
c --------- Inner Detector variables
       integer pmt_id(MXGATE), nphoto_gate(MXGATE)
       real timing_gate(MXGATE),width_photons(MXGATE)
       integer num_gate(MAXPM)
#ifdef NEUTRON
       integer list_gate(MAXPM,500)
#else
       integer list_gate(MAXPM,50)
#endif

       integer num_gate_all

       common /DSMULTI/pmt_id,nphoto_gate,timing_gate,width_photons,
     & list_gate,num_gate,num_gate_all

c --------- outer Detector variables
       integer pmt_id_anti(mxgate)
       integer nphoto_gate_anti(mxgate)
       real    timing_gate_anti(mxgate)
       real    width_photons_anti(mxgate)
#ifdef NEUTRON
       integer list_gate_anti(maxpma,500)
#else
       integer list_gate_anti(maxpma,50)
#endif
       integer num_gate_anti(maxpma)
       integer num_gate_all_anti

       common /dsmultia/pmt_id_anti,nphoto_gate_anti,timing_gate_anti,width_photons_anti,
     & list_gate_anti,num_gate_anti,num_gate_all_anti

