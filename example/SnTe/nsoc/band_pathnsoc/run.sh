bsub -q q_sw_share -J top_nsoc_out200  -n 20 -cgsp 64  -b -share_size 6144 -host_stack 512 ../../../../../../Application/VASP/swvasp >vasp.out
