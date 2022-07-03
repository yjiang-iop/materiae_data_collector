# materiae_data_collector
A script for collecting material data for the materiae datebase.

- Assume the material is computed using vasp, with results in folders scf_(n)soc/, wave_(n)soc/, band_(n)soc/, dos_(n)soc/.
   If only part of the computation is done, eg, only soc data, this script will only collect data from existing folders.

   In the wave_(n)soc/ folder, the wavefunction at high-symmetry points is computed, 
   and the output of symtopo should be saved in a file named symtopo_output.

- This script is used to collect basic info of the material, generating output files mat_data.npy and bandplot.svg,
   stored in folder 'materiae_data'.

   Files that used: scf_(n)soc/OUTCAR, POSCAR, wave_(n)soc/symtopo_output, 
   band_(n)soc/EIGENVAL, POSCAR, KPOINTS, dos_(n)soc/DOSCAR

   Required package: phonopy, spglib (To installed: pip install phonopy /spglib)

- If one has the mp_id and icsd_ids, i.e., the id of materials project website and ICSD database, 
   one can store them in a file named 'mp_icsd_id'.