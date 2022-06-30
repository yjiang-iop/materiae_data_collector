import os
import numpy as np
from func.read_vasp_data import *
from func.read_doscar import *
from func.plot_band import plot_band, read_EIGENVAL
from func.bz import get_BZ

'''
1. Assume the material is computed using vasp, with results in folders scf_(n)soc/, wave_(n)soc/, band_(n)soc/, dos_(n)soc/.
   If only part of the computation is done, eg, only soc data, this script will collect data from existing folder.

2. In the wave_(n)soc/ folder, the wavefunction at high-symmetry points is computed, 
   and the output of symtopo is saved in a file named symtopo_output.

3. This script is used to collect basic info of the material, generating output files mat_data.npy and bandplot.svg.
   Info that needs to be input manually: mp_id, icsd_ids.

4. Files that used: scf_(n)soc/OUTCAR, wave_(n)soc/POSCAR, symtopo_output, 
   band_(n)soc/EIGENVAL, POSCAR, KPOINTS, dos_(n)soc/DOSCAR

'''


def collect_data(home_path):
    # this script should be exculated in the folder that contains the following folders:   
    # if the folders have different names, one can change the following dir_name: 
    dir_name = {'scf_soc': 'scf_soc', 'wave_soc': 'wave_soc', 'band_soc': 'band_soc', 'dos_soc': 'dos_soc',
                'scf_nsoc': 'scf_nsoc', 'wave_nsoc': 'wave_nsoc', 'band_nsoc': 'band_nsoc', 'dos_nsoc': 'dos_nsoc'}

    os.chdir(home_path)
    assert dir_name['scf_soc'] in os.listdir() and dir_name['wave_soc'] in os.listdir(), (
        'Please exculated this script in the folder that contains the following folders: scf/, wave/, band/, and dos/')
    scf_soc_path, wave_soc_path  = home_path +'/'+ dir_name['scf_soc'], home_path +'/'+ dir_name['wave_soc']
    band_soc_path, dos_soc_path  = home_path +'/'+ dir_name['band_soc'], home_path +'/'+ dir_name['dos_soc']
    scf_nsoc_path, wave_nsoc_path  = home_path +'/'+ dir_name['scf_nsoc'], home_path +'/'+ dir_name['wave_nsoc']
    band_nsoc_path, dos_nsoc_path  = home_path +'/'+ dir_name['band_nsoc'], home_path +'/'+ dir_name['dos_nsoc']

    data = {'formula': None, 'elements_num': None, 'nelem': None, 'mp_id': None, 'icsd_ids': None,
            'nelec': None, 'nsite': None, 'spacegroup': None,  
            'prim_cif_str': None, 'conv_cif_str': None, 'BZ': None, 
            'soc_efermi': None, 'soc_dos': None, 'soc_fermi_dos': None, 'soc_dos_gap': None, 
            'soc_indirect_gap': None, 'soc_band': None, 
            'soc_topo_class': None, 'soc_ind_group': None, 'soc_sym_ind': None, 
            'soc_dgn_IR': None, 'soc_dgn_dim': None, 'soc_dgn_irrep': None, 'soc_dgn_hsp': None, 'soc_dgn_hsl': None,
            'nsoc_efermi': None, 'nsoc_dos': None, 'nsoc_fermi_dos': None, 'nsoc_dos_gap': None, 
            'nsoc_indirect_gap': None, 'nsoc_band': None, 
            'nsoc_topo_class': None, 'nsoc_ind_group': None, 'nsoc_sym_ind': None, 
            'nsoc_dgn_IR': None, 'nsoc_dgn_dim': None, 'nsoc_dgn_irrep': None, 'nsoc_dgn_hsp': None, 'nsoc_dgn_hsl': None}
    
    # read info from wave_soc/POSCAR
    poscar_data = read_structure(scf_soc_path)
    data['nelem'] = poscar_data['nelem']
    data['formula'] = poscar_data['formula']
    data['elements_num'] = poscar_data['elements_num']
    data['nsite'] = poscar_data['nsite']
    data['spacegroup'] = poscar_data['spacegroup']
   #data['poscar_str'] = poscar_data['poscar_str'] 
    data['prim_cif_str'] = poscar_data['prim_cif_str']
    data['conv_cif_str'] = poscar_data['conv_cif_str']
    data['BZ'] = get_BZ(poscar_data['latt_vec'], poscar_data['spacegroup'])

    # read info from scf_soc/OUTCAR
    os.chdir(home_path)
    if os.path.exists(scf_soc_path):
        outcar_data = read_OUTCAR(scf_soc_path)
        Ef_scf_soc = outcar_data['Ef']
        data['nelec'] = outcar_data['NELECT']
        data['soc_efermi'] = Ef_scf_soc
        data['soc_indirect_gap'] = outcar_data['indirect_gap']

    os.chdir(home_path)
    if os.path.exists(scf_nsoc_path):
        outcar_data = read_OUTCAR(scf_nsoc_path)
        Ef_scf_nsoc = outcar_data['Ef']
        data['nelec'] = outcar_data['NELECT']
        data['nsoc_efermi'] = Ef_scf_nsoc
        data['nsoc_indirect_gap'] = outcar_data['indirect_gap']

    # read info from wave_(n)soc/symtopo_output, i.e., topological classification data
    os.chdir(home_path)
    if os.path.exists(wave_soc_path):
        assert os.path.exists('%s/symtopo_output'%wave_soc_path), ('Please save the output of symtopo to a file named symtopo_output.')
        symtopo_data = read_symtopo_output(wave_soc_path)
        data['soc_topo_class'], data['soc_ind_group'], data['soc_sym_ind'] = symtopo_data['topo_class'], symtopo_data['ind_group'], symtopo_data['sym_ind']
        data['soc_dgn_IR'], data['soc_dgn_dim'], data['soc_dgn_irrep'] = symtopo_data['dgn_IR'], symtopo_data['dgn_dim'], symtopo_data['dgn_IRdim']
        data['soc_dgn_hsp'], data['soc_dgn_hsl'] = symtopo_data['dgn_hsp'], symtopo_data['dgn_hsl']

    os.chdir(home_path)
    if os.path.exists(wave_nsoc_path):
        assert os.path.exists('%s/symtopo_output'%wave_nsoc_path), ('Please save the output of symtopo to a file named symtopo_output.')
        symtopo_data = read_symtopo_output(wave_nsoc_path)
        data['nsoc_topo_class'], data['nsoc_ind_group'], data['nsoc_sym_ind'] = symtopo_data['topo_class'], symtopo_data['ind_group'], symtopo_data['sym_ind']
        data['nsoc_dgn_IR'], data['nsoc_dgn_dim'], data['nsoc_dgn_irrep'] = symtopo_data['dgn_IR'], symtopo_data['dgn_dim'], symtopo_data['dgn_IRdim']
        data['nsoc_dgn_hsp'], data['nsoc_dgn_hsl'] = symtopo_data['dgn_hsp'], symtopo_data['dgn_hsl']

    # read info from dos_(n)soc/DOSCAR
    os.chdir(home_path)
    if os.path.exists(dos_soc_path):
        dos_energy, dos, totdos = read_DOSCAR(dos_soc_path)
        dos_data = save_dos_data(dos_energy, dos, totdos, Ef_scf_soc, data['spacegroup'], data['formula'], dos_soc_path)
        Ef_dos, fermi_dos, dos_gap = find_fermi_dos(dos_energy, dos, totdos, data['nelec'])
        data['soc_fermi_dos'] = fermi_dos
        data['soc_dos_gap'] = dos_gap
        data['soc_efermi'] = Ef_dos
        data['soc_dos'] = dos_data

    os.chdir(home_path)
    if os.path.exists(dos_nsoc_path):
        dos_energy, dos, totdos = read_DOSCAR(dos_nsoc_path)
        dos_data = save_dos_data(dos_energy, dos, totdos, Ef_scf_nsoc, data['spacegroup'], data['formula'], dos_nsoc_path)
        Ef_dos, fermi_dos, dos_gap = find_fermi_dos(dos_energy, dos, totdos, data['nelec'])
        data['nsoc_fermi_dos'] = fermi_dos
        data['nsoc_dos_gap'] = dos_gap
        data['nsoc_efermi'] = Ef_dos
        data['nsoc_dos'] = dos_data

    # read info from band_(n)soc/EIGENVAL, POSCAR, KPOINTS
    os.chdir(home_path)
    if os.path.exists(band_soc_path):
        os.chdir(band_soc_path)
        band_energy = read_EIGENVAL()
        plot_band(band_soc_path, data['soc_efermi'], data['nelec'], soc=True)
        data['soc_band'] = True

    os.chdir(home_path)
    if os.path.exists(band_nsoc_path):
        os.chdir(band_nsoc_path)
        band_energy = read_EIGENVAL()
        plot_band(band_nsoc_path, data['nsoc_efermi'], data['nelec'], soc=False)
        data['nsoc_band'] = True

    os.chdir(home_path)
    for key, value in data.items():
        print(key, ':', value)
    np.save('mat_data.npy', data)
    return data


def load_data(path):
    data = np.load('%s/mat_data.npy'%path, allow_pickle=True)


if __name__ == '__main__':
    path = os.getcwd() + '/example/SnTe'
    collect_data(path)

   #data = load_data(path)



