import os
import numpy as np
from func.read_vasp_data import *
from func.read_doscar import *
from func.plot_band import plot_band, read_EIGENVAL

'''
1. Assume the material is computed using symtopo, with folders scf/, wave/, band/, dos/.
   In the wave/ folder, the symtopo output is saved in a file symtopo_output.
2. collect basic info of the material.
   Info that needs to be input manually: mp_id, icsd_ids

'''


def collect_data(home_path):
    # this script should be exculated in the folder that contains the following folders: scf/, wave/, band/, and dos/  
    os.chdir(home_path)
    assert 'scf' in os.listdir() and 'band' in os.listdir(), (
        'Please exculated this script in the folder that contains the following folders: scf/, wave/, band/, and dos/')
    scf_path, wave_path,  = home_path + '/scf', home_path + '/wave'

    data = {'formula': None, 'elements_num': None, 'nelem': None, 'mp_id': None, 'icsd_ids': None,
            'nelec': None, 'nsite': None, 'spacegroup': None, 
            'prim_cif_str': None, 'conv_cif_str': None, 'BZ': None, 
            'efermi': None, 'dos': None, 'fermi_dos': None, 'dos_gap': None, 
            'indirect_gap': None, 'band': None, 
            'topo_class': None, 'ind_group': None, 'sym_ind': None, 
            'dgn_IR': None, 'dgn_dim': None, 'dgn_irrep': None, 'dgn_hsp': None,
            'dgn_hsl': None}
    
    # read info from POSCAR
    poscar_data = read_structure(wave_path)
    data['nelem'] = poscar_data['nelem']
    data['formula'] = poscar_data['formula']
    data['elements_num'] = poscar_data['elements_num']
    data['nsite'] = poscar_data['nsite']
    data['spacegroup'] = poscar_data['spacegroup']
    # FIXME: not sure what is 'prim_cif_str', 'conv_cif_str', check database

    # read info from OUTCAR
    os.chdir(home_path)
    outcar_data = read_OUTCAR(scf_path)
    Ef_scf = outcar_data['Ef']
    soc = outcar_data['soc']
    data['nelec'] = outcar_data['NELECT']
    data['efermi'] = Ef_scf
    data['indirect_gap'] = outcar_data['indirect_gap']

    # read info from wave, i.e., topological classification data
    os.chdir(home_path)
    assert os.path.exists('%s/symtopo_output'%wave_path), ('Please save the output of symtopo to a file named symtopo_output.')
    symtopo_data = read_symtopo_output(wave_path)
    data['topo_class'], data['ind_group'], data['sym_ind'] = symtopo_data['topo_class'], symtopo_data['ind_group'], symtopo_data['sym_ind']
    data['dgn_IR'], data['dgn_dim'], data['dgn_irrep'] = symtopo_data['dgn_IR'], symtopo_data['dgn_dim'], symtopo_data['dgn_irrep']
    data['dgn_hsp'], data['dgn_hsl'] = symtopo_data['dgn_hsp'], symtopo_data['dgn_hsl']
    # FIXME: not sure what is 'dng_IR' and 'dgn_irrep', check database

    # read info from DOSCAR
    os.chdir(home_path)
    if os.path.exists('dos'):
        dos_path = home_path + '/dos'
        dos_energy, dos, totdos = read_DOSCAR(dos_path)
        save_dos_data(dos_energy, dos, totdos, Ef_scf, data['spacegroup'], data['formula'], home_path)
        Ef_dos, fermi_dos, dos_gap = find_fermi_dos(dos_energy, dos, totdos, data['nelec'])
        data['fermi_dos'] = fermi_dos
        data['dos_gap'] = dos_gap
        data['efermi'] = Ef_dos
        assert abs(Ef_scf - Ef_dos) < 2, ('Large difference of Ef from scf and dos! Results may be wrong.', Ef_scf, Ef_dos)
       #data['dos'] = dos  # FIXME: not sure what is data['dos'], check database

    # read info from band
    os.chdir(home_path)
    if os.path.exists('band'):
        band_path = home_path + '/band'
        os.chdir(band_path)
        band_energy = read_EIGENVAL()
        plot_band(band_path, data['efermi'], data['nelec'], soc)

    os.chdir(home_path)
    print(data)
    np.save('mat.npy', data)
    return data



if __name__ == '__main__':
    path = os.getcwd() + '/example/SnTe'
    collect_data(path)



