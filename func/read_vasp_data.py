import os
import numpy as np
from numpy.linalg import norm
from math import acos, pi
from functools import reduce

def read_structure(path):
    os.chdir(path)
    try:
        os.system('phonopy --symmetry > phonopy_output')
    except:
        raise Exception('Fail to read POSCAR! Please install phonopy using \'pip install phonopy\'')

    data = read_POSCAR(path + '/PPOSCAR')
    data_conv = read_POSCAR(path + '/BPOSCAR')

    prim_cif_str = convert_poscar2cif(data['latt_vec'], data['atom_position'], data['elem'], data['elements_num'])
    conv_cif_str = convert_poscar2cif(data_conv['latt_vec'], data_conv['atom_position'], data_conv['elem'], data_conv['elements_num'])
    data['latt_vec_conv'] = data_conv['latt_vec']
    data['prim_cif_str'] = prim_cif_str
    data['conv_cif_str'] = conv_cif_str
        
    with open('phonopy_output', 'r') as f:
        f = f.readlines()
        data['nsite'] = sum(['Wyckoff' in line for line in f])
        assert 'space_group_number' in f[2]
        data['spacegroup'] =  int(f[2].strip().split()[-1])

   #os.system('rm PPOSCAR BPOSCAR phonopy_output')
    return data

def read_POSCAR(file_name):
    data = {}
    with open(file_name, 'r') as f:
        f = f.readlines()
        prim_vec = np.zeros((3,3))
        lines = f[2:5]
        for irow, line in enumerate(f[2:5]):
            prim_vec[irow] = np.array([float(num) for num in line.split()])
        data['latt_vec'] = prim_vec
        data['elem'] = f[5].split()
        data['nelem'] = len(data['elem'])
        data['elements_num'] = [int(n) for n in f[6].split()]
        assert len(data['elem']) == len(data['elements_num'])
        data['formula'] = reduce(lambda a, b: a+b, [e + str(n) if n > 1 else e for e, n in zip(data['elem'], data['elements_num'])])
        tot_atom_num = sum(data['elements_num'])

        atom_position = np.zeros((tot_atom_num, 3))
        for irow in range(tot_atom_num):
            atom_position[irow] = [float(num) for num in f[8 + irow].split()]
        data['atom_position'] = atom_position
        data['poscar_str'] = f
    return data

def convert_poscar2cif(latt_vec, atom_position, elem, elem_num):
    cell_length = [norm(v) for v in latt_vec]
    cell_angle = np.zeros(3)
    a, b, c = latt_vec 
    cell_angle[0] = acos((b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) / (cell_length[1] * cell_length[2])) * 180 / pi
    cell_angle[1] = acos((a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / (cell_length[0] * cell_length[2])) * 180 / pi
    cell_angle[2] = acos((b[0] * a[0] + b[1] * a[1] + b[2] * a[2]) / (cell_length[0] * cell_length[1])) * 180 / pi

    cif_str = 'data\n'
    cif_str += '_cell_length_a  %.10f\n'%(cell_length[0])
    cif_str += '_cell_length_a  %.10f\n'%(cell_length[1])
    cif_str += '_cell_length_a  %.10f\n'%(cell_length[2])
    cif_str += '_cell_angle_alpha  %.5f\n'%(cell_angle[0])
    cif_str += '_cell_angle_beta   %.5f\n'%(cell_angle[1])
    cif_str += '_cell_angle_gamma  %.5f\n\n'%(cell_angle[2])
    cif_str += '_symmetry_space_group_name_H-M    \'P 1\'\n'
    cif_str += 'loop_\n'
    cif_str += '_atom_site_label\n'
    cif_str += '_atom_site_type_symbol\n'
    cif_str += '_atom_site_fract_x\n'
    cif_str += '_atom_site_fract_y\n'
    cif_str += '_atom_site_fract_z\n'
    cif_str += '_atom_site_occupancy\n'
    
    cnt = 1
    elem_list = []
    for ele, ne in zip(elem, elem_num):
        elem_list.extend([ele] * ne)
    assert len(elem_list) == len(atom_position)
    name_list = [ele + str(i+1) for i, ele in enumerate(elem_list)]

    for name, ele, pos in zip(name_list, elem_list, atom_position):
        cif_str += '%s  %s  %2.15f  %2.15f %2.15f  1.0\n'%(name, ele, pos[0], pos[1], pos[2])

   #print(cif_str)
    return cif_str


def read_OUTCAR(path, tol_fermi=1e-6):

    os.chdir(path)
    with open('OUTCAR','r') as f:
        fr = f.readlines()
        out=[]
        nkp=[]
        valence = []
        conduction = []
        direct_gaps = []
        fsemi = False
        for cnt,line in enumerate(fr):
            if line[6:28]=='direct lattice vectors':
                prim_vec = np.array([[float(fr[cnt + 1][4:16]), float(fr[cnt + 1][16:29]), float(fr[cnt + 1][29:42])],
                                     [float(fr[cnt + 2][4:16]), float(fr[cnt + 2][16:29]), float(fr[cnt + 2][29:42])],
                                     [float(fr[cnt + 3][4:16]), float(fr[cnt + 3][16:29]), float(fr[cnt + 3][29:42])]])
            if line[94:100]=='NBANDS':
                nb = int(line[101:108])
            if 'NKPTS' in line:
                nk = int(line[29:39])
            if 'NELECT' in line:
                sline=line.split()
                NELECT = sline[2]
                N_occu = int(float(NELECT))
            if 'LSORBIT' in line:
                soc = True if 'T' in line.split()[2] else False
            if 'E-fermi' in line:
                ssline=line.split()
                Efermi = float(ssline[2])
                cnt +=1
                for nn in range(nk):
                    cnt +=3
                    num = 0
                    for mm in range(nb): 
                        cnt +=1
                        num +=1
                        if soc and num==N_occu:
                            sline = fr[cnt].split()
                            e1 = sline[1]
                            sline = fr[cnt+1].split()
                            e2 = sline[1]
                            valence.append(float(e1))
                            conduction.append(float(e2))
                            if (float(e2)-float(e1) < tol_fermi):
                                fsemi = True
                                nkp.append(nn+1)
                        elif not soc and N_occu%2==0 and num==N_occu/2:
                            sline = fr[cnt].split()
                            e1 = sline[1]
                            sline = fr[cnt+1].split()
                            e2 = sline[1]
                            valence.append(float(e1))
                            conduction.append(float(e2))
                            if (float(e2)-float(e1) < tol_fermi):
                                fsemi = True
                                nkp.append(nn+1)
        if (N_occu%2==0):
            v = valence.index(max(valence))
            c = conduction.index(min(conduction))
            gap = conduction[c]-valence[v]
            direct_gaps = [conduction[i]-valence[i] for i in range(len(valence))]
            direct_gap = np.min(direct_gaps)
        else:
            gap = 0.0
            direct_gaps = 0
            direct_gap =  0

    data = {}
    data['NELECT'] = N_occu
    data['Ef'] = Efermi
    data['soc'] = soc
    data['is_gapless'] = fsemi
    data['num_kpts'] = nk
    data['prim_vec'] = prim_vec
    data['indirect_gap'] = gap
    data['direct_gap'] = direct_gap
    return data 

def read_symtopo_output(path):
    os.chdir(path)
    data = {'topo_class': None, 'ind_group': None, 'sym_ind': None, 
            'dgn_IR': None, 'dgn_dim': None, 'dgn_IRdim': None, 'dgn_hsp': None,
            'dgn_hsl': None} 

    with open('symtopo_output', 'r') as f:
        f = f.readlines()

        for iline, line in enumerate(f):
            if 'classification' in line:
                data['topo_class'] = line.strip().split()[-1].strip('()')
                if 'Indicator' in f[iline + 1] and 'None' not in f[iline + 1]:
                    data['ind_group'] = f[iline + 1].split()[2]
                    data['sym_ind'] = [int(n) for n in f[iline + 1].split(':')[1].strip().split(',')]

                if data['topo_class'] == 'insulator':
                    data['topo_class'] = 'Triv_Ins'

                elif data['topo_class'] == 'HSLSM':
                    data['dgn_hsl'] = [s.strip().strip('\'') for s in f[iline + 1].split(':')[1].strip().split(',')]

                elif data['topo_class'] == 'HSPSM':
                    n_dgn_kpts = int(f[iline + 2].split(':')[-1])
                    dgn_hsp, dgn_IR, dgn_dim, dgn_IRdim = [], [], [], []
                    for ik in range(n_dgn_kpts):
                        line = f[iline + 4 + ik].split()
                        dgn_hsp.append(line[0])
                        dgn_IR_tmp = line[2:-1]
                        if len(dgn_IR_tmp) == 1:
                            dgn_IR_str = dgn_IR_tmp[0]
                        elif len(dgn_IR_tmp) == 3:
                            dgn_IR_str = dgn_IR_tmp[0] + '_' + dgn_IR_tmp[-1]
                        else:
                            raise ValueError('Wrong format!', line)
                        dgn_IR.append(dgn_IR_str)
                        dgn_dim.append(int(line[-1]))
                        dgn_IRdim.append('%s(%s)'%(dgn_IR_str, line[-1]))
                    data['dgn_hsp'], data['dgn_IR'], data['dgn_dim'], data['dgn_IRdim'] = dgn_hsp, dgn_IR, dgn_dim, dgn_IRdim
    return data

def read_mp_icsd_id(home_path):
    os.chdir(home_path)
    if os.path.exists('mp_icsd_id'):
        with open('mp_icsd_id', 'r') as f:
            f = f.readlines()
            assert 'mp_id:' in f[0] and 'icsd_id' in f[1], 'The two ids should be written as: mp_id: XXX'
            mp_id = f[0].strip().split(':')[1].strip()
            icsd_id = f[1].strip().split(':')[1].strip(' []').split(',')
            return mp_id, icsd_id
    else:
        return None, None


if __name__ == '__main__':
    read_symtopo_output('./test')





