import os
import numpy as np
from functools import reduce

def read_structure(path):
    os.chdir(path)
    os.system('phonopy --symmetry > phonopy_output')
    data = {}
    with open('PPOSCAR', 'r') as f:
        f = f.readlines()
        prim_vec = np.zeros((3,3))
        lines = f[2:5]
        for irow, line in enumerate(f[2:5]):
            prim_vec[irow] = np.array([float(num) for num in line.split()])
        elem = f[5].split()
        data['nelem'] = len(elem)
        elem_num = [int(n) for n in f[6].split()]
        assert len(elem) == len(elem_num)
        data['formula'] = reduce(lambda a, b: a+b, [e + str(n) if n > 1 else e for e, n in zip(elem, elem_num)])
        data['elements_num'] = elem_num
        tot_atom_num = sum(elem_num)

        atom_position = np.zeros((tot_atom_num, 3))
        for irow in range(tot_atom_num):
            atom_position[irow] = [float(num) for num in f[8 + irow].split()]
        data['poscar_str'] = f
        
    with open('phonopy_output', 'r') as f:
        f = f.readlines()
        data['nsite'] = sum(['Wyckoff' in line for line in f])
        assert 'space_group_number' in f[2]
        data['spacegroup'] =  int(f[2].strip().split()[-1])

    os.system('rm PPOSCAR BPOSCAR phonopy_output')
    return data

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

if __name__ == '__main__':
    read_symtopo_output('./test')





