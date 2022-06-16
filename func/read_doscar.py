import numpy as np
import os

def read_DOSCAR(path):
    os.chdir(path)
    with open('DOSCAR', 'r') as doscar:
        dosline = doscar.readlines()
        num_dospoint = int(float(dosline[5].split()[2])) # number of points in doscar
        energy = np.zeros(num_dospoint)
        dos = np.zeros(num_dospoint)
        totdos = np.zeros(num_dospoint)
 
        for cnt, line in enumerate(dosline):
            if cnt==6:
                assert len(line.split())==3
                for ii in range(num_dospoint):
                    line = dosline[cnt + ii].split()
                    assert len(line)==3
                    energy[ii] = float(line[0])
                    
                    # read dos
                    dos_tmp = line[1].split('E')
                    if dos_tmp[1][0] == '+':
                        dos[ii] = float(dos_tmp[0]) * 10**(float(dos_tmp[1][1:]))
                    elif dos_tmp[1][0] == '-':
                        dos[ii] = float(dos_tmp[0]) * 10**(-1 * float(dos_tmp[1][1:]))
                    
                    # read totdos
                    dos_tmp = line[2].split('E')
                    if dos_tmp[1][0] == '+':
                        totdos[ii] = float(dos_tmp[0]) * 10**(float(dos_tmp[1][1:]))
                    elif dos_tmp[1][0] == '-':
                        totdos[ii] = float(dos_tmp[0]) * 10**(-1 * float(dos_tmp[1][1:]))
                break
       #print(energy)
       #print(dos)
       #print(totdos)
    return energy, dos, totdos 

def save_dos_data(energy, dos, totdos, Ef_scf, sg, mat_formula, save_path):
    assert energy[0] < Ef_scf and energy[-1]> Ef_scf, (energy, Ef_scf) 
    energy = np.array(energy) - Ef_scf
    
    num_dos_points = len(energy)
    with open(save_path + '/dos.dat', 'w') as doscar:
        doscar.write('%s   %s   mp-num    N_points: %s \n\n'%(sg, mat_formula, num_dos_points))
        doscar.write('%15s    %15s    %15s \n\n'%('energy', 'dos', 'total_dos'))
        for ii in range(len(energy)):
            doscar.write('%15.7f     %15.7f    %15.7f     \n'%(energy[ii], dos[ii], totdos[ii]))
    
    dos_data = np.array([[e, d] for e, d in zip(energy, dos)])
    return dos_data

def find_fermi_dos(energy, dos, totdos, NELECT, delta=1e-4, error=1e-3):
    # find the dos value at fermi level, as well as the more accuracy fermi energy from doscar

   #energy, dos, totdos = read_DOSCAR()
    delta_n = NELECT * delta
    lower_n = NELECT - delta_n
    upper_n = NELECT + delta_n
    num_dospoint = len(dos)
    
    for i in range(num_dospoint):
        if totdos[i] <= lower_n and totdos[i+1] > lower_n:
            ratio1 = (lower_n - totdos[i]) / (totdos[i+1] - totdos[i])
            lower_energy = energy[i] + ratio1 * (energy[i+1] - energy[i])

        if totdos[i] <= upper_n and totdos[i+1] > upper_n:
            ratio2 = (upper_n - totdos[i]) / (totdos[i+1] - totdos[i])
            upper_energy = energy[i] + ratio2 * (energy[i+1] - energy[i])
            break

    E_fermi = (upper_energy + lower_energy) / 2

    for i in range(num_dospoint):
        if energy[i] <= E_fermi and energy[i+1] > E_fermi:
            fermi_index = i
            ratio = (E_fermi - energy[i]) / (energy[i+1] - energy[i])
            fermi_dos = dos[i] + ratio * (dos[i+1] - dos[i])
            fermi_dos = fermi_dos if fermi_dos > error else 0
            break
    
    # find gap 
    if fermi_dos < error:
        i = fermi_index
        while dos[i] < error:
            i -= 1
        lower_index = i
        if lower_index == fermi_index:
            ratio = (error - fermi_dos) / (dos[lower_index] - fermi_dos)
            lower_gap = ratio * (E_fermi - energy[lower_index])
        elif lower_index < fermi_index:
            ratio = (error - dos[lower_index + 1]) / (dos[lower_index] - dos[lower_index + 1])
            lower_gap = ratio * (energy[lower_index + 1] - energy[lower_index]) + (E_fermi - energy[lower_index + 1])
        else:
            raise Exception('Wrong lower_index!', lower_index,fermi_index)

        i = fermi_index + 1
        while dos[i] < error:
            i += 1
        
        upper_index = i
        if upper_index == fermi_index + 1:
            ratio = (error - fermi_dos) / (dos[upper_index] - fermi_dos)
            upper_gap = ratio * (energy[upper_index] - E_fermi)
        elif upper_index > fermi_index + 1:
            ratio = (error - dos[upper_index - 1]) / (dos[upper_index] - dos[upper_index - 1])
            upper_gap = ratio * (energy[upper_index] - energy[upper_index - 1]) + (energy[upper_index - 1] - E_fermi)
        else:
            raise Exception('Wrong upper_index!', upper_index,fermi_index)       

        gap = lower_gap + upper_gap
    else:
        gap = 0

    return E_fermi, fermi_dos, gap
   