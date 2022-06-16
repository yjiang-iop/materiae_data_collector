import os
import re
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.switch_backend('agg')

def calculate_reciprocal_lattice_vectors(a):
    a1, a2, a3 = a[0], a[1], a[2]
    norm = np.dot(a1, np.cross(a2, a3))
    b1 = np.cross(a2, a3) / norm
    b2 = np.cross(a3, a1) / norm
    b3 = np.cross(a1, a2) / norm
    return 2*np.pi*np.array([b1, b2, b3])

def read_direct_lattice_vectors():
    with open('POSCAR','r') as f:
        f = f.readlines()
        a = []
        lines = f[2:5]
        for line in lines:
            a.append([float(i) for i in line.split()])
    return a

def trans_BZ_vec(vec, BZ_latt):
    # transform a reciprocal vec to cartesian lattice
    return np.dot(vec, BZ_latt) 

def read_KPOINTS():
    # transform HSP to cartesian lattice, in order to compute correct distance
    BZ_latt = calculate_reciprocal_lattice_vectors(read_direct_lattice_vectors())
    # only for line-mode kpoint file
    with open('KPOINTS','r') as file:
        lines = file.readlines()
        assert lines[2].split()[0][0].lower() == 'l', 'KPOINTS file not line-mode!'
        assert lines[3].split()[0][0].lower() == 'r', 'KPOINTS file not rec coordinates!'
        assert lines[1].split()[0].isdigit(), lines[1]
        intersections = int(lines[1].split()[0]) # how many points between two kpoints
        
        count = 1
        kcoord_list = []
        klabel_list = []
        for cnt, line in enumerate(lines):
            if cnt <4:
                continue
            else:
                if len(line.split()) >=4:
                    line = line.split()
                    kcoord = [float(line[ii]) for ii in range(3)]
                    kcoord = trans_BZ_vec(kcoord, BZ_latt)
 
                    if line[3][0].isalpha():
                        klabel = line[3]
                    elif len(line) > 4 and line[4][0].isalpha():
                        klabel = line[4]
                    else:
                        raise Exception(line)
 
                    if count % 2 == 1 and count != 1:
                        if np.linalg.norm(kcoord - kcoord_list[-1][0]) < 1e-6:
                            assert klabel == klabel_list[-1][0], klabel
                        else:
                            assert klabel != klabel_list[-1][0], (klabel, klabel_list)
                            kcoord_list[-1].append(kcoord)
                            klabel_list[-1].append(klabel)
                    else:
                        kcoord_list.append([kcoord])
                        klabel_list.append([klabel])
 
                    count += 1
    return kcoord_list, klabel_list, intersections


def read_EIGENVAL():
    with open('EIGENVAL','r') as file:
        lines=file.readlines()
        total_line=len(lines) # total number of lines in OUTCAR

        for cnt,line in enumerate(lines):
            if cnt < 5:
                continue
            if cnt == 5:
                line = line.split()
                nb = round(float(line[0]))

                nk = round(float(line[1]))
                totband = round(float(line[2]))
                assert len(lines[6].split())==0 and len(lines[7].split())==4
                start_line = 7
                break
        
        energy = np.zeros((nk, totband))
        current_k = 0

        for cnt in range(start_line, total_line):
             if len(lines[cnt].split())==4 and current_k < nk:
                cnt += 1
                for i in range(totband):
                    tmp_line = lines[cnt + i].split()
                    assert len(tmp_line) in [2,3] #and int(tmp_line[0])==i+1 # when nb>1000, ***
                    energy[current_k][i]=float(lines[cnt+i].split()[1])
                
                current_k += 1
    return energy
           
def read_bandxmgr():
    # return:
    # kpoint_position: several point on the x axis to place klabel
    # klabel_list: list of kpoint label
    # xcoord: discrete points on the x axis
    # band_info: eigenvalue of all bands, has minus Efermi from scf-OUTCAR
    
    kcoord_list, klabel_list, intersections = read_KPOINTS()
    num_kpoints = len(klabel_list)

    with open('band.xmgr','r') as file:
        lines=file.readlines()
        
        kpoint_position = []
        num_line = 0
        for ik in range(num_kpoints):            
            kpoint_position.append(float(lines[num_line].split()[0]))
            num_line += 4
        
        num_andsymbol = sum([ line.split()[0].strip('\n') == '&' for line in lines])
        num_bands = num_andsymbol - (num_kpoints - 1) * 2        

        band_info = np.zeros((num_bands, intersections * (num_kpoints-1)))
        xcoord = np.zeros(intersections * (num_kpoints-1))
        start = num_kpoints *4 - 1
        for iband in range(num_bands):
            print('line num',start - 1)
            assert lines[start - 1].split()[0].strip('\n') == '&'

            for cnt in range(intersections * (num_kpoints-1)):
                line_num = cnt + start
                tmp = lines[line_num].split()
                
                band_info[iband, cnt] = float(tmp[1])
                if iband == 0:
                    xcoord[cnt] = float(tmp[0])

            start += intersections * (num_kpoints - 1) + 1

    return kpoint_position, klabel_list, xcoord, band_info 
            
            
def label_to_tex(klabel):
    if klabel == 'GAMMA' or klabel == 'GM':
        tex = '$\Gamma$'

    elif 'SIGMA' in klabel:
        if '_' in klabel:
            klabel = klabel.split('_')
            tex = '$\Sigma_%s$'%klabel[1]
        else:
            tex = '$\Sigma$'

    elif 'LAMBDA' in klabel:
        if '_' in klabel:
            klabel = klabel.split('_')
            tex = '$\Lambda_%s$'%klabel[1]
        else:
            tex = '$\Lambda$'

    elif 'DELTA' in klabel:
        if '_' in klabel:
            klabel = klabel.split('_')
            tex = '$\Delta_%s$'%klabel[1]
        else:
            tex = '$\Delta$'

    else:
        if '_' in klabel:
            klabel = klabel.split('_')
            tex = '$%s_%s$'%(klabel[0], klabel[1])
        else:
            tex = '$%s$'%klabel

    return tex


def klabel_to_tex(klabels):
    newlabel = []
    for label in klabels:
        if len(label) == 2:
            newlabel.append(label_to_tex(label[0]) + '|' + label_to_tex(label[1]))
        elif len(label) == 1:
            newlabel.append(label_to_tex(label[0]))
        else:
            raise Exception('Wrong klabel!',label)
    
    return newlabel


def adjust_kdistance(xcoord, xlabel):
    length = xcoord[-1]
    mv_index = []
    xcoordnew = np.zeros(len(xcoord))

    for i in range(len(xcoord) - 1):
        current_len = xcoord[i+1] - xcoord[i]
        if '|' in xlabel[i] and '|' in xlabel[i+1]:
            tol = 0.05
        elif '|' in xlabel[i] or '|' in xlabel[i+1]:
            tol = 0.04
        else:
            tol = 0.03

        if  current_len < length * tol:
            mv_index.append(i)
            move_len = length * tol - current_len

            if i-1 not in mv_index and i != 0:
                xcoordnew[i+1] = xcoord[i+1] + move_len / 2
                xcoordnew[i]  = xcoord[i] - move_len /2
            
            elif i-1 in mv_index and  i!= 0:
                move_len1 = length * tol - (xcoord[i] - xcoord[i-1])
                xcoordnew[i-1] = xcoord[i-1] - move_len1
                xcoordnew[i] = xcoord[i]
                xcoordnew[i+1] = xcoord[i+1] + move_len

            else: # i == 0
                xcoordnew[i+1] = xcoord[i+1] + move_len
        else:
            xcoordnew[i+1] = xcoord[i+1]

    return xcoordnew


def plot_helper_settings(ylim, reciprocal_point_locations, reciprocal_point_labels, vertical_line_locations):
    plt.xlim(reciprocal_point_locations[0], reciprocal_point_locations[-1])
    ax = plt.gca()
    ax.xaxis.set_ticks(reciprocal_point_locations)
    
    reciprocal_point_labels = klabel_to_tex(reciprocal_point_labels)

    ax.xaxis.set_ticklabels(reciprocal_point_labels)
    ax.tick_params(tick1On = False)
    plt.axhline(0, ls='--', c='k', alpha=0.5)
    if ylim:
        plt.ylim(ylim)

    for kp_end_point in range(len(reciprocal_point_locations)):
        plt.axvline(vertical_line_locations[kp_end_point], ls='--', c='k', alpha=0.5)
    plt.ylabel('$E - E_f$ (eV)')
    plt.xlabel('Wave Vector')
   #mat_info = reduce(lambda x,y: str(x) + '-' + str(y), mat_info[:-4])
   #plt.xlabel(mat_info)

    ymajorLocator = MultipleLocator(1)
    ax.yaxis.set_major_locator(ymajorLocator)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
  #     plt.legend(loc=0, fontsize='small')
    fig = plt.gcf()
    fig.set_tight_layout(True)
    plt.draw()


def plot_band(band_path, E_fermi, nb, soc=True):
    os.chdir(band_path)
    kpoint_position, klabel_list, intersections = read_KPOINTS()
    num_kpoints = len(klabel_list)
    klabel_list_tex = klabel_to_tex(klabel_list)
    
    # read energy from EIGENVAL
    band_info = read_EIGENVAL()
    
    # set plot window
    band_num = nb if soc else nb // 2
    upper_gap_energy = np.min(band_info[:, band_num])
    lower_gap_energy = np.max(band_info[:, band_num - 1])
    plot_window = [-4, 4]
    if upper_gap_energy - E_fermi > 2:
        plot_window[1] += int(upper_gap_energy - E_fermi)
    if lower_gap_energy - E_fermi < -2:
        plot_window[0] += int(lower_gap_energy - E_fermi)

    # move fermi level to 0
    band_info = band_info - E_fermi
    band_info = band_info.T 
    xcoord = []
    kdistance = [0]
    for ik in range(num_kpoints - 1):
        dis = np.linalg.norm(np.array(kpoint_position[ik][-1]) - np.array(kpoint_position[ik+1][0]))

        x = np.linspace(kdistance[-1], kdistance[-1] + dis, intersections)
        xcoord = xcoord + x.tolist()
        kdistance.append(kdistance[-1] + dis)

    xcoord = np.array(xcoord)
    kdistance_adjusted = adjust_kdistance(kdistance, klabel_list_tex)
    
    # determine the discontinue point on bandplot
    continue_list = [[]]
    for ik in range(num_kpoints):
        if len(klabel_list[ik]) == 1:
            continue_list[-1].append(ik)
        else:
            continue_list[-1].append(ik)
            continue_list.append([ik])
    
    plt.figure()
    for iband in range(band_info.shape[0]):
        if iband < band_num:
            band_color = 'green'
        else:
            band_color = 'red'

        for i,ilist in enumerate(continue_list):
            plt.plot(xcoord[intersections*ilist[0]:intersections*ilist[-1]], band_info[iband, intersections*ilist[0]:intersections*ilist[-1]], \
                    color = band_color, linewidth = 1)
    
    plot_helper_settings(plot_window, kdistance_adjusted, klabel_list, kdistance)
    plt.savefig('bandplot.svg',format = 'svg')
    plt.close()



if __name__ == '__main__':
    
    soc = True 
    soc_str = 'soc' if soc else 'nsoc'
    band_path = '.'
    os.chdir(path)

    NELECT = int(float(os.popen('grep NELECT OUTCAR').read().split()[2]))
    nb = NELECT if soc else round(NELECT/2)
    E_fermi = float(os.popen('grep E-fermi ../scf_%s/OUTCAR'%soc_str).read().split()[2])

    plot_band(band_path, E_fermi, NELECT, soc)

