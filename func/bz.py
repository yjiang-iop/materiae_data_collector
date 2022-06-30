import numpy
from func import brillouinzone
from pprint import pprint

all_kpoints = {}
for i in range(1, 231):
    path = './func/kpoints/%d' % i
    with open(path) as file:
        all_kpoints[i] = {}
        for line in file.readlines()[3:]:
            line = line.split()
            k_vec = [float(s) for s in line[:3]]
            k_label = line[4]
            all_kpoints[i][k_label] = k_vec

def get_BZ(prim_lattice, spg_num):
    # Calculate reciprocal lattice from primitive lattice
    rec_lattice = numpy.array(get_reciprocal_cell_rows(prim_lattice))
    b1, b2, b3 = rec_lattice
    a1, a2, a3 = prim_lattice
    # Calculate BZ faces from reciprocal lattice
    faces_data = brillouinzone.get_BZ(b1=b1, b2=b2, b3=b3)
    # Get high symmetry points and lines from database
    points_rel, lines = get_HSP_HSL(spg_num)
    # Convert kpoints to absolute values
    points_abs = rel_to_abs(points_rel, rec_lattice) 
    BZ_data = {'b1': b1.tolist(),
            'b2': b2.tolist(),
            'b3': b3.tolist(),
            'a1': a1,
            'a2': a2,
            'a3': a3,
            'kpoints': points_abs,
            'faces_data': faces_data,
            'path': lines}
    #print('rec_lattice', rec_lattice)
    #print('points_rel', points_rel)
    #print('points_abs', points_abs)
    return BZ_data

def get_reciprocal_cell_rows(real_space_cell):
    r"""
    Given the cell in real space (3x3 matrix, vectors as rows, 
    return the reciprocal-space cell where again the G vectors are 
    rows, i.e. satisfying 
    ``dot(real_space_cell, reciprocal_space_cell.T)`` = :math:`2 \pi I`,
    where :math:`I` is the :math:`3\times 3` identity matrix.

    :return: the :math:`3\times 3` list of reciprocal lattice vectors where each row is 
        one vector.
    """
    reciprocal_space_columns = 2. * numpy.pi * \
        numpy.linalg.inv(real_space_cell)
    return (reciprocal_space_columns.T).tolist()

def rel_to_abs(points_rel, rec_lattice):
    b1, b2, b3 = rec_lattice
    # Bring all points to the first BZ
    points_rel0 = {}
    def gen_coor(x):
        if x == 0:
            return [x]
        if x > 0:
            return [x-1, x]
        else:
            return [x, x+1]

    def align(p):
        if p[0] > 0:
            return p
        elif p[0] < 0:
            return -p
        if p[1] > 0:
            return p
        elif p[1] < 0:
            return -p
        if p[2] < 0:
            return -p
        return p

    for label, p in points_rel.items():
        p = numpy.array(p)
        c1s = gen_coor(p[0])
        c2s = gen_coor(p[1])
        c3s = gen_coor(p[2])
        ps_rel = numpy.array([[c1, c2, c3] for c1 in c1s for c2 in c2s for c3 in c3s])
        M = numpy.array([b1, b2, b3])
        ps_abs = numpy.dot(ps_rel, M)
        min_i = numpy.argmin([numpy.dot(p, p) for p in ps_abs])
        points_rel0[label] = align(ps_rel[min_i])

    ## Convert to absolute
    points_abs = {
        k: (v[0] * b1 + v[1] * b2 + v[2] * b3).tolist()
        for k, v in points_rel0.items()
    }
    return points_abs

def get_HSP_HSL(spg_num):
    #TODO: investigate has_inversion_symmetry and time_revsersal_symmetry
    points_rel = all_kpoints[spg_num]
    lines = []
    return points_rel, lines

if __name__ == '__main__':
    pprint(all_kpoints)

