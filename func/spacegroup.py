import spglib

def get_lattice_type(number):
    if number in [1, 2]:
        return 'Triclinic'
    if number in range(3, 16):
        return 'Monoclinic'
    if number in range(16, 75):
        return 'Orthohombic'
    if number in range(75, 143):
        return 'Tetragonal'
    if number in range(143, 168):
        return 'Trigonal'
    if number in range(168, 195):
        return 'Hexagonal'
    if number in range(195, 231):
        return 'Cubic'
    else:
        raise 'Wrong value for space group number!'


def get_spacegroup(SG_number):
    # Get spg dictionary
    for hall_number in range(1, 531):
        spg_type = spglib.get_spacegroup_type(hall_number)
        spg_number = spg_type['number']
        if spg_number == SG_number:
            schoenflies = spg_type['schoenflies']
            international = spg_type['international_short']
            pointgroup_schoenflies = spg_type['pointgroup_schoenflies']
            pointgroup_international = spg_type['pointgroup_international']
            lattice_type = get_lattice_type(SG_number)
            spg_dict = {
                    'number': SG_number,
                    'schoenflies': schoenflies,
                    'international': international,
                    'pointgroup_schoenflies': pointgroup_schoenflies,
                    'pointgroup_international': pointgroup_international,
                    'lattice_type': lattice_type
                    }

    return spg_dict


if __name__ == '__main__':
    from pprint import pprint
    pprint(spg_dict)

