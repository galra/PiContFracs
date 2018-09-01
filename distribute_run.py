from configfile import ConfigParser
from configfile import CONFIG_PARAMS_TYPES as LEGAL_CONFIG_PARAMS
from copy import deepcopy
import math
import os
import sys
from pprint import pprint

class Parameters:
    pass


def main(split_num=4, configfile='config.ini'):
    if not math.log(split_num, 2).is_integer():
        raise ValueError('Only split numbers of a power of 2 are allowed (2, 4, 8, ...)')
    split_num = int(math.log(split_num, 2))
    params = Parameters()
    config_parser = ConfigParser(configfile=configfile)
    config = config_parser.get_config()
    # Load variables from config
    for config_variable in LEGAL_CONFIG_PARAMS:
        # print('Setting up %s' % config_variable)
        setattr(params, config_variable, config[config_variable])
        print('%s is set up to %s ' % (config_variable, getattr(params, config_variable)))

    distro_dirname = '%s_dist_%i' % (os.path.splitext(configfile)[0], 2**split_num)
    if not os.path.isdir(distro_dirname):
        os.mkdir(distro_dirname)

    params_list = split_params([params.a_coeffs_range, params.b_coeffs_range], split_num)
    pprint(params_list)
    for i, (a_coeffs, b_coeffs) in enumerate(params_list):
        results_config = ConfigParser()
        results_config.add_section('Setup')
        for config_variable in LEGAL_CONFIG_PARAMS:
            if config_variable == 'a_coeffs_range':
                results_config.set('Setup', config_variable, str(a_coeffs))
            elif config_variable == 'b_coeffs_range':
                results_config.set('Setup', config_variable, str(b_coeffs))
            else:
                results_config.set('Setup', config_variable, str(getattr(params, config_variable)))
        results_config_filename = os.path.splitext(configfile)
        results_config_filename = results_config_filename[0] + '_%d' % i + results_config_filename[1]
        results_config_filename = os.path.join(distro_dirname,
                                               results_config_filename)
        with open(results_config_filename, 'w') as results_config_file:
            results_config.write(results_config_file)

def split_params(ab_coeffs, split_num):
    if split_num == 0:
        return [ab_coeffs]

    a_coeffs, b_coeffs = ab_coeffs
    ab_coeffs_1, ab_coeffs_2 = deepcopy(ab_coeffs), deepcopy(ab_coeffs)
    a_first_poly, b_first_poly = a_coeffs[0], b_coeffs[0]
    i = 0
    while i < len(a_first_poly) and a_first_poly[i][1] - a_first_poly[i][0] <= 1:
        i += 1
    if i < len(a_first_poly):
        ab_index = 0
        ab_new_bound = int((a_first_poly[i][1] + a_first_poly[i][0])/2)
    else:
        i = 0
        while i < len(b_first_poly) and b_first_poly[i][1] - b_first_poly[i][0] <= 1:
            i += 1
        if i < len(b_first_poly):
            ab_index = 1
            ab_new_bound = int((b_first_poly[i][1] + b_first_poly[i][0])/2)
        else:
            raise ValueError('Non-dividable parameters')
    # ab_coeffs_1[ab_index][0][i] = a/b_first_poly[i]
    ab_coeffs_1[ab_index][0][i][1] = ab_new_bound
    ab_coeffs_2[ab_index][0][i][0] = ab_new_bound
    return split_params(ab_coeffs_1, split_num-1) + split_params(ab_coeffs_2, split_num-1)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3
        main(sys.argv[1], sys.argv[2])
    else:
        print('%s [split_num] [config_file]' % sys.argv[0])
