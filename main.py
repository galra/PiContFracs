#!/usr/bin/env python3
"""Integrates the different modules to run a search for a continued fraction that converges to a target constant."""
import enum_params
import time
import datetime
from gen_consts import gen_pi_const, gen_e_const, gen_feig_consts, gen_euler_masch_const, gen_percolation_consts, \
                            gen_zeta_consts
from postprocfuncs import POSTPROC_FUNCS, EVALUATED_POSTPROC_FUNCS
import dill as pickle
from latex import generate_latex
import os
from configfile import ConfigParser
from configfile import CONFIG_PARAMS_TYPES as LEGAL_CONFIG_PARAMS
import sys
import pdb
from utils import MeasureRuntime


class Parameters(dict):
    """Stub class for parameters, parsed parameters can be defined as instance's members."""
    def __getattr__(self, item):
        try:
            return self.__getitem__(item)
        except KeyError as exc:
            raise AttributeError(item) from exc

    def __setattr__(self, key, value):
        self.__setitem__(key, value)


# example - this is for testing a specific set of parameters
# rhs_polys_enumer = enum_params.BasicEnumPolyParams()
# rhs_polys_enumer.polys_generator(range_a=[[6,7], [-1,0], [1,2]], range_b=[[5,6],[4,5],[-3,-2]])

def main(configfile='config.ini'):
    """Supported consts: pi, e, feig(0-3), euler_masch, percolation (0-1), zeta (2-20).
    For feig, i=0,1,2,3 is required.
    For percolation, i=0,1 is required.
    For zeta, i=2,3,4,...,19,20 is required."""
    # Load configuration file
    params = Parameters()
    config_parser = ConfigParser(configfile=configfile)
    config = config_parser.get_config()
    # Load variables from config
    for config_variable in LEGAL_CONFIG_PARAMS:
        # print('Setting up %s' % config_variable)
        setattr(params, config_variable, config[config_variable])
        print('%s is set up to %s ' % (config_variable, getattr(params, config_variable)))
    print('\n-------------------------------------------------------\n')
    gen_real_feig_const = lambda: gen_feig_consts(params.i)
    gen_real_percolation_const = lambda: gen_percolation_consts(params.i)
    gen_real_zeta_conts = lambda: gen_zeta_consts(params.i)

    consts_generators = {'e': gen_e_const,
                         'pi': gen_pi_const,
                         'feig': gen_real_feig_const,
                         'euler_masch': gen_euler_masch_const,
                         'percolation': gen_real_percolation_const,
                         'zeta': gen_real_zeta_conts}
    if params.const in consts_generators:
        target_generator = consts_generators[params.const]
    else:
        raise ValueError('Invalid const.')

    measure_runtime = MeasureRuntime()
    measure_runtime.start_measure()
    if params.const == 'feig':
        params.const = 'feig, %d' % i
    if params.const == 'percolation':
        params.const = 'percolation, %d' % i

    # Either loads the hashtable from previous runs or loads it
    if params.hashtable_file_operation in ['use', 'expand']:
        with open(params.hashtable_file, 'rb') as input_file:
            mitm = pickle.load(input_file)
        print('Loaded mitm object and hashtable from %s. Runtime: %s ' %
              (params.hashtable_file, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
        if params.hashtable_file_operation == 'expand':
            mitm.redefine_settings(target_generator=target_generator, target_name=params.const,
                                   postproc_funcs=EVALUATED_POSTPROC_FUNCS,
                                   postproc_funcs_filter=params.postproc_funcs_filter,
                                   ab_poly_class=params.ab_polys_type)
            print('Updated mitm object. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
    elif params.hashtable_file_operation == 'generate':
        mitm = enum_params.MITM(target_generator=target_generator, target_name=params.const,
                                postproc_funcs=EVALUATED_POSTPROC_FUNCS,
                                postproc_funcs_filter=params.postproc_funcs_filter,
                                hashtable_prec=params.hashtable_precision,
                                num_of_iterations=params.hashtable_num_of_iterations)
        print('Finished creating mitm object. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
    else:
        raise ValueError("hashtable_file_operation must be either 'use', 'expand' or 'generate'.")
    if params.hashtable_file_operation in ['expand', 'generate']:
        mitm.build_hashtable(range_a=params.a_coeffs_range, range_b=params.b_coeffs_range)
        print('Finished building hashtable. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
        # TODO: Any new parameters should be added and saved here?
        mitm.dec_hashtable['parameters'] = {'target_generator': target_generator, 'target_name': params.const,
                                            'POSTPROC_FUNCS': EVALUATED_POSTPROC_FUNCS}
        with open(params.hashtable_file, 'wb') as output_file:
            pickle.dump(mitm, output_file, protocol=pickle.HIGHEST_PROTOCOL)
        print('Stored hashtable as %s. Runtime: %s ' % (params.hashtable_file,
                                                        str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    if params.gen_hashtable_only:
        return
    mitm.find_clicks(params.lhs_type, params.lhs_params)
    print('Finished finding clicks. Number of clicks: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    mitm.reset_hashtable()
    mitm.filter_uniq_params()
    print('Finished filtering unique. Number of clicks: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    mitm.filter_only_exp_convergence(params.print_surprising_nonexp_contfracs)
    print('Finished fast filtering exponential convergence. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    mitm.filter_clicks_by_approach_type(whitelist=['exp', 'super_exp'])
    print('Finished full filtering exponential convergence. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    mitm.refine_clicks(accuracy=10, num_of_iterations=300, print_clicks=False)
    print('Finished refining clicks, 8 digits accuracy, 2000 iterations. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    print('---REFINING CAN BE CANCELLED NOW---')
    try:
        mitm.refine_clicks(accuracy=20, num_of_iterations=350, print_clicks=False)
        print('Finished refining clicks, 19 digits accuracy, 40000 iterations. Number of clicks left: %d. Runtime: %s ' %
              (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
        mitm.filter_uniq_params()
        mitm.filter_uniq_params()
        mitm.filter_uniq_params()
        print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
              (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
        mitm.filter_integer_roots_numerators()
        print('Finished filtering parameters with integer numerators roots. Number of unique parameters: %d. Runtime: %s ' %
              (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    except KeyboardInterrupt as e:
        print('Canceled refining clicks, 14 digits accuracy, 40000 iterations. Runtime: %s ' %
              (str(datetime.timedelta(seconds=measure_runtime.measure_time()))))

    if not os.path.isdir('results'):
        os.mkdir('results')
    os.mkdir(os.path.join('results', '%s') % time.strftime('%d%m%y_%H%M'))
    export_filename = os.path.join('results', '%s', 'results') % time.strftime('%d%m%y_%H%M')
    mitm.export_to_csv(export_filename + '.csv')
    print('Finished saving results. Filename: %s.csv. Runtime: %s ' %
          (export_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))

    # Next generate a PDF of the results (module is not finished)
    eqns = mitm.get_results_as_eqns(POSTPROC_FUNCS)
    generate_latex(export_filename, eqns)
    print('Generated PDF of results. Filename: %s.pdf. Runtime: %s ' %
          (export_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    # Save the configuration file to the results directory
    results_config = ConfigParser()
    results_config.add_section('Setup')
    for config_variable in LEGAL_CONFIG_PARAMS:
        results_config.set('Setup', config_variable, str(getattr(params, config_variable)))
    results_config_filename = os.path.split(export_filename)
    results_config_filename = os.path.join(results_config_filename[0],
                                           results_config_filename[1].replace('results', 'config.ini'))
    with open(results_config_filename, 'w') as results_config_file:
        results_config.write(results_config_file)
    print('Generated config file of the results. Filename: %s. Runtime: %s ' %
          (results_config_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))

# try:
if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()
# except:
#     pdb.pm()

