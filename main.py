#!/usr/bin/env python3

# TODO: check for integer roots for p \in Z_3,4[x]
# TODO: fix the convergence params (usecase: the e contfracs: 'results_5303_050418.csv'
# this is a temporary execution file
import enum_params
from decimal import Decimal as dec
import time
import datetime
from gen_real_consts import gen_real_pi, gen_real_e, gen_real_feig, gen_real_euler_masch, gen_real_percolation
from decmath import sin as dec_sin
import dill as pickle  # Makes sure pickle also supports Lambda functions
from latex import generate_latex
import os
from configfile import ConfigParser
import configfile

LEGAL_CONFIG_PARAMS = configfile.CONFIG_PARAMS_TYPES

class MeasureRuntime():
    def __init__(self):
        self.times = []

    def start_measure(self):
        self.times.append(time.time())

    def measure_time(self):
        self.times.append(time.time())
        return self.times[-1] - self.times[-2]


class Parameters:
    pass

# example - this is the standard continuous fraction with sixes. it works.
# mitm = enum_params.MITM(postproc_func=lambda x:x)
# mitm.build_hashtable(range_a = [[[6,7], [0,1], [0,1]]], range_b=[[[1,2],[-4,-3],[4,5]]])
# mitm.find_clicks(params.u_range=[0], params.l_range=[1], params.c_range=range(-4,4), params.d_range=[1])
# mitm.refine_clicks()
# print(mitm.filtered_params)

# example - this is for testing a specific set of parameters
# bep = enum_params.BasicEnumPolyParams()
# bep.polys_generator(range_a=[[6,7], [-1,0], [1,2]], range_b=[[5,6],[4,5],[-3,-2]])

# example - this is for interlace
# main(a_coeffs_range=[[[2, 3]], [[5, 6]]], b_coeffs_range=[[[1, 2]], [[-9, -7]], [[1,2]]], a_interlace=2, b_interlace=3)

# default is a,b in Z_2[x]
def safe_inverse(x):
    if x.is_zero():
        return dec('inf')
    else:
        return 1/x

def safe_sqrt(x):
    if x < 0:
        return dec('NaN')
    else:
        return x**dec('0.5')

# poly_coeffs_range=3, ulcd_params=3, const='e', a_poly_size=3, b_poly_size=3, a_interlace=1, b_interlace=1,
         # print_surprising_nonexp_contfracs=False, a_coeffs_range=None, b_coeffs_range=None, params.u_range=None, params.l_range=None,
         # params.c_range=None, params.d_range=None, i=0
def main(configfile='config.ini'):
    """supported consts: pi, e, feig(0-3), euler_masch, percolation (0-1). for feig, i=0,1,2,3 is required.
    for percolation, i=0,1 is required"""

    # Load configuration file
    params = Parameters()
    config_parser = ConfigParser(configfile=configfile)
    config = config_parser.get_config()
    # Load variables from config
    for config_variable in LEGAL_CONFIG_PARAMS:
        # print('Setting up %s' % config_variable)
        setattr(params, config_variable, config[config_variable])
        print('%s is set up to %s ' % (config_variable, getattr(params, config_variable)))
    # poly_coeffs_range = config['poly_coeffs_range']
    # ulcd_params = config['ulcd_params']
    # const = config['const']
    # a_poly_size = config['a_poly_size']
    # b_poly_size = config['b_poly_size']
    # a_interlace = config['a_interlace']
    # b_interlace = config['b_interlace']
    # print_surprising_nonexp_contfracs = config['print_surprising_nonexp_contfracs']
    # a_coeffs_range = config['a_coeffs_range']
    # b_coeffs_range = config['b_coeffs_range']
    # params.u_range = config['params.u_range']
    # params.l_range = config['params.l_range']
    # params.c_range = config['params.c_range']
    # params.d_range = config['params.d_range']
    # i = config['i']

    # if not a_coeffs_range:
    #     a_coeffs_range = poly_coeffs_range
    # if not b_coeffs_range:
    #     b_coeffs_range = poly_coeffs_range
    if not params.u_range:
        params.u_range = params.ulcd_range
    if not params.l_range:
        params.l_range = params.ulcd_range
    if not params.c_range:
        params.c_range = params.ulcd_range
    if not params.d_range:
        params.d_range = params.ulcd_range

    gen_real_feig_const = lambda: gen_real_feig(i)
    gen_real_percolation_const = lambda: gen_real_percolation(i)

    consts_generators = {'e': gen_real_e,
                         'pi': gen_real_pi,
                         'feig': gen_real_feig_const,
                         'euler_masch': gen_real_euler_masch,
                         'percolation': gen_real_percolation_const}
    if params.const in consts_generators:
        target_generator = consts_generators[params.const]
    else:
        raise ValueError('Invalid const.')
    
    # our "defult"
    postproc_funcs = ['safe_inverse', 'lambda x: x', 'lambda x: x**2', 'lambda x: safe_inverse(x**2)', 'safe_sqrt']

    # trying sqrt on e
    # postproc_funcs = ['safe_inverse', 'lambda x: x', 'lambda x: x.sqrt()', 'lambda x: safe_inverse(x.sqrt())']
 
    # added to capture the percolation constants
    # postproc_funcs = ['lambda x: dec_sin(x/18)', 'lambda x: 2*dec_sin(x/18)', 'lambda x: dec_sin(2*safe_inverse(x)/9)', 'lambda x: 2*dec_sin(2*safe_inverse(x)/9)']

    # added just for fun:
    # postproc_funcs = ['lambda x: x**3', 'lambda x: safe_inverse(x**3)', 'lambda x: x**4', 'lambda x: safe_inverse(x**4)', 'lambda x: x**5', 'lambda x: safe_inverse(x**5)']

    evaluated_postproc_funcs = [ eval(ppf) for ppf in postproc_funcs ]
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
            mitm.redefine_settings(target_generator=target_generator, target_name=params.const, a_poly_size=params.a_poly_size,
                                   b_poly_size=params.b_poly_size, num_of_a_polys=params.a_interlace, num_of_b_polys=params.b_interlace,
                                   postproc_funcs=evaluated_postproc_funcs)
            print('Updated mitm object. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
    elif params.hashtable_file_operation == 'generate':
        mitm = enum_params.MITM(target_generator=target_generator, target_name=params.const, a_poly_size=params.a_poly_size,
                                b_poly_size=params.b_poly_size, num_of_a_polys=params.a_interlace, num_of_b_polys=params.b_interlace,
                                postproc_funcs=evaluated_postproc_funcs)
        print('Finished creating mitm object. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
    if params.hashtable_file_operation in ['expand', 'generate']:
        # a,b polynoms coefficients will be enumerated in [-2,2]
        # one can either set enum_range to set a uniform boundary to all the coefficients,
        # or set a different range to the a's coefficients and b's coefficients.
        # the given value should be either int (then the range will be [-a,a], enumeration includes both edges), or a 2-elements tuple/list
        # of the form [a,b] where a<b. enumeration includes only lower edge (b isn't included)
        mitm.build_hashtable(enum_range=params.poly_coeffs_range, range_a=params.a_coeffs_range, range_b=params.b_coeffs_range)
        print('Finished building hashtable. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
        mitm.dec_hashtable['parameters'] = {'target_generator': target_generator, 'target_name': params.const,
                                            'a_poly_size': params.a_poly_size, 'b_poly_size': params.b_poly_size,
                                            'num_of_a_polys': params.a_interlace, 'num_of_b_polys': params.b_interlace,
                                            'postproc_funcs': evaluated_postproc_funcs}
        with open(params.hashtable_file, 'wb') as output_file:
            pickle.dump(mitm, output_file, protocol=pickle.HIGHEST_PROTOCOL)
        print('Stored hashtable as %s. Runtime: %s ' % (params.hashtable_file,
                                                        str(datetime.timedelta(seconds=measure_runtime.measure_time()))))

    # for finding clicks, we enumerate u,l,c,d: (u/pi+pi/l+c)*1/d
    # here a range should e either an int (then the enumeration is over [-i,i]), or an iterable of any type
    # (e.g. list, range object etc.)
    mitm.find_clicks(u_range=params.u_range, l_range=params.l_range, c_range=params.c_range, d_range=params.d_range)
    mitm.delete_hashtable()
    mitm.filter_uniq_params()
    print('Finished finding clicks. Number of clicks: %d. Runtime: %s ' %
          (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    mitm.filter_only_exp_convergence(params.print_surprising_nonexp_contfracs)
    print('Finished fast filtering exponential convergence. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    mitm.filter_clicks_by_approach_type()
    print('Finished full filtering exponential convergence. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    mitm.refine_clicks(accuracy=10, num_of_iterations=4000, print_clicks=False)
    print('Finished refining clicks, 8 digits accuracy, 2000 iterations. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    # mitm.filter_uniq_params()
    # print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
    #       (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    # mitm.refine_clicks(accuracy=12, num_of_iterations=10000)
    # print('Finished refining clicks, 12 digits accuracy, 10000 iterations. Number of clicks left: %d. Runtime: %s ' %
    #       (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    # mitm.filter_uniq_params()
    # print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
    #       (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    print('---REFINING CAN BE CANCELLED NOW---')
    try:
        mitm.refine_clicks(accuracy=20, num_of_iterations=2000, print_clicks=False)
        print('Finished refining clicks, 19 digits accuracy, 40000 iterations. Number of clicks left: %d. Runtime: %s ' %
              (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
        mitm.filter_uniq_params()
        mitm.filter_uniq_params()
        mitm.filter_uniq_params()
        print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
              (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    except KeyboardInterrupt as e:
        print('Canceled refining clicks, 14 digits accuracy, 40000 iterations. Runtime: %s ' %
              (str(datetime.timedelta(seconds=measure_runtime.measure_time()))))

    if not os.path.isdir('results'):
        os.mkdir('results')
    os.mkdir(os.path.join('results', '%s') % time.strftime('%d%m%y_%H%M'))
    export_filename = os.path.join('results', '%s', 'results') % time.strftime('%d%m%y_%H%M')
    mitm.export_to_csv(export_filename + '.csv', postproc_funcs)
    print('Finished saving results. Filename: %s.csv. Runtime: %s ' %
          (export_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))

    # Next generate a PDF of the results (module is not finished)
    eqns = mitm.get_results_as_eqns(postproc_funcs)
    generate_latex(export_filename + '.pdf', eqns)
    print('Generated PDF of results. Filename: %s.pdf. Runtime: %s ' %
          (export_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    # Save the configuration file to the results directory
    results_config = ConfigParser()
    results_config.add_section('Setup')
    # results_config.add_section('Data')
    for config_variable in LEGAL_CONFIG_PARAMS:
        results_config.set('Setup', config_variable, str(getattr(params, config_variable)))
    # for config_variable in LEGAL_CONFIG_PARAMS[-2:]:
    #     results_config.set('Data', config_variable, getattr(params, config_variable))
    results_config_filename = os.path.split(export_filename)
    results_config_filename = os.path.join(results_config_filename[0],
                                           results_config_filename[1].replace('results', 'config.ini'))
    with open(results_config_filename, 'w') as results_config_file:
        results_config.write(results_config_file)
    print('Generated config file of the results. Filename: %s. Runtime: %s ' %
          (results_config_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time()))))


if __name__ == '__main__':
    main()
