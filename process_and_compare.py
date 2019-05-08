"""This is a legacy file, and left here solely to be rewritten if needed. If one needs to convert csv to pdfs, load old
 results from a csv to ContFrac objects etc, updating this file may save some effort."""

import gen_consts
import cont_fracs
import csv
import re
import sys
from decimal import Decimal as dec
import decimal
import math
import matplotlib.pyplot as plt
import scipy.stats
import os
from gen_consts import gen_pi_const, gen_e_const, gen_feig_consts, gen_euler_masch_const, gen_percolation_consts
import enum_params
from postprocfuncs import EVALUATED_POSTPROC_FUNCS, POSTPROC_FUNCS
from lhs_evaluators import ULCDEvaluator
import json
from latex import generate_latex

dc = decimal.getcontext().prec=150

def old_enum_csv_to_pdf(csv_path, constant=None):
    csvfile_in = open(csv_path)
    csvreader = csv.reader(csvfile_in)
    contfrac_params = []
    postproc_funcs = []

    consts_generators = {'e': gen_e_const,
                         'pi': gen_pi_const,
                         'feig': gen_feig_consts,
                         'euler_masch': gen_euler_masch_const,
                         'percolation': gen_percolation_consts}

    if constant:
        if constant in consts_generators:
            constant_gen = consts_generators[constant]
        elif constant.startswith('feig'):
            i = int(constant.split(',')[1].strip())
            constant_gen = lambda: consts_generators['feig'](i)
        elif constant.startswith('percolation'):
            i = int(constant.split(',')[1].strip())
            constant_gen = lambda: consts_generators['percolation'](i)

    csv_iter = enumerate(csvreader)
    _, row = next(csv_iter)
    postproc_funcs = row[1].strip('[').strip(']').replace("'", "").split(', ')
    evaluated_postproc_funcs = [ EVALUATED_POSTPROC_FUNCS[POSTPROC_FUNCS.index(ppf)] for ppf in postproc_funcs]
    if not constant:
        if len(row) > 2 and row[3]:
            constant = row[3]
            if constant in consts_generators:
                constant_gen = consts_generators[constant]
            elif constant.startswith('feig'):
                _feig_index = int(constant.split(',')[1].strip())
                constant_gen = lambda: consts_generators['feig'](_feig_index)
            elif constant.startswith('percolation'):
                _percolation_index = int(constant.split(',')[1].strip())
                constant_gen = lambda: consts_generators['percolation'](_percolation_index)
        else:
            constant = 'pi'
            constant_gen = consts_generators[constant]
    _, row = next(csv_iter)

    mitm = enum_params.MITM(target_generator=constant_gen, target_name=constant,
                            postproc_funcs=evaluated_postproc_funcs)
    for i,row in enumerate(csvreader):
        a_params = json.loads(row[0].replace('(', '[').replace(',)', ']').replace(')', ']'))
        b_params = json.loads(row[1].replace('(', '[').replace(',)', ']').replace(')', ']'))
        if not isinstance(a_params[0], list):
            a_params = [a_params]
        if not isinstance(b_params[0], list):
            b_params = [b_params]
        u, l, c, d = [ int(c) for c in row[3:7] ]
        ulcd_obj = ULCDEvaluator((u, l, c, d), constant_gen())
        post_func_ind = int(row[2])
        ab = (a_params, b_params)
        mitm.filtered_params.append((ab, ulcd_obj, post_func_ind, None))

    pdf_path = csv_path.strip('.csv')
    eqns = mitm.get_results_as_eqns(postproc_funcs)
    generate_latex(pdf_path, eqns)
    print('Generated PDF of results. Filename: %s.pdf.' % (pdf_path))

# def load_enum_csv_as_cont_fracs_csv(csv_path, constant=None):
#     csvfile_in = open(csv_path)
#     csvreader = csv.reader(csvfile_in)
#     contfrac_params = []
#     postproc_funcs = []
#
#     consts_generators = {'e': gen_e_const,
#                          'pi': gen_pi_const,
#                          'feig': gen_feig_consts,
#                          'euler_masch': gen_euler_masch_const}
#
#     if constant in consts_generators:
#         constant = consts_generators[constant]()
#     elif constant.startswith('feig'):
#         i = int(constant.split(',')[1])
#         constant = consts_generators['feig'](i)
#
#     # WARNING: This block is a HUGE security issue!
#     for i,row in enumerate(csvreader):
#         if i == 0:
#             # postproc_funcs = eval(row[1])
#             # postproc_funcs = [ eval(p) for p in postproc_funcs ]
#             postproc_funcs = [safe_inverse, lambda x: x, lambda x: x**(type(x)(0.5))]
#             continue
#         if i == 1:
#             if not constant:
#                 constant = row[3]
#                 if constant in consts_generators:
#                     constant = consts_generators[constant]()
#                 elif constant.startswith('feig'):
#                     i = int(constant.split(',')[1])
#                     constant = consts_generators['feig'](i)
#
#         a_params = eval(row[0])
#         b_params = eval(row[1])
#         u, l, c, d = [ int(c) for c in row[3:7] ]
#         gen_ulcd_func = lambda u, l, c, d: lambda x: (u/x + x/l + c) / d
#         gen_func = lambda x, y: lambda z: x(y(z))
#         target_formula = gen_func(postproc_funcs[int(row[2])], gen_ulcd_func(u, l, c, d))
#         contfrac_params.append({'a_params': a_params, 'b_params': b_params, 'target_formula': target_formula,
#                            'target_x': constant})
#
#     return contfrac_params

def load_cont_fracs(csv_path):
    csvfile = open(csv_path)
    csvreader = csv.reader(csvfile)
    contfrac_params = []
    for i,row in enumerate(csvreader):
        csv_a_params = row[0]
        a_params = [ [ int(n) for n in p ] for p in re.findall('\\[(\\-?[0-9]+), ?(\\-?[0-9]+), ?(\\-?[0-9]+)\\]', csv_a_params) ]
        if not a_params:
            print('Error loading row number %d, incorrect a parameters' % (i+1), file=sys.stderr)
            continue
        csv_b_params = row[1]
        b_params = [ [ int(n) for n in p ] for p in re.findall('\\[(\\-?[0-9]+), ?(\\-?[0-9]+), ?(\\-?[0-9]+)\\]', csv_b_params) ]
        if not b_params:
            print('Error loading row number %d, incorrect b parameters' % (i+1), file=sys.stderr)
            continue
        csv_target_formula = row[2].strip()
        if not re.match('[0-9+\\-*/x ()]+', csv_target_formula):
            print('Error loading row number %d, incorrect target formula' % (i+1), file=sys.stderr)
            continue
        try:
            target_formula = eval('lambda x: ' + csv_target_formula)
        except:
            print('Error loading row number %d, incorrect target formula, evaluation error' % (i+1), file=sys.stderr)
            continue
        csv_target_x = row[3].strip().lower()
        if csv_target_x == 'pi':
            target_x = gen_consts.gen_pi_const()
        elif csv_target_x == 'e':
            target_x = gen_consts.gen_e_const()
        elif re.match('^\\-?[0-9]+(\\.[0-9]+)?$', csv_target_x):
            target_x = dec(csv_target_x)
        else:
            print('Error loading row number %d, incorrect target x' % (i+1), file=sys.stderr)
            continue
        contfrac_params.append({'a_params': a_params, 'b_params': b_params, 'target_formula': target_formula,
                           'target_x': target_x})

    return contfrac_params


def build_contfrac_from_params(params):
    return cont_fracs.ContFrac(a_coeffs=params['a_params'], b_coeffs=params['b_params'],
                               target_val=params['target_formula'](params['target_x']))


def build_contfrac_errors(contfrac, iters=2000, show_progress=False):
    """Receives ContFrac object"""
    errors = []
    contfrac.gen_iterations(0)
    for i in range(iters):
        contfrac.add_iterations(1)
        errors.append((i+1, contfrac.compare_result()))
        if show_progress and i % 500 == 0:
            print('\r%d' % i, end='')
    contfrac.gen_iterations(0)
    return errors


def plot_errors(errors, title=None, show_plot=True, lin_fit=None, save_fig_path=None, matplot_interactive=True,
                figsize=(15, 7)):
    """errors = [(i,err(i))]"""
    effective_zero = dec('1E-140')
    if matplot_interactive:
        plt.ion()
    else:
        plt.ioff()

    errors = [ (i, math.log(i), math.log(err)) for i, err in errors if abs(err)> effective_zero ]
    errors_y = [e[2] for e in errors ]
    # fig = plt.figure(figsize=figsize)
    if not title:
        title = 'Errors plot'
    plt.suptitle(title)

    plt.subplot(121)
    plt.title('Exponential')
    errors_x = [ e[0] for e in errors ]
    if lin_fit:
        intercept = lin_fit['exponential']['intercept']
        slope = lin_fit['exponential']['slope']
        r_square = lin_fit['exponential']['r_value']**2
        linfit_y = [ intercept + slope * x for x in errors_x ]
        plt.plot(errors_x, errors_y, 'r.', label='error')
        plt.plot(errors_x, linfit_y, 'g-', label=r'$lin.\ fit.\ R^2=%f\ slope=%f$' % (r_square, slope))
        plt.legend(prop={'size': 6})
    else:
        plt.plot(errors_x, errors_y, 'r.')
    plt.ylabel(r'$\log(error)$')
    plt.xlabel(r'$iter-num$')

    plt.subplot(122)
    plt.title('Polynomial')
    errors_x = [ e[1] for e in errors ]
    if lin_fit:
        intercept = lin_fit['polynomial']['intercept']
        slope = lin_fit['polynomial']['slope']
        r_square = lin_fit['polynomial']['r_value']**2
        linfit_y = [ intercept + slope * x for x in errors_x ]
        plt.plot(errors_x, errors_y, 'r.', label='error')
        plt.plot(errors_x, linfit_y, 'g-', label=r'$lin.\ fit.\ R^2=%f\ slope=%f$' % (r_square, slope))
        plt.legend(prop={'size': 6})
    else:
        plt.plot(errors_x, errors_y, 'r.')
    plt.ylabel(r'$\log(error)$')
    plt.xlabel(r'$\log(iter-num)$')

    if show_plot:
        plt.show()
    if save_fig_path:
        plt.savefig(save_fig_path, dpi=200)

    if not show_plot:
        plt.clf()
        plt.close()


def find_linear_fitting(errors, initial_cutoff=500):
    effective_zero = dec('1E-50')
    linear_fitting = {}
    while errors[initial_cutoff][1] < 10**-15 and initial_cutoff > 0:
        initial_cutoff >>= 1
    errors = errors[initial_cutoff:]

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress([ x for x,err in errors
                                                                           if abs(err)> effective_zero ],
                                                                         [ math.log(err) for x, err in errors
                                                                           if abs(err)> effective_zero ])
    linear_fitting['exponential'] = {'slope': slope, 'intercept': intercept, 'r_value': r_value,
                                     'p_value': p_value, 'std_err': std_err}

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress([ math.log(x) for x, err in errors
                                                                           if abs(err)> effective_zero ],
                                                                         [ math.log(err) for x, err in errors
                                                                           if abs(err)> effective_zero ])
    linear_fitting['polynomial'] = {'slope': slope, 'intercept': intercept, 'r_value': r_value,
                                     'p_value': p_value, 'std_err': std_err}
    return linear_fitting


def analyze_contfracs_csv(csv_path, filetype='enum_res'):
    """filetype in ['enum_res', 'cont_frac']"""
    dirname = os.path.dirname(csv_path)
    filename = os.path.basename(csv_path)
    filename_no_ext = os.path.splitext(filename)[0]
    if filetype == 'enum_res':
        cfs_params = load_enum_csv_as_cont_fracs_csv(csv_path)
    elif filetype == 'cont_frac':
        cfs_params = load_cont_fracs(csv_path)
    else:
        raise ValueError('Invalid value for filetype: %s' % str(filetype))
    conv_types = []
    conv_params = []
    for i,cf_params in enumerate(cfs_params):
        cf = build_contfrac_from_params(cf_params)
        errors = build_contfrac_errors(cf, iters=2000, show_progress=False)
        lf = find_linear_fitting(errors)
        output_fig_path = os.path.join(dirname, '%s-%d.png' % (filename_no_ext, i+1))
        fig_title = '%s - %dth row' % (filename_no_ext, i+1)
        plot_errors(errors, title=fig_title,
                    show_plot=False, save_fig_path=output_fig_path, lin_fit=lf)
        conv_types.append(cf.is_convergence_fast())
        cf.estimate_approach_type_and_params()
        conv_params.append(cf.get_approach_type_and_params())
    return conv_types, conv_params

# the exponential_threshold=1.25 fits for an exponential growth parameter of ~ 1.12 = sqrt(1.25)
# see the convergence analysis pdf for more details.
def is_poly_or_exp_approach(contfrac, iters=5000, initial_cutoff=1500, iters_step=500, exponential_threshold=1.1,
                            find_poly_parameter=False):
    """Returns 'exp', 'poly2sympoly', 'undefined', 'fast' and 'mixed', as a tuple of (string,num): (approach_type, approach_parameter)
or ('poly2sympoly', (approach_parameter, R**2))."""
    if iters_step < 6:
        ValueError('iters_step should be at least 4')

    approach_type = None
    approach_parameter = 0

    delta_pair = []
    delta_odd = []
    contfrac.gen_iterations(initial_cutoff)
    res_0 = contfrac.get_result()
    contfrac.add_iterations(1)
    res_1 = contfrac.get_result()
    contfrac.add_iterations(1)
    res_2 = contfrac.get_result()
    contfrac.add_iterations(1)
    res_3 = contfrac.get_result()
    contfrac.add_iterations(1)
    res_4 = contfrac.get_result()
    contfrac.add_iterations(1)
    res_5 = contfrac.get_result()
    delta_pair.append((initial_cutoff, abs(res_2 - res_0)))
    delta_pair.append((initial_cutoff + 2, abs(res_4 - res_2)))
    delta_odd.append((initial_cutoff + 1, abs(res_3 - res_1)))
    delta_odd.append((initial_cutoff + 3, abs(res_5 - res_3)))

    for i in range(initial_cutoff+iters_step, iters+1, iters_step):
        # -3 for the iterations of res_1, res_2, res_3 that were already executed
        contfrac.add_iterations(iters_step - 5)
        res_0 = contfrac.get_result()
        contfrac.add_iterations(1)
        res_1 = contfrac.get_result()
        contfrac.add_iterations(1)
        res_2 = contfrac.get_result()
        contfrac.add_iterations(1)
        res_3 = contfrac.get_result()
        contfrac.add_iterations(1)
        res_4 = contfrac.get_result()
        contfrac.add_iterations(1)
        res_5 = contfrac.get_result()
        delta_pair.append((i, abs(res_2 - res_0)))
        delta_pair.append((i + 2, abs(res_4 - res_2)))
        delta_odd.append((i + 1, abs(res_3 - res_1)))
        delta_odd.append((i + 3, abs(res_5 - res_3)))
        # if show_progress and i % 500 == 0:
        #     print('\r%d' % i, end='')

    pair_diminish = False
    odd_diminish = False
    if len(delta_pair) > 3 and all([ p[1].is_zero() for p in delta_pair[-3:] ]):
        pair_diminish = True
    if len(delta_odd) > 3 and all([ p[1].is_zero() for p in delta_odd[-3:] ]):
        odd_diminish = True
    # if one diminishes and the other isn't, return 'undefined'
    if pair_diminish ^ odd_diminish:
        approach_type = 'undefined'
    elif pair_diminish and odd_diminish:
        approach_type = 'fast'

    # if approach_type:
    #     return (approach_type, approach_parameter)

    pair_ratio = [ (delta_pair[i][0], delta_pair[i][1] / delta_pair[i+1][1])
                   for i in range(0, len(delta_pair), 2) if delta_pair[i][1] != 0 and delta_pair[i+1][1] != 0 ]
    odd_ratio = [ (delta_odd[i][0], delta_odd[i][1] / delta_odd[i+1][1])
                  for i in range(0, len(delta_odd), 2) if delta_odd[i][1] != 0 and delta_odd[i+1][1] != 0 ]

    if len(pair_ratio) == 0:
        return (approach_type, approach_parameter)

    mean_pair_ratio = sum([ p for i, p in pair_ratio] ) / len(pair_ratio)
    mean_pair_ratio_avg_square_error = sum([ (r-mean_pair_ratio)**2 for i, r in pair_ratio ]) / len(pair_ratio)
    mean_odd_ratio = sum([ p for i, p in odd_ratio ]) / len(odd_ratio)
    mean_odd_ratio_avg_square_error = sum([ (r-mean_odd_ratio)**2 for i, r in odd_ratio ]) / len(odd_ratio)
    relative_pair_sq_err = mean_pair_ratio_avg_square_error / mean_pair_ratio
    relative_odd_sq_err = mean_odd_ratio_avg_square_error / mean_odd_ratio
    if relative_odd_sq_err > 0.5 or relative_pair_sq_err > 0.5:
        approach_type = 'undefined'
    else:
        is_pair_exp = mean_pair_ratio > exponential_threshold
        is_odd_exp = mean_odd_ratio > exponential_threshold
        print('mean_pair_ratio', mean_pair_ratio)
        print('mean_odd_ratio', mean_odd_ratio)
        # in case one is exponential and the other isn't return 'mixed'
        if is_pair_exp ^ is_odd_exp:
            approach_type = 'mixed'
        elif is_pair_exp and is_odd_exp:
            approach_type = 'exp'
            approach_parameter = min(mean_pair_ratio**type(mean_pair_ratio)(0.5),
                                     mean_odd_ratio**type(mean_odd_ratio)(0.5))
        else:
            approach_type = 'poly2sympoly'

    if approach_type != 'poly2sympoly' or not find_poly_parameter:
        return (approach_type, approach_parameter)

#     We're requested to find the poly2sympoly parameter
    log_x_pair = [ math.log(i) for i, d in delta_pair ]
    log_y_pair = [ math.log(d) for i, d in delta_pair ]
    slope_pair, intercept_pair, r_value_pair, p_value_pair, std_err_pair = scipy.stats.linregress(log_x_pair,
                                                                                                  log_y_pair)
    plt.ion()
    plt.plot(log_x_pair, log_y_pair)
    plt.show()

    log_x_odd = [ math.log(i) for i, d in delta_odd ]
    log_y_odd = [ math.log(d) for i, d in delta_odd ]
    slope_odd, intercept_odd, r_value_odd, p_value_odd, std_err_odd = scipy.stats.linregress(log_x_odd, log_y_odd)

    approach_parameter = (min(abs(slope_pair), abs(slope_odd))-1, min(r_value_pair**2, r_value_odd**2))
    return (approach_type, approach_parameter)
