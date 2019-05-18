import gen_consts
from cont_fracs import ContFrac
from mpmath import mpf as dec, log10

import matplotlib.pyplot as plt
import scipy.stats

pi = gen_consts.gen_pi_const(100)
e = gen_consts.gen_e_const(100)

# unknown. super exponential
target1 = 1/(e/2 - 1)
a1 = ((2, 2),)
b1 = ((0, 4),)
cf1 = ContFrac(a1, b1, target_val=target1, logging=True)

# unknown. exponential
target2 = 4/pi
a2 = ((1, 3),)
b2 = ((0, 3, -2),)
cf2 = ContFrac(a2, b2, target_val=target2, logging=True)

# known.exponential
target3 = 4 / pi
a3 = ((1, 2),)
b3 = ((0, 0, 1),)
cf3 = ContFrac(a3, b3, target_val=target3, logging=True)

# known, polynomial
target4 = pi+3
a4 = ((6,),)
b4 = ((1, -4, 4),)
cf4 = ContFrac(a4, b4, target_val=target4, logging=True)

def plot_errors(errors, title=None, show_plot=True, lin_fit=None, save_fig_path=None, matplot_interactive=True,
                figsize=(15, 7), exp_or_poly='poly', iters_boundaries=[0, -1]):
    """errors = [(i,err(i))]"""
    effective_zero = dec('1E-140')
    if matplot_interactive:
        plt.ion()
    else:
        plt.ioff()

    errors = [ (i, log10(i), log10(dec(err))) for i, err in errors if abs(err)> effective_zero ]
    errors_y = [e[2] for e in errors[:iters_boundaries[1]] ]
    # fig = plt.figure(figsize=figsize)
    if not title:
        title = 'Errors plot'
    plt.suptitle(title)

    if exp_or_poly == 'exp':
        plt.title('Exponential')
        errors_x = [ e[0] for e in errors[:iters_boundaries[1]] ]
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
        plt.ylabel(r'$\log_{10}(error)$')
        plt.xlabel(r'$iter-num$')

    elif exp_or_poly == 'poly':
        plt.title('Polynomial')
        errors_x = [ e[1] for e in errors[:iters_boundaries[1]] ]
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
        plt.ylabel(r'$\log_{10}(error)$')
        plt.xlabel(r'$\log_{10}(iter-num)$')
    else:
        raise ValueError('exp_or_poly = %s' % exp_or_poly)

    if show_plot:
        plt.show()
    if save_fig_path:
        plt.savefig(save_fig_path, dpi=600, format='svg')

    if not show_plot:
        plt.clf()
        plt.close()


def find_linear_fitting(errors, iters_boundaries=[0, -1]):
    effective_zero = dec('1E-100')
    linear_fitting = {}
    # while errors[initial_cutoff][1] < 10**-15 and initial_cutoff > 0:
    #     initial_cutoff >>= 1
    errors = errors[iters_boundaries[0] : iters_boundaries[1]]

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress([ float(x) for x,err in errors
                                                                           if abs(err)> effective_zero ],
                                                                         [ float(log10(err)) for x, err in errors
                                                                           if abs(dec(err))> effective_zero ])
    linear_fitting['exponential'] = {'slope': slope, 'intercept': intercept, 'r_value': r_value,
                                     'p_value': p_value, 'std_err': std_err}

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress([ float(log10(dec(x))) for x, err in errors
                                                                           if abs(err)> effective_zero ],
                                                                         [ float(log10(dec(err))) for x, err in errors
                                                                           if abs(err)> effective_zero ])
    linear_fitting['polynomial'] = {'slope': slope, 'intercept': intercept, 'r_value': r_value,
                                    'p_value': p_value, 'std_err': std_err}
    return linear_fitting


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


def analyze_contfrac(cf, output_fig, exp_or_poly, iters_boundaries):
    iters = 2000
    if iters_boundaries[1] > 2000:
        iters = iters_boundaries[1]
    errors = build_contfrac_errors(cf, iters=iters, show_progress=False)
    lf = find_linear_fitting(errors, iters_boundaries=iters_boundaries)
    plot_errors(errors, title='',
                show_plot=True, save_fig_path=output_fig, lin_fit=lf, exp_or_poly=exp_or_poly,
                iters_boundaries=iters_boundaries)
    return lf

lf = analyze_contfrac(cf1, output_fig='errors_cf1.svg', exp_or_poly='exp', iters_boundaries=[30, 65])
lf = analyze_contfrac(cf2, output_fig='errors_cf2.svg', exp_or_poly='exp', iters_boundaries=[60, 300])
lf = analyze_contfrac(cf3, output_fig='errors_cf3.svg', exp_or_poly='exp', iters_boundaries=[30, 130])
lf = analyze_contfrac(cf4, output_fig='errors_cf4.svg', exp_or_poly='poly', iters_boundaries=[30, 130])
