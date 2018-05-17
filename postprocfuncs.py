from decimal import Decimal as dec
from decmath import sin as dec_sin

def safe_inverse(x):
    if x.is_zero():
        return dec('inf')
    else:
        return 1/x

def safe_sqrt(x):
    if x < 0:
        return dec('NaN')
    else:
        return x.sqrt()

        # poly_coeffs_range=3, ulcd_params=3, const='e', a_poly_size=3, b_poly_size=3, a_interlace=1, b_interlace=1,
        # print_surprising_nonexp_contfracs=False, a_coeffs_range=None, b_coeffs_range=None, params.u_range=None, params.l_range=None,
        # params.c_range=None, params.d_range=None, i=0


# our "defult"
POSTPROC_FUNCS = ['safe_inverse', 'lambda x: x', 'lambda x: x**2', 'lambda x: safe_inverse(x**2)', 'safe_sqrt']

# trying sqrt on e
# POSTPROC_FUNCS = ['safe_inverse', 'lambda x: x', 'lambda x: x.sqrt()', 'lambda x: safe_inverse(x.sqrt())']

# added to capture the percolation constants
# POSTPROC_FUNCS = ['lambda x: dec_sin(x/18)', 'lambda x: 2*dec_sin(x/18)', 'lambda x: dec_sin(2*safe_inverse(x)/9)', 'lambda x: 2*dec_sin(2*safe_inverse(x)/9)']

# added just for fun:
# POSTPROC_FUNCS = ['lambda x: x**3', 'lambda x: safe_inverse(x**3)', 'lambda x: x**4', 'lambda x: safe_inverse(x**4)', 'lambda x: x**5', 'lambda x: safe_inverse(x**5)']

POSTPROC_FUNCS_LATEX = {'safe_inverse': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'lambda x: x': lambda lhs, rhs: (lhs, rhs),
                        'lambda x: x**2': lambda lhs, rhs: (r'{{ {0} }}^2'.format(x), rhs),
                        'lambda x: x**3': lambda lhs, rhs: (r'{{ {0} }}^3'.format(x), rhs),
                        'lambda x: x**4': lambda lhs, rhs: (r'{{ {0} }}^4'.format(x), rhs),
                        'lambda x: x**5': lambda lhs, rhs: (r'{{ {0} }}^5'.format(x), rhs),
                        'lambda x: safe_inverse(x**2)': lambda lhs, rhs: (r'{{ {0} }}^2'.format(lhs), r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'lambda x: safe_inverse(x**3)': lambda lhs, rhs: (r'{{ {0} }}^3'.format(lhs), r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'lambda x: safe_inverse(x**4)': lambda lhs, rhs: (r'{{ {0} }}^4'.format(lhs), r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'lambda x: safe_inverse(x**5)': lambda lhs, rhs: (r'{{ {0} }}^5'.format(lhs), r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'safe_sqrt': lambda lhs, rhs: (r'\sqrt{ {0} }'.format(lhs), rhs),
                        'lambda x: safe_inverse(x.sqrt())': lambda lhs, rhs: (r'\sqrt{{ {0} }}'.format(lhs), r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'lambda x: dec_sin(x)': lambda lhs, rhs: (r'\sin\left( {0} \right)'.format(x), rhs),
                        'lambda x: dec_cos(x)': lambda lhs, rhs: (r'\cos\left( {0} \right)'.format(x), rhs),
                        'lambda x: dec_tan(x)': lambda lhs, rhs: (r'\tan\left( {0} \right)'.format(x), rhs),
                        'lambda x: dec_cot(x)': lambda lhs, rhs: (r'\cot\left( {0} \right)'.format(x), rhs),
                        'lambda x: x.exp()': lambda lhs, rhs: (r'\exp\left( {0} \right)'.format(x), rhs),
                        'lambda x: x.ln()': lambda lhs, rhs: (r'\ln\left( {0} \right)'.format(x), rhs),
                        }
EVALUATED_POSTPROC_FUNCS = [ eval(ppf) for ppf in POSTPROC_FUNCS ]