"""Defines, implements and provides required additional data (e.g. latex representation) for the post-process functions
to be applied on the result continued fractions."""

from decimal import Decimal as dec
from decmath import sin as dec_sin
from decmath import cos as dec_cos

def safe_inverse(x):
    """Returns 1/x if x != 0, otherwise returns infinity."""
    if x.is_zero():
        return dec('inf')
    else:
        return 1/x

def safe_sqrt(x):
    """Returns sqrt(x) if x >= 0, otherwise returns NaN"""
    if x < 0:
        return dec('NaN')
    else:
        return x.sqrt()


# All the supported postproc functions. New ones should be added at the end, to allow backward compatibility.
# IMPORTANT: after every new postproc func should appear it's inverse: safe_inverse(new_postproc_func).
POSTPROC_FUNCS = ['lambda x: x',
                  'safe_inverse',
                  'lambda x: x**2',
                  'lambda x: safe_inverse(x**2)',
                  'lambda x: x**3',
                  'lambda x: safe_inverse(x**3)',
                  'lambda x: x**4',
                  'lambda x: safe_inverse(x**4)',
                  'lambda x: x**5',
                  'lambda x: safe_inverse(x**5)',
                  'safe_sqrt',
                  'lambda x: safe_inverse(safe_sqrt(x))',
                  'lambda x: dec_sin(x)',
                  'lambda x: safe_inverse(dec_sin(x))',
                  'lambda x: dec_cos(x)',
                  'lambda x: safe_inverse(dec_cos(x))',
                  'lambda x: dec_tan(x)',
                  'lambda x: safe_inverse(dec_tan(x))',
                  'lambda x: dec_cot(x)',
                  'lambda x: safe_inverse(dec_cot(x))',
                  'lambda x: x.exp()',
                  'lambda x: safe_inverse(x.exp())',
                  'lambda x: x.ln()',
                  'lambda x: safe_inverse(x.ln())',]
                  # 'lambda x: dec_sin(x/18)',
                  # 'lambda x: dec_cos(x/18)',
                  # 'lambda x: 2*dec_sin(x/18)',
                  # 'lambda x: 2*dec_cos(x/18)',
                  # 'lambda x: dec_sin(2*safe_inverse(x)/9)',
                  # 'lambda x: dec_cos(2*safe_inverse(x)/9)',
                  # 'lambda x: 2*dec_sin(2*safe_inverse(x)/9)',
                  # 'lambda x: 2*dec_cos(2*safe_inverse(x)/9)']

# A list of pair indices of inverse postproc functions, for easier validation
INVERSE_POSTPROC_PAIRS = zip(POSTPROC_FUNCS[0::2], POSTPROC_FUNCS[1::2])
INVERSE_POSTPROC_PAIRS = [ (POSTPROC_FUNCS.index(f1), POSTPROC_FUNCS.index(f2)) for f1, f2 in INVERSE_POSTPROC_PAIRS ]
INVERSE_POSTPROC_PAIRS += [ funcs_pair[::-1] for funcs_pair in INVERSE_POSTPROC_PAIRS ]

# Latex representations of the funcs in POSTPROC_FUNCS. They are expected to be synchronized with the same keys.
# ** Please keep them with the same order too **
# TODO: change this to a joint object, e.g. a single list with Funcs objects, each with func and latex repr as members.
#  This solution may be problematic due to the way of usage. Requires examination.
POSTPROC_FUNCS_LATEX = {'lambda x: x': lambda lhs, rhs: (lhs, rhs),
                        'safe_inverse': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ {0} }}'.format(rhs)),
                        'lambda x: x**2': lambda lhs, rhs: (lhs, r'\left( {0} \right)^2'.format(rhs)),
                        'lambda x: safe_inverse(x**2)': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \left( {0} \right)^2 }}'.format(rhs)),
                        'lambda x: x**3': lambda lhs, rhs: (lhs, r'\left( {0} \right)^3'.format(rhs)),
                        'lambda x: safe_inverse(x**3)': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \left( {0} \right)^3 }}'.format(rhs)),
                        'lambda x: x**4': lambda lhs, rhs: (lhs, r'\left( {0} \right)^4'.format(rhs)),
                        'lambda x: safe_inverse(x**4)': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \left( {0} \right)^4 }}'.format(rhs)),
                        'lambda x: x**5': lambda lhs, rhs: (lhs, r'\left( {0} \right)^5'.format(rhs)),
                        'lambda x: safe_inverse(x**5)': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \left( {0} \right)^5 }}'.format(rhs)),
                        'safe_sqrt': lambda lhs, rhs: (r'\sqrt{ {0} }'.format(lhs), rhs),
                        'lambda x: safe_inverse(safe_sqrt(x))': lambda lhs, rhs: (lhs, r'\frac{{1}}{\sqrt{ {0} }}'.format(rhs)),
                        'lambda x: dec_sin(x)': lambda lhs, rhs: (lhs, r'\sin\left( {0} \right)'.format(rhs)),
                        'lambda x: safe_inverse(dec_sin(x))': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \sin\left( {0} \right) }}'.format(rhs)),
                        'lambda x: dec_cos(x)': lambda lhs, rhs: (lhs, r'\cos\left( {0} \right)'.format(rhs)),
                        'lambda x: safe_inverse(dec_cos(x))': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \cos\left( {0} \right) }}'.format(rhs)),
                        'lambda x: dec_tan(x)': lambda lhs, rhs: (lhs, r'\tan\left( {0} \right)'.format(rhs)),
                        'lambda x: safe_inverse(dec_tan(x))': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \tan\left( {0} \right) }}'.format(rhs)),
                        'lambda x: dec_cot(x)': lambda lhs, rhs: (lhs, r'\cot\left( {0} \right)'.format(rhs)),
                        'lambda x: safe_inverse(dec_cot(x))': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \cot\left( {0} \right) }}'.format(rhs)),
                        'lambda x: x.exp()': lambda lhs, rhs: (rhs, r'\exp\left( {0} \right)'.format(rhs)),
                        'lambda x: safe_inverse(x.exp())': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \exp\left( {0} \right) }}'.format(rhs)),
                        'lambda x: x.ln()': lambda lhs, rhs: (rhs, r'\ln\left( {0} \right)'.format(lhs)),
                        'lambda x: safe_inverse(x.ln())': lambda lhs, rhs: (lhs, r'\frac{{1}}{{ \ln\left( {0} \right) }}'.format(rhs)),
                        # 'lambda x: dec_sin(x/18)': lambda lhs, rhs: (lhs, r'\sin\left( \frac{{ {0} }}{{18}} \right)'.format(rhs)),
                        # 'lambda x: dec_cos(x/18)': lambda lhs, rhs: (lhs, r'\cos\left( \frac{{ {0} }}{{18}} \right)'.format(rhs)),
                        # 'lambda x: 2*dec_sin(x/18)': lambda lhs, rhs: (lhs, r'2\sin\left( \frac{{ {0} }}{{18}} \right)'.format(rhs)),
                        # 'lambda x: 2*dec_cos(x/18)': lambda lhs, rhs: (lhs, r'2\cos\left( \frac{{ {0} }}{{18}} \right)'.format(rhs)),
                        # 'lambda x: dec_sin(2*safe_inverse(x)/9)': lambda lhs, rhs: (lhs, r'\sin\left( \frac{{2}}{{9 \cdot {0}}} \right)'.format(rhs)),
                        # 'lambda x: dec_cos(2*safe_inverse(x)/9)': lambda lhs, rhs: (lhs, r'\cos\left( \frac{{2}}{{9 \cdot {0}}} \right)'.format(rhs)),
                        # 'lambda x: 2*dec_sin(2*safe_inverse(x)/9)': lambda lhs, rhs: (lhs, r'2\sin\left( \frac{{2}}{{9 \cdot {0}}} \right)'.format(rhs)),
                        # 'lambda x: 2*dec_cos(2*safe_inverse(x)/9)': lambda lhs, rhs: (lhs, r'2\cos\left( \frac{{2}}{{9 \cdot {0}}} \right)'.format(rhs)),
                        }

# Evaluates all the postproc funcs from strings to functions
EVALUATED_POSTPROC_FUNCS = [ eval(ppf) for ppf in POSTPROC_FUNCS ]

# Asserting all postproc funcs have a latex representation, and vice versa
assert(all([ ppf in POSTPROC_FUNCS_LATEX for ppf in POSTPROC_FUNCS ] +
           [ ppf_latex in POSTPROC_FUNCS for ppf_latex in POSTPROC_FUNCS_LATEX ]))