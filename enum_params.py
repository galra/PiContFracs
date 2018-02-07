import itertools
import cont_fracs
from basic_algo import set_precision
from decimal import Decimal as dec


class BasicEnumPolyParams:
    def __init__(self, a_poly_size=3, b_poly_size=3, num_of_iterations=15, threshold=None, prec=1000):
        self.a_poly_size = a_poly_size
        self.b_poly_size = b_poly_size
        self.num_of_iterations = num_of_iterations
        self.good_params = []
        if threshold is None:
            self.threshold = dec(10)**dec(-4)
        else:
            self.threshold = threshold
        self._prec = prec
        set_precision(prec)

    def reset_precision(self):
        set_precision(self._prec)

    def reinitialize_good_params(self):
        self.good_params = []

    def pis_generator(self, enum_range=None, range_a=None, range_b=None, prec=None):
        if not (enum_range or range_a or range_b):
            raise ValueError("No range was given")
        if isinstance(enum_range, list):
            range_a = [enum_range for i in range(self.a_poly_size) ]
            range_b = [enum_range for i in range(self.b_poly_size) ]
        elif enum_range:
            range_a = [ [-enum_range, enum_range] for i in range(self.a_poly_size) ]
            range_b = [ [-enum_range, enum_range] for i in range(self.b_poly_size) ]

        pi_cont_frac = cont_fracs.PiContFrac([ 0 for i in range_a ], [ 0 for i in range_b ])
        a_params_iterator = itertools.product(*[ range(*r) for r in range_a ])
        for pa in a_params_iterator:
            b_params_iterator = itertools.product(*[ range(*r) for r in range_b ])
            for pb in b_params_iterator:
                pi_cont_frac.reinitialize(pa, pb)
                pi_cont_frac.gen_iterations(self.num_of_iterations)
                yield (pi_cont_frac, pa, pb)

    def enum_params(self, enum_range=None, range_a=None, range_b=None, show_progress=True):
        params_iter = self.pis_generator(enum_range, range_a, range_b)
        iter_num = 1
        for pi_cont_frac, pa, pb in params_iter:
            if pi_cont_frac.is_pi_valid() and pi_cont_frac.compare_result() < self.threshold:
                    self.good_params.append({'a': pa, 'b': pb})
            if iter_num % 100 == 0 and show_progress:
                print('\r%d' % (iter_num), end='')
            iter_num += 1
        if show_progress:
            print('')
