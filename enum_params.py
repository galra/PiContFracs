from basic_algo import set_precision
import cont_fracs
from decimal_hashtable import DecimalHashTable
from decimal import Decimal as dec
import itertools
from gen_real_pi import gen_real_pi
import time

class BasicEnumPolyParams:
    def __init__(self, a_poly_size=3, b_poly_size=3, num_of_iterations=100, threshold=None, prec=100):
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
            range_a = [ [-enum_range, enum_range+1] for i in range(self.a_poly_size) ]
            range_b = [ [-enum_range, enum_range+1] for i in range(self.b_poly_size) ]

        pi_cont_frac = cont_fracs.PiContFrac([0 for i in range_a], [0 for i in range_b])
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

class MITM:
    def __init__(self, postproc_func=lambda x:x, a_poly_size=3, b_poly_size=3, num_of_iterations=100,
                 threshold=None, prec=50):
        self.bep = BasicEnumPolyParams(a_poly_size=a_poly_size, b_poly_size=b_poly_size,
                                       num_of_iterations=num_of_iterations, threshold=threshold, prec=prec)
        self.postproc_func = postproc_func
        self.dec_hashtable = DecimalHashTable(6)
        self.filtered_params = []

    def build_hashtable(self, enum_range=None, range_a=None, range_b=None, prec=None, ):
        pg = self.bep.pis_generator(enum_range=enum_range, range_a=range_a, range_b=range_b, prec=prec)
        self._iter2hashtalbe(pg)

    def _iter2hashtalbe(self, itr):
        for pi_cont_frac, pa, pb in itr:
            k = self.postproc_func(pi_cont_frac.get_pi())
            if k not in self.dec_hashtable:
                self.dec_hashtable[k] = []
            self.dec_hashtable[k].append((pa, pb))

    def find_clicks(self, u_range, l_range, c_range, d_range):
        if isinstance(u_range, int):
            u_range = range(-u_range, u_range+1)
        if isinstance(l_range, int):
            l_range = range(-l_range, l_range+1)
        if isinstance(c_range, int):
            c_range = range(-c_range, c_range+1)
        if isinstance(d_range, int):
            d_range = range(-d_range, d_range+1)
        real_pi = gen_real_pi()
        filtered_params = []

        for u,l,c,d in itertools.product(u_range, l_range, c_range, d_range):
            if d == 0 or l == 0:
                continue
            else:
                r = (u/real_pi + real_pi/l + c) / d
                if r in self.dec_hashtable:
                    filtered_params.extend([ (ab, (u,l,c,d)) for ab in self.dec_hashtable[r] ])

        self.filtered_params = filtered_params

    def refine_clicks(self, accuracy=10, num_of_iterations=3000, print_clicks=True):
        refined_params = []
        real_pi = gen_real_pi()
        # pi_cont_frac = cont_fracs.PiContFrac()
        for ab, ulcd in self.filtered_params:
            pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            # pi_cont_frac.reinitialize(a_coeffs=ab[0], b_coeffs=ab[1])
            pi_cont_frac.gen_iterations(num_of_iterations)
            u, l, c, d = ulcd
            rhs = self.postproc_func(pi_cont_frac.get_pi())
            lhs = (u/real_pi + real_pi/l + c) / d
            if self.compare_dec_with_accuracy(rhs, lhs, accuracy):
                if print_clicks
                    print(ab)
                    print(ulcd)
                    print(rhs)
                    print(lhs)
                    print('')
                refined_params.append((ab, ulcd))

        self.filtered_params = refined_params

    @staticmethod
    def build_pi_from_params(params):
        real_pi = gen_real_pi()
        ab, ulcd = params
        u,l,c,d = ulcd
        pa, pb = ab
        pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=pa, b_coeffs=pb)
        pi_cont_frac.gen_iterations(3000)
        pi_expression = (u/real_pi + real_pi/l + c) / d
        return (pi_cont_frac, pi_cont_frac.get_pi(), pi_expression)

    @staticmethod
    def compare_dec_with_accuracy(d1, d2, accuracy):
        # +1 for decimal dot
        accuracy += 1
        return d1.to_eng_string()[:accuracy+1] == d2.to_eng_string()[:accuracy+1]

