import itertools
from functools import wraps
from basic_algo import set_precision
import cont_fracs
from decimal import Decimal as dec


def _len_decorator(func):
    @wraps(func)
    def wrapper(self, enum_range=None, range_a=None, range_b=None, prec=None):
        range_a, range_b = self._convert_ranges_to_range_a_range_b(enum_range, range_a, range_b)
        gen_len = 1
        for ar in range_a:
            for r in ar:
                gen_len *= (r[1] - r[0])
        for br in range_b:
            for r in br:
                gen_len *= (r[1] - r[0])

        return func(self, enum_range, range_a, range_b, prec), gen_len
    return wrapper


class BasicEnumPolyParams:
    def __init__(self, a_poly_size=3, b_poly_size=3, num_of_a_polys=1, num_of_b_polys=1, num_of_iterations=100,
                 enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True, avoid_zero_b=True,
                 threshold=None, prec=100):
        self._a_poly_size = a_poly_size
        self._b_poly_size = b_poly_size
        self._num_of_a_polys = num_of_a_polys
        self._num_of_b_polys = num_of_b_polys
        self.num_of_iterations = num_of_iterations
        self._enum_only_exp_conv = enum_only_exp_conv
        self._avoid_int_roots = avoid_int_roots
        self._should_gen_contfrac = should_gen_contfrac
        self._avoid_zero_b = avoid_zero_b
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

    def _convert_ranges_to_range_a_range_b(self, enum_range, range_a, range_b):
        if not (enum_range or range_a or range_b):
            raise ValueError("No range was given")
        if isinstance(enum_range, list):
            if not range_a:
                # range_a = [enum_range for i in range(self.a_poly_size) ]
                range_a = [[enum_range] * self._a_poly_size] * self._num_of_a_polys
            if not range_b:
                # range_b = [enum_range for i in range(self.b_poly_size) ]
                range_b = [[enum_range] * self._b_poly_size] * self._num_of_b_polys
        elif enum_range:
            if not range_a:
                range_a = [ [[-enum_range, enum_range+1] for i in range(self._a_poly_size)] ] * self._num_of_a_polys
            if not range_b:
                range_b = [ [[-enum_range, enum_range+1] for i in range(self._b_poly_size)] ] * self._num_of_b_polys

        return range_a, range_b

    @_len_decorator
    def polys_generator(self, enum_range=None, range_a=None, range_b=None, prec=None):
        """enum_range - a number for the range [-enum_range, enum_range] or a specific range of the format [first, last+1].
range_a/range_b - should be of the format [first, last+1].
    if two switching polynomials are used, this should be [[first1, last1+1], [first2, last2+1]]"""
        range_a, range_b = self._convert_ranges_to_range_a_range_b(enum_range, range_a, range_b)
        # else:
        #     if self._num_of_a_polys == 1:
        #         range_a = [range_a]
        #     if self._num_of_b_polys == 1:
        #         range_b = [range_b]

        # self.time_measure = 0
        cont_frac = cont_fracs.ContFrac([0 for i in range_a], [0 for i in range_b])
        a_params_iterator = itertools.product(*[ itertools.product(*[ range(*r) for r in ra ]) for ra in range_a ])
        for pas in a_params_iterator:
            b_params_iterator = itertools.product(*[itertools.product(*[ range(*r) for r in rb ]) for rb in range_b ])
            for pbs in b_params_iterator:
                for pb in pbs:
                    if self._avoid_int_roots and self._does_have_integer_roots(pb):
                        continue
                if len(pas) > 1 and all([ pas[0] == p for p in pas ]):
                    pas = (pas[0],)
                if len(pbs) > 1 and all([ pbs[0] == p for p in pbs ]):
                    pbs = (pbs[0],)
                if len(pas) == 1 and len(pbs) == 1 and self._enum_only_exp_conv and \
                   self._polynom_degree(pbs[0]) > 2 * self._polynom_degree(pas[0]):
                    continue
                if self._should_gen_contfrac:
                    cont_frac.reinitialize(pas, pbs)
                    if self._avoid_zero_b:
                        try:
                            cont_frac.gen_iterations(self.num_of_iterations)
                        except cont_fracs.ZeroB:
                            continue
                    yield (cont_frac, pas, pbs)
                else:
                    yield(pas, pbs)

    @staticmethod
    def _does_have_integer_roots(poly):
        mask = [ i for i,e in enumerate(poly) if e != 0 ]
        if len(mask) == 0:
            return True
        elif len(mask) == 1:
            return False
        elif len(mask) == 2 and mask[1] == mask[0] + 1:
            return True
        if len(poly) == 3:
            if poly[2] == 0:
                if (poly[1] / poly[0]).is_integer():
                    return True
            else:
                disc = (poly[1]**2-4*poly[0]*poly[2])
                if disc < 0:
                    return False
                disc **= 0.5
                if ((disc - poly[1])/(2*poly[0])).is_integer():
                    return True
                if ((-disc - poly[1])/(2*poly[0])).is_integer():
                    return True
        return False

    @staticmethod
    def _polynom_degree(p):
        deg = len(p)
        for i in p[::-1]:
            if i != 0:
                break
            deg -= 1
        return deg

    def enum_params(self, enum_range=None, range_a=None, range_b=None, show_progress=True):
        params_iter = self.polys_generator(enum_range, range_a, range_b)
        iter_num = 1
        for cont_frac, pa, pb in params_iter:
            if cont_frac.is_result_valid() and cont_frac.compare_result() < self.threshold:
                    self.good_params.append({'a': pa, 'b': pb})
            if iter_num % 100 == 0 and show_progress:
                print('\r%d' % (iter_num), end='')
            iter_num += 1
        if show_progress:
            print('')