import itertools
from functools import wraps
import cont_fracs
from decimal import Decimal as dec
import decimal
from scipy.special import binom

def set_precision(prec):
    """Sets the precision of the current decimal context"""
    decimal.getcontext().prec=prec

def _len_decorator(func):
    """Used as a decorator to estimated the length of a generator with signature:
    gen(range_a, range_b, prec)
    Instead of returning a new generator 'gen',
    the tuple (gen, estimated_len(gen)) is returned."""
    @wraps(func)
    def wrapper(self, range_a=None, range_b=None, prec=None):
        gen_len = 1
        for ar in range_a:
            for r in ar:
                gen_len *= (r[1] - r[0])
        for br in range_b:
            for r in br:
                gen_len *= (r[1] - r[0])

        return func(self, enum_range, range_a, range_b, prec), gen_len
    return wrapper


class NormalMetaClass(type):
    """Provides the value of str(BasicEnumPolyParams) => 'normal'"""
    def __str__(self):
        return 'normal'


class BasicEnumPolyParams(metaclass=NormalMetaClass):
    """Enumerates over a,b polynomials for continued fractions generation. The generated polynomials are taken as is,
    unlike the non-basic enum poly params.
    Other classes may inherit from this one and implement the static method manipulate_poly as a generator, to generate
    from each "regular polynom" the required polynomials (or yield a single transformed polynomial)."""
    def __init__(self, num_of_iterations=300, enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True,
                 avoid_zero_b=True, threshold=None, prec=80):
        """
        init
        :param self: self
        :param num_of_iterations: if should_gen_contfrac=True, then how many iterations should be promoted.
        :param enum_only_exp_conv: skip contfracs that surely won't converge (at least) exponentially
        :param avoid_int_roots: avoid contfracs for which b (or any of its interlaces) have an integer root.
                                it will affect only b polynomial of degree<=3.
                                notice that if interlace is used, this may skip valid contfracs (the i root value may be
                                assigned to another interlace polynomial instead of the one with the root).
        :param should_gen_contfrac: generate ContFrac object and return (cont_frac, pas, pbs) instead of (pas, pbs)
        :param avoid_zero_b: raise exception ZeroB and cancel if finite contfrac (b_i=0 for some i) is achieved
        :param threshold: threshold for considering contfrac_res == target_val
                          if abs(contfrac_res-target_val) < threshold.
        :param prec: decimal precision to be used for calculations
        :return: nothing
        """
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

    @_len_decorator
    def polys_generator(self, range_a=None, range_b=None):
        """range_a/range_b - for example: [ [[], []], [[], [], []], [[], [], [], []], [[], []] ] is a 4-interlace with
                             degrees of 2,3,4,2 . [m n] means running on coefficients between m to n-1."""
        cont_frac = cont_fracs.ContFrac([0]*range_a, [0]*range_b, avoid_zero_b=self._avoid_zero_b)
        a_params_iterator = itertools.product(*[ itertools.product(*[ range(*r) for r in ra ]) for ra in range_a ])

        for pas_premanipulate in a_params_iterator:
            pas_manipulated_gen = self.manipulate_poly(pas_premanipulate)
            for pas in pas_manipulated_gen:
                b_params_iterator = itertools.product(*[itertools.product(*[ range(*r) for r in rb ]) for rb in range_b ])
                for pbs_premanipulate in b_params_iterator:
                    pbs_manipulated_gen = self.manipulate_poly(pbs_premanipulate)
                    for pbs in pbs_manipulated_gen:
                        for pb in pbs:
                            if self._avoid_int_roots and self._does_have_integer_roots(pb):
                                continue
                        # in the case of an interlace in which all the interlace-polynomials are identical,
                        # squeeze it to a single polynomial with no interlace
                        if len(pas) > 1 and all([ pas[0] == p for p in pas ]):
                            pas = (pas[0],)
                        if len(pbs) > 1 and all([ pbs[0] == p for p in pbs ]):
                            pbs = (pbs[0],)
                        # if no interlace and only exponential convergence contfracs should be enumerated, make sure
                        # that 2*deg(a) >= deg(b)
                        if len(pas) == 1 and len(pbs) == 1 and self._enum_only_exp_conv and \
                           self._polynom_degree(pbs[0]) > 2 * self._polynom_degree(pas[0]):
                            continue
                        # generate contfrac and return everything / return the polynomials
                        if self._should_gen_contfrac:
                            cont_frac.reinitialize(pas, pbs)
                            try:
                                cont_frac.gen_iterations(self.num_of_iterations)
                            except cont_fracs.ZeroB:
                                continue
                            yield (cont_frac, pas, pbs)
                        else:
                            yield(pas, pbs)

    @staticmethod
    def _does_have_integer_roots(poly):
        """For a poly of deg(poly)<=3, check if it has integer roots. Returns a boolean."""
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
        """Finds deg(p). It may be different than len(p). E.g. deg([1, 1, 0, 0]) is 1 and not 3."""
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

    @staticmethod
    def manipulate_poly(poly):
        yield poly


class IndexedMetaClass(NormalMetaClass):
    """Provides the value of str(IndexedParameterEnumPolyParams) => 'indexed'"""
    def __str__(self):
        return 'indexed'


class IndexedParameterEnumPolyParams(BasicEnumPolyParams, metaclass=IndexedMetaClass):
    def __init__(self, a_poly_size=3, b_poly_size=3, num_of_a_polys=1, num_of_b_polys=1, num_of_iterations=300,
                 enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True, avoid_zero_b=True,
                 threshold=None, prec=80, special_params=None):
        super().__init__(a_poly_size=a_poly_size, b_poly_size=b_poly_size, num_of_a_polys=num_of_a_polys,
                         num_of_b_polys=num_of_b_polys, num_of_iterations=num_of_iterations,
                         enum_only_exp_conv=enum_only_exp_conv, avoid_int_roots=avoid_int_roots,
                         should_gen_contfrac=should_gen_contfrac, avoid_zero_b=avoid_zero_b,
                         threshold=threshold, prec=prec)

    @staticmethod
    def manipulate_poly(poly):
        poly_new = []
        for p in poly:
            p_new = [ 0 ] * len(p)
            deg_p = len(p)-1
            # p_new += sum_i (x+i)**deg_p
            for i, a_i in enumerate(p):
                if a_i == 0:
                    continue
                # p_new += (x+i)**deg_p
                for j in range(deg_p+1):
                    p_new[j] += a_i * int(binom(deg_p, j)) * i**j
            poly_new.append(p_new)
        yield poly_new


class SparseMetaClass(NormalMetaClass):
    """Provides the value of str(SparseParameterEnumPolyParams) => 'sparse'"""
    def __str__(self):
        return 'sparse'


class SparseParameterEnumPolyParams(BasicEnumPolyParams, metaclass=IndexedMetaClass):
    def __init__(self, a_poly_size=3, b_poly_size=3, num_of_a_polys=1, num_of_b_polys=1, num_of_iterations=300,
                 enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True, avoid_zero_b=True,
                 threshold=None, prec=80, special_params=None):
        # for the calculation of the sparse polynom by "n over k" options
        self.n, self.k = special_params
        super().__init__(a_poly_size=a_poly_size, b_poly_size=b_poly_size, num_of_a_polys=num_of_a_polys,
                         num_of_b_polys=num_of_b_polys, num_of_iterations=num_of_iterations,
                         enum_only_exp_conv=enum_only_exp_conv, avoid_int_roots=avoid_int_roots,
                         should_gen_contfrac=should_gen_contfrac, avoid_zero_b=avoid_zero_b,
                         threshold=threshold, prec=prec)

    def manipulate_poly(poly):
        for poly_template in itertools.combinations(range(self.n), self.k):
            yield self._apply_poly_template_on_coeffs_list(poly, )

    def _apply_poly_template_on_coeffs_list(self, coeffs_list, poly_template):
        if len(coeffs_list) > poly_template[-1]:
            print('Warning: coeffs_list=%s\npoly_template=%s' % (str(coeffs_list), str(poly_template)))
        res_poly = [ 0 ] * poly_template[-1]
        for i,v in zip(poly_template, coeffs_list):
            res_poly[i] = v
