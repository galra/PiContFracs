from basic_enum_params import BasicEnumPolyParams
import cont_fracs
from decimal import Decimal as dec
import itertools
from functools import reduce
import operator
import math
# from flint import fmpz_poly as poly
import sympy
from postprocfuncs import POSTPROC_FUNCS

x_sym = sympy.symbols('x')
def poly(coeffs):
    return sympy.poly('+'.join('%d*x**%d' % (c, i) for i, c in enumerate(coeffs)), x_sym)

INVERSE_POSTPROC_PAIRS = [('safe_inverse', 'lambda x: x'),
                          ('lambda x: x**2', 'lambda x: safe_inverse(x**2)'),
                          ('lambda x: x**3', 'lambda x: safe_inverse(x**3)'),
                          ('lambda x: x**4', 'lambda x: safe_inverse(x**4)'),
                          ('lambda x: x**5', 'lambda x: safe_inverse(x**5)'),
                          ('safe_sqrt', 'lambda x: safe_inverse(safe_sqrt(x))')]
INVERSE_POSTPROC_PAIRS = [ (POSTPROC_FUNCS.index(f1), POSTPROC_FUNCS.index(f2)) for f1, f2 in INVERSE_POSTPROC_PAIRS ]
INVERSE_POSTPROC_PAIRS += [ funcs_pair[::-1] for funcs_pair in INVERSE_POSTPROC_PAIRS ]

class LHSEvaluator:
    def __init__(self, lhs_evaluator_params, target_constant=None):
        if isinstance(lhs_evaluator_params, LHSEvaluator):
            self.params, self.target_constant = lhs_evaluator_params.get_params(), lhs_evaluator_params.get_target_const()
        else:
            self.params = lhs_evaluator_params
            self.target_constant = target_constant
        self.reinit_params(self.params)

    def __str__(self):
        return 'Params %s, Res %s' % (str(self.params), self.val.to_eng_string())

    def __repr__(self):
        return self.__str__()

    def get_val(self):
        return self.val

    def get_params(self):
        return self.params

    def get_target_const(self):
        return self.target_constant

    def is_equiv(self, params):
        return False

    def get_latex_exp(self, target_constant_name):
        raise NotImplementedError('get_latex_exp is not implemeted!')

    def canonalize_params(self):
        pass


class LHSEnumerator:
    def __init__(self, params, target_constant=None):
        pass

    def __len__(self):
        return self._iter_len


class ULCDEvaluator(LHSEvaluator):
    def __init__(self, params, target_constant=None):
        super().__init__(params, target_constant)

    def reinit_params(self, params):
        self.params = params
        u, l, c, d = self.params
        self.val = (u / self.target_constant + self.target_constant / l + c) / d

    def flip_sign(self):
        u, l, c, d = self.params
        d *= -1
        self.params = (u, l, c, d)
        self.val = (u / self.target_constant + self.target_constant / l + c) / d

    def add_int(self, i):
        u, l, c, d = self.params
        c += d * i
        self.params = (u, l, c, d)
        self.val = (u / self.target_constant + self.target_constant / l + c) / d

    @staticmethod
    def is_equiv(params1, params2):
        ab1, ulcd1_obj, post_func_ind1, convergence_info1 = params1
        ab2, ulcd2_obj, post_func_ind2, convergence_info2 = params2
        if not isinstance(ulcd1_obj, type(ulcd2_obj)):
            return False
        ulcd1 = ulcd1_obj.get_params()
        ulcd2 = ulcd2_obj.get_params()

        if post_func_ind1 != 0 or post_func_ind2 != 0:
            return (ab1 == ab2 and ulcd1 == ulcd2 and post_func_ind1 == post_func_ind2)

        pa1, pb1 = ab1
        pa2, pb2 = ab2
        if len(pa1) != len(pa2) or len(pb1) != len(pb2):
            return False

        # if we're here, then params1 and params2 are single-element tuples
        if pb1 == pb2 and (pa1 == pa2 or pa1 == tuple(( tuple(( -k for k in p )) for p in pa2 ))):
            u1, l1, c1, d1 = ulcd1
            u2, l2, c2, d2 = ulcd2
            for s in [1, -1]:
                if u1 == u2 * s and l1 == l2 * s and c1 == c2 * s and abs(d1) == abs(d2):
                    return True
            ulcd_ratio = abs(d1 / d2)
            if u1 == u2 * ulcd_ratio and l1 == l2 / ulcd_ratio and c1 == c2 * ulcd_ratio:
                return True
            if u1 == -u2 * ulcd_ratio and l1 == -l2 / ulcd_ratio and c1 == -c2 * ulcd_ratio:
                return True
            # if ulcd1 == ulcd2 and (pa1 == pa2 or list(pa1) == [ -x for x in pa2 ]):
            #     return True

        return False

    def get_latex_exp(self, target_constant_name):
        latex_exp = []
        u, l, c, d = self.params

        if u != 0:
            latex_exp.append(r'\frac{{ {0} }}{{ {1} }}'.format(u, target_constant_name))
        if abs(l) != 1:
            latex_exp.append(r'\frac{{ {0} }}{{ {1} }}'.format(target_constant_name, l))
        elif l == -1:
            latex_exp.append(r'\left(-1\right) \cdot {0}'.format(target_constant_name))
        else:
            latex_exp.append(target_constant_name)
        if c != 0:
            latex_exp.append(str(c))
        latex_exp = '+'.join(latex_exp)
        if abs(d) != 1:
            latex_exp = r'\frac{{ 1 }}{{ {0} }} \left( {1} \right)'.format(d, latex_exp)
        elif d == -1:
            latex_exp = r'- \left( {1} \right)'.format(d, latex_exp)

        return latex_exp


class ULCDEnumerator(LHSEnumerator):
    def __init__(self, lhs_evaluator_params, target_constant):
        super().__init__(lhs_evaluator_params, target_constant)
        self.u_range, self.l_range, self.c_range, self.d_range = lhs_evaluator_params
        if isinstance(self.u_range, int):
            self.u_range = range(-self.u_range, self.u_range+1)
        if isinstance(self.l_range, int):
            self.l_range = range(-self.l_range, self.l_range+1)
        if isinstance(self.c_range, int):
            self.c_range = range(-self.c_range, self.c_range+1)
        if isinstance(self.d_range, int):
            self.d_range = range(1, self.d_range+1)
        self.target_value = target_constant

        self._iter_len = len(self.u_range) * len(self.l_range) * len(self.c_range) * len(self.d_range)

    def generator(self):
        ulcd_evaluator = ULCDEvaluator((1, 1, 1, 1), self.target_value)
        for u, l, c, d in itertools.product(self.u_range, self.l_range, self.c_range, self.d_range):
            if d == 0 or l == 0:
                continue
            elif abs(c) >= abs(d):
                continue
            ulcd_evaluator.reinit_params((u, l, c, d))
            yield ulcd_evaluator


class RationalFuncEvaluator(LHSEvaluator):
    def __init__(self, lhs_evaluator_params, target_constant=None):
        super().__init__(lhs_evaluator_params, target_constant)

    def reinit_params(self, params):
        self.numerator_p, self.denominator_p, self.added_int = params
        if not isinstance(self.added_int, int):
            raise TypeError('added_int has to be an integer')
        self._calc_val()

    def _calc_val(self):
        self.numerator = RationalFuncEvaluator.array_to_polynom([ int(i) for i in self.numerator_p.coeffs() ],
                                                                self.target_constant)
        self.denominator = RationalFuncEvaluator.array_to_polynom([ int(i) for i in self.denominator_p.coeffs() ],
                                                                  self.target_constant)

        # self.numerator = dec(str(self.numerator_p(self.target_constant))[2:-1])
        # self.denominator = dec(str(self.denominator_p(self.target_constant))[2:-1])


        if not self.denominator.is_normal() or not self.numerator.is_normal():
            self.val = dec('nan')
        else:
            # self.val = dec(self.numerator) / dec(self.denominator) + self.added_int
            self.val = self.numerator / self.denominator + self.added_int

    def flip_sign(self):
        self.denominator_p *= -1
        self.update_params()
        self._calc_val()

    def update_params(self):
        self.params = (self.numerator_p, self.denominator_p, self.added_int)
        self._calc_val()

    def add_int(self, n):
        self.added_int += n
        self.update_params()

    @staticmethod
    def is_equiv(params1, params2):
        ab1, ratio_func1_obj, post_func_ind1, convergence_info1 = params1
        ab2, ratio_func2_obj, post_func_ind2, convergence_info2 = params2
        if not isinstance(ratio_func1_obj, type(ratio_func2_obj)):
            return False

        # check if one is the inverse postproc func of the other, or if they are the same
        is_inverse = False
        if (post_func_ind1, post_func_ind2) in INVERSE_POSTPROC_PAIRS:
            is_inverse = True
        elif post_func_ind1 != post_func_ind2:
            return False

        numerator1_p, denominator1_p, _ = ratio_func1_obj.get_params()
        quotient1, rem1 = divmod(numerator1_p, denominator1_p)
        quotient1[0] = 0
        numerator2_p, denominator2_p, _ = ratio_func2_obj.get_params()
        if is_inverse:
            numerator2_p, denominator2_p = denominator2_p, numerator2_p
        quotient2, rem2 = divmod(numerator2_p, denominator2_p)
        quotient2[0] = 0
        # if ab1 == ab2 and f1/g1 = f2/g2 + C
        pa1, pb1 = ab1
        pa2, pb2 = ab2
        if pb1 == pb2 and (pa1 == pa2 or pa1 == tuple(( tuple(( -k for k in p )) for p in pa2 ))):
            if quotient1 == quotient2 and rem1 == rem2:
                return True
        return False
        # a1_ps, b1_ps = [ [ poly(list(k)) for k in cf_polys ] for cf_polys in ab1 ]
        # a2_ps, b2_ps = [ [ poly(list(k)) for k in cf_polys ] for cf_polys in ab2 ]
        # if (
        #             ((a1_ps == [ -p for p in a2_ps ]) and
        #              (all([ p.degree() < 2 for p in a1_ps ]) and all([ p.degree() < 2 for p in a2_ps ])) or
        #              (a1_ps == a2_ps)) and
        #             ((b1_ps == [ -p for p in b2_ps ]) and
        #              (all([ p.degree() < 2 for p in b1_ps ]) and all([ p.degree() < 2 for p in b2_ps ])) or
        #              (b1_ps == b2_ps))
        #    ):
        #     if ((numerator1_p == numerator2_p or numerator1_p == -numerator2_p) and
        #             (denominator1_p == denominator2_p or denominator1_p == -denominator2_p)):
        #         return True
        # return False

    def get_latex_exp(self, target_constant_name):
        latex_exp_parts = []
        for poly in [self.numerator_p, self.denominator_p]:
            poly_elements = []
            for i, c in enumerate([ int(i) for i in poly.coeffs() ]):
                if c == 0:
                    continue
                if i == 0:
                    poly_elements.append(str(c))
                    continue
                if i == 1:
                    if c == 1:
                        poly_elements.append(r'{0}'.format(target_constant_name))
                    elif c == -1:
                        poly_elements.append(r'-{0}'.format(target_constant_name))
                    else:
                        poly_elements.append(r'{0} {1}'.format(c, target_constant_name))
                    continue
                if c == 1:
                    poly_elements.append(r'{0}^{1}'.format(target_constant_name, i))
                else:
                    poly_elements.append(r'{0} {1}^{2}'.format(c, target_constant_name, i))
            poly_elements[0] = poly_elements[0].replace('{0}^{1}'.format(target_constant_name, 0), '')
            latex_exp_parts.append( '+'.join(poly_elements))
        if latex_exp_parts[1].strip() == '1':
            latex_exp = latex_exp_parts[0]
        else:
            latex_exp = r'\frac{{ {0} }}{{ {1} }}'.format(latex_exp_parts[0], latex_exp_parts[1])
        if self.added_int:
            latex_exp += '+ {0}'.format(self.added_int)
        return latex_exp

    def canonalize_params(self):
        numerator_p, denominator_p, added_int = self.params
        quotient, remainder = divmod(numerator_p, denominator_p)
        # params[2] = added int
        self.added_int += int(remainder[0])
        remainder[0] = 0
        # params[0] = numerator , params[1] = denominator
        self.old_numerator_denominator_p = (self.numerator_p, self.denominator_p)
        self.numerator_p = quotient * denominator_p + remainder
        self.denominator_p = denominator_p
        self.update_params()

    @staticmethod
    def array_to_polynom(coeffs, x):
        return cont_fracs.ContFrac._array_to_polynom(coeffs, x)

    # @staticmethod
    # def _normalize_poly(poly):
        # while poly and poly[-1] == 0:
        #     poly.pop()
        # if poly == []:
        #     poly.append(0)
    #
    # @staticmethod
    # def poly_divmod(num, den):
    #     #Create normalized copies of the args
    #     num = list(num[:])
    #     RationalFuncEvaluator._normalize_poly(num)
    #     den = list(den[:])
    #     RationalFuncEvaluator._normalize_poly(den)
    #
    #     if len(num) >= len(den):
    #         #Shift den towards right so it's the same degree as num
    #         shiftlen = len(num) - len(den)
    #         den = [0] * shiftlen + den
    #     else:
    #         return [0], num
    #
    #     quot = []
    #     divisor = float(den[-1])
    #     for i in range(shiftlen + 1):
    #         #Get the next coefficient of the quotient.
    #         mult = num[-1] / divisor
    #         quot = [mult] + quot
    #
    #         #Subtract mult * den from num, but don't bother if mult == 0
    #         #Note that when i==0, mult!=0; so quot is automatically normalized.
    #         if mult != 0:
    #             d = [mult * u for u in den]
    #             num = [u - v for u, v in zip(num, d)]
    #
    #         num.pop()
    #         den.pop(0)
    #
    #     RationalFuncEvaluator._normalize_poly(num)
    #     return quot, num


class RationalFuncEnumerator(LHSEnumerator):
    def __init__(self, lhs_evaluator_params, target_value):
        super().__init__(lhs_evaluator_params, target_value)
        self._lhs_rational_numerator, self._lhs_rational_denominator = lhs_evaluator_params
        self.target_value = target_value
        numerator_num_of_options = reduce(operator.mul, [ c[1] - c[0] for c in self._lhs_rational_numerator])
        denominator_num_of_options = reduce(operator.mul, [ c[1] - c[0] for c in self._lhs_rational_denominator])
        self._iter_len = numerator_num_of_options * denominator_num_of_options

    def generator(self):
        generated_funcs = []
        numerator_iterator = itertools.product(*[ range(*c_range) for c_range in self._lhs_rational_numerator ])
        for numerator_coeffs in numerator_iterator:
            if not any(numerator_coeffs):
                continue
            numerator_p = poly(list(numerator_coeffs))
            sign_coeff = 1
            if numerator_p.coeffs()[-1] < 0:
                numerator_p *= -1
                sign_coeff = -1
            denominator_iterator = itertools.product(*[ range(*c_range) for c_range in self._lhs_rational_denominator ])
            for denominator_coeffs in denominator_iterator:
                if not any(denominator_coeffs):
                    continue
                denominator_p = poly(list(denominator_coeffs))
                denominator_p *= sign_coeff
                # continue if a constant
                if (numerator_p.degree() < 2) and (denominator_p.degree() < 2):
                    continue
                # continue if zero or NaN numerator/denominator
                numerator_val = RationalFuncEvaluator.array_to_polynom(numerator_coeffs, self.target_value)
                denominator_val = RationalFuncEvaluator.array_to_polynom(denominator_coeffs, self.target_value)
                if not numerator_val.is_normal() or not denominator_val.is_normal():
                    continue
                # continue if they are linearly dependent
                quotient, rem = divmod(numerator_p, denominator_p)
                if not rem:
                    numerator_p = quotient
                    denominator_p = poly([1])
                else:
                    quotient, rem = divmod(denominator_p, numerator_p)
                    if not rem:
                        denominator_p = quotient
                        numerator_p = poly([1])

                gcd_p = numerator_p.gcd(denominator_p)
                if gcd_p != 1:
                    numerator_p = divmod(numerator_p, gcd_p)[0]
                    denominator_p = divmod(denominator_p, gcd_p)[0]
                func_params = (numerator_p, denominator_p)
                if (numerator_p.degree() < 2) and (denominator_p.degree() < 2):
                    continue
                if func_params in generated_funcs:
                    continue
                if numerator_p.degree() == -1 or denominator_p.degree() == -1:
                    continue
                generated_funcs.append(func_params)
                yield RationalFuncEvaluator((numerator_p, denominator_p, 0), self.target_value)

    # @staticmethod
    # def _are_polys_linearly_dependent(p1, p2):
    #     # 1st line: p1 is a monom, 2nd line p2 is a monom, 3rd line is linear dependency
    #     if (len(p1) == 1 or not any(p1[1:])) and (len(p2) == 1 or not any(p2[1:])) and \
    #             not all([p1[0] % p2[0], p2[0] % p1[0]]):
    #         return True
    #
    #     gcd_p1 = math.gcd(p1[0], p1[1])
    #     for c in p1[2:]:
    #         gcd_p1 = math.gcd(gcd_p1, c)
    #
    #     gcd_p2 = math.gcd(p2[0], p2[1])
    #     for c in p2[2:]:
    #         gcd_p2 = math.gcd(gcd_p2, c)
    #
    #     p1 = [ divmod(c, gcd_p1)[0] for c in p1 ]
    #     p2 = [ divmod(c, gcd_p2)[0] for c in p2 ]
    #     min_len = min(len(p1), len(p2))
    #     if any(p1[min_len:]) or any(p2[min_len:]):
    #         return False
    #     if all([ c1 == c2 for c1, c2 in zip(p1, p2) ]):
    #         return True
    #
    #     return False