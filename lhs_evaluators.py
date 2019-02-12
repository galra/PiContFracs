"""Implementation of enumeration over different LHS types, functions to be applied on the target constant."""

from decimal import Decimal as dec
import itertools
from functools import reduce
import operator
# from flint import fmpz_poly as poly
import sympy
from postprocfuncs import INVERSE_POSTPROC_PAIRS
from utils import MathOperations

x_sym = sympy.symbols('x')
def poly(coeffs):
    """Converts a polynomial of the form [a_0, a_1, a_2, ...] to  'a_0+a_1*x+a_2*x^2+...' sympy.poly object."""
    return sympy.poly('+'.join('%d*x**%d' % (c, i) for i, c in enumerate(coeffs)), x_sym)


class LHSEvaluator:
    """Evaluates a specific LHS with given parameters."""
    def __init__(self, lhs_evaluator_params, target_constant=None):
        """Initalizes. Calls the LHS type-specific class reinit_params method at the end. See their help for more info
        about the lhs_evaluator_params.
        Parameters:
            lhs_evaluator_params - either an LHSEvaluator instance to copy, or parameters for the specific type of LHS.
            target_constant - the sought target value (that will be substituted to the LHS equation)."""
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
        """Returns the evaluated value of the LHS."""
        return self.val

    def get_params(self):
        """Returns the LHS parameters (the structure differs between LHS types)."""
        return self.params

    def get_target_const(self):
        """Returns the sought target constant."""
        return self.target_constant

    def flip_sign(self):
        """Flips the sign of the LHS by changing the parameters appropriately."""
        raise NotImplementedError()

    def add_int(self, n):
        """Adds an integer to the LHS by changning the parameters appropriately."""
        raise NotImplementedError()

    def is_equiv(self, params):
        """Returns True if the LHS is equivalent to the LHS generated/represented by params."""
        return False

    def get_latex_exp(self, target_constant_name):
        """Returns a latex expression of the LHS.
        Parameters:
            target_constant_name - the name of the substituted target constant."""
        raise NotImplementedError('get_latex_exp is not implemeted!')

    def canonalize_params(self):
        """Turns the LHS parameters into canonical form."""
        raise NotImplementedError()


class LHSEnumerator:
    """Enumerates over different LHSs and create their instances as LHSEvaluator."""
    def __init__(self, params, target_constant=None):
        pass

    def __len__(self):
        return self._iter_len


class ULCDMetaClass(type):
    """Provides the value of str(ULCDClass) => 'ulcd'"""
    def __str__(self):
        return 'ulcd'


class ULCDEvaluator(LHSEvaluator, metaclass=ULCDMetaClass):
    """Evaluates a ULCD LHS, for which LHS_value = (u/target_val + target_val/l + c) / d"""
    def __init__(self, params, target_constant=None):
        """See LHSEvaluator.__init__"""
        super().__init__(params, target_constant)

    def reinit_params(self, params):
        """Reinitializes the parameters and recalculates the value.
        Parameters:
            params - (u, l, c, d)"""
        self.params = params
        u, l, c, d = self.params
        self.val = (u / self.target_constant + self.target_constant / l + c) / d

    def flip_sign(self):
        """Flips the sign of the LHS by changing the parameters appropriately."""
        u, l, c, d = self.params
        d *= -1
        self.params = (u, l, c, d)
        self.val = (u / self.target_constant + self.target_constant / l + c) / d

    def add_int(self, i):
        """Adds an integer to the LHS by changning the parameters appropriately."""
        u, l, c, d = self.params
        c += d * i
        self.params = (u, l, c, d)
        self.val = (u / self.target_constant + self.target_constant / l + c) / d

    @staticmethod
    def is_equiv(params1, params2):
        """Checks if params1 and params2 are equivalent.
        params = (ab, ulcd_obj, post_func_ind, convergence_info)
        ab = (a_poly, b_poly)

        Only non-interlace comparisons are supported (otherwise False is returned automatically).
        Checks for signs equivalences and if ulcd params are redundant, i.e. (u, l, c, d) and (2u, l/2, 2c, 2d)"""
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
        """Returns a latex expression of the LHS.
        Parameters:
            target_constant_name - the name of the substituted target constant."""
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


class ULCDEnumerator(LHSEnumerator, metaclass=ULCDMetaClass):
    """Enumerates over ULCD LHSs. See help for ULCDEvaluator for more info about ULCD."""
    def __init__(self, lhs_evaluator_params, target_constant):
        """Initializes.
        Parameters:
            lhs_evaluator_params - a range for each of the parameters u, l, c, d in the format:
                        [u_range, l_range, c_range, d_range]
                        where each of the u, l, c ranges is in one of the formats (e.g. u_range):
                            u_range = n ==> the range is -n to n (including end points)
                            u_range = iter ==> the options for u are those that iter results with. iter must have len().
                            u_range = [u_min, u_max, <u_step>] ==> the range is u_min to u_max-1 <with u_step steps>.
                        for d the options are the same, with the difference that:
                            d_range = n ==> the range is 1 to d (including end points)
            target_constant - the sought constant to be substitued into the ULCD expression"""
        super().__init__(lhs_evaluator_params, target_constant)
        self.u_range, self.l_range, self.c_range, self.d_range = lhs_evaluator_params
        if isinstance(self.u_range, int):
            self.u_range = range(-self.u_range, self.u_range+1)
        elif (isinstance(self.u_range, list) or isinstance(self.u_range, tuple)) and len(self.u_range) in [2, 3]:
            self.u_range = range(*self.u_range)
        if isinstance(self.l_range, int):
            self.l_range = range(-self.l_range, self.l_range+1)
        elif (isinstance(self.l_range, list) or isinstance(self.l_range, tuple)) and len(self.l_range) in [2, 3]:
            self.l_range = range(*self.l_range)
        if isinstance(self.c_range, int):
            self.c_range = range(-self.c_range, self.c_range+1)
        elif (isinstance(self.c_range, list) or isinstance(self.c_range, tuple)) and len(self.c_range) in [2, 3]:
            self.c_range = range(*self.c_range)
        if isinstance(self.d_range, int):
            self.d_range = range(1, self.d_range+1)
        elif (isinstance(self.d_range, list) or isinstance(self.d_range, tuple)) and len(self.d_range) in [2, 3]:
            self.d_range = range(*self.d_range)
        self.target_value = target_constant

        self._iter_len = len(self.u_range) * len(self.l_range) * len(self.c_range) * len(self.d_range)

    def generator(self):
        """Generates all the ULCD LHSs (as ULCDEvaluator), from the cartesian product of the parameters ranges.
        The generator assumes that the results are EQUIVALENT UP TO AN INTEGER CONSTANT, therefore different LHSs with
        a difference of an integer are skipped!"""
        # Creates a single instance of ULCDEvaluator that will be reinitialized to save CPU time.
        ulcd_evaluator = ULCDEvaluator((1, 1, 1, 1), self.target_value)
        for u, l, c, d in itertools.product(self.u_range, self.l_range, self.c_range, self.d_range):
            # skip illegal options
            if d == 0 or l == 0:
                continue
            # Assuming that c range includes small numbers, and that the results are compared up to a constant, skip
            # instances for which abs(c/d) > 1
            elif abs(c) >= abs(d):
                continue
            ulcd_evaluator.reinit_params((u, l, c, d))
            yield ulcd_evaluator


class RationalFuncMetaClass(type):
    """Provides the value of str(ULCDClass) => 'ulcd'"""
    def __str__(self):
        return 'rationalfunc'


class RationalFuncEvaluator(LHSEvaluator, metaclass=RationalFuncMetaClass):
    def __init__(self, lhs_evaluator_params, target_constant=None):
        """See LHSEvaluator.__init__"""
        super().__init__(lhs_evaluator_params, target_constant)

    def reinit_params(self, params):
        """Initializes and calculates the value.
        Parameters:
            params - (numerator_poly, denominator_poly, added_int) the added int is an integer to be added to the result
                     and the polynomials are of the format: [a_0, a_1, a_2, ...] <=> a_0 + a_1*x + a_2*x^2 + ..."""
        # WARNING: as far as I see, only the first two arguments will be fed to params from the config file. Bug?
        self.numerator_p, self.denominator_p, self.added_int = params
        if not isinstance(self.added_int, int):
            raise TypeError('added_int has to be an integer')
        self._calc_val()

    def _calc_val(self):
        """Calculates the value of the LHS based on the LHS parameters (numerator, denominator, target constant and
        added int)."""
        self.numerator = MathOperations.subs_in_polynom([ int(i) for i in self.numerator_p.coeffs() ],
                                                        self.target_constant)
        self.denominator = MathOperations.subs_in_polynom([ int(i) for i in self.denominator_p.coeffs() ],
                                                          self.target_constant)

        # self.numerator = dec(str(self.numerator_p(self.target_constant))[2:-1])
        # self.denominator = dec(str(self.denominator_p(self.target_constant))[2:-1])

        if not self.denominator.is_normal() or not self.numerator.is_normal():
            self.val = dec('nan')
        else:
            # self.val = dec(self.numerator) / dec(self.denominator) + self.added_int
            self.val = self.numerator / self.denominator + self.added_int

    def flip_sign(self):
        """Flips the sign of the LHS by changing the parameters appropriately."""
        self.denominator_p *= -1
        self.update_params()
        self._calc_val()

    def add_int(self, n):
        """Adds an integer to the LHS by changning the parameters appropriately."""
        self.added_int += n
        self.update_params()

    def update_params(self):
        """Updates the self.params to the new, current parameters (numerator, denominator and added_int)
        and recalculates the LHS value."""
        self.params = (self.numerator_p, self.denominator_p, self.added_int)
        self._calc_val()

    @staticmethod
    def is_equiv(params1, params2):
        """Checks if params1 and params2 are equivalent.
        params = (ab, ratio_func_obj, post_func_ind, convergence_info)
        ab = (a_poly, b_poly)

        Validates that both are RationalFuncEvaluator, with the same postproc_func or inverse (1/...) ones, and checks
        if p1/q1 == p2/q2 up to a constant, and a1 == +-a2, b1 == b2"""
        ab1, ratio_func1_obj, post_func_ind1, convergence_info1 = params1
        ab2, ratio_func2_obj, post_func_ind2, convergence_info2 = params2
        # Make sure both LHSs are RationalFuncEvaluator instances
        if not isinstance(ratio_func1_obj, RationalFuncEvaluator) or \
                not isinstance(ratio_func2_obj, RationalFuncEvaluator):
            return False

        # Validate that the postproc funcs are the same, or one is the inverse (1/...) of the other
        is_inverse = False
        if (post_func_ind1, post_func_ind2) in INVERSE_POSTPROC_PAIRS:
            is_inverse = True
        elif post_func_ind1 != post_func_ind2:
            return False

        # find the reduced form, and check if p1/q1==p2/q2 in the reduced form up to a constant, and if b1==b2, a1==+-a2
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
        """Returns a latex expression of the LHS.
        Parameters:
            target_constant_name - the name of the substituted target constant."""
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
        """Turns the LHS parameters into canonical form by moving from p=q*k+r <=> p/q=k+r/q to the same expression,
        with the free coefficient k_0 set to 0.
        A reduction to (p/gcd(p, q)) / (q/gcd(p, q)) isn't dont because p, q are expected to be supplied in reduced form
        in advance, and so RationalFuncEnumerator does."""
        # TODO: fix the above documentation.
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

    # Old method. Left as a backup. Delete if not used in a while
    # 30/01/2019
    # TODO: Delete if needed
    # @staticmethod
    # def array_to_polynom(coeffs, x):
    #     """Syntax sugar (to avoid VERY long expressions) to call _array_to_polynom). For more info about the parameters,
    #     check the help for cont_fracs.ContFrac._array_to_polynom."""
    #     return cont_fracs.ContFrac._array_to_polynom(coeffs, x)

    # Old methods for lists polynomials instead of sympy.Poly. Delete if it's not needed in a while.
    # 30/01/2019
    # TODO: Delete if needed
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


class RationalFuncEnumerator(LHSEnumerator, metaclass=RationalFuncMetaClass):
    """Enumerates over rational functions."""
    def __init__(self, ratiofunc_evaluator_params, target_constant):
        """Initializes.
        Parameters:
            ratiofunc_evaluator_params - a range for the parameters of the numerator/denominator polynomials:
                [ [[p_0_min, p_0_max], [p_1_min, p_1_max], ...], [[q_0_min, p_0_max], [q_1_min, q_1_max], ...] ]
                where p(x) = p_0 + p_1*x + p_2*x^2 + ... and similarly for q(x)=q_0 + ...
                The range is from (e.g. for p_0) p_0_min to p_1_max-1.
            target_constant - the sought constant to be substitued into the ULCD expression."""
        super().__init__(ratiofunc_evaluator_params, target_constant)
        self._lhs_rational_numerator, self._lhs_rational_denominator = ratiofunc_evaluator_params
        self.target_constant = target_constant
        numerator_num_of_options = reduce(operator.mul, [ c[1] - c[0] for c in self._lhs_rational_numerator])
        denominator_num_of_options = reduce(operator.mul, [ c[1] - c[0] for c in self._lhs_rational_denominator])
        self._iter_len = numerator_num_of_options * denominator_num_of_options

    def generator(self):
        """Generates rational functions as instances of RationalFuncEvaluator from the parameters range supplied during
        initialization.
        """
        generated_funcs = []
        numerator_iterator = itertools.product(*[ range(*c_range) for c_range in self._lhs_rational_numerator ])
        for numerator_poly_coeffs in numerator_iterator:
            # skip the zero polynomial
            if not any(numerator_poly_coeffs):
                continue
            numerator_p = poly(list(numerator_poly_coeffs))
            sign_coeff = 1
            # switch the numerator to a 'fixed' polynomial, with positive leading coefficient.
            if numerator_p.coeffs()[-1] < 0:
                numerator_p *= -1
                sign_coeff = -1
            denominator_iterator = itertools.product(*[ range(*c_range) for c_range in self._lhs_rational_denominator ])
            for denominator_poly_coeffs in denominator_iterator:
                # skip the zero polynomial
                if not any(denominator_poly_coeffs):
                    continue
                denominator_p = poly(list(denominator_poly_coeffs))
                # switch sign if p had to be 'fixed'
                denominator_p *= sign_coeff
                # Skip the constant functions p_0/q_0
                if (numerator_p.degree() < 2) and (denominator_p.degree() < 2):
                    continue
                # continue if zero or NaN numerator/denominator
                numerator_val = MathOperations.subs_in_polynom(numerator_poly_coeffs, self.target_constant)
                denominator_val = MathOperations.subs_in_polynom(denominator_poly_coeffs, self.target_constant)
                # normal: finite, non-zero, not ridiculously small (~E-9999999999999...) and not NaN
                if not numerator_val.is_normal() or not denominator_val.is_normal():
                    continue

                # Simplify p, q if they aren't mutually prime
                # Try to divide them by one another.
                # Redundant because the gcd is calculated right after. Left here in case problems will arise.
                # Please delete later if note.
                # 30/01/2019
                # TODO: delete if needed.
                # quotient, rem = divmod(numerator_p, denominator_p)
                # if not rem:
                #     numerator_p = quotient
                #     denominator_p = poly([1])
                # else:
                #     quotient, rem = divmod(denominator_p, numerator_p)
                #     if not rem:
                #         denominator_p = quotient
                #         numerator_p = poly([1])
                #
                gcd_p = numerator_p.gcd(denominator_p)
                if gcd_p != 1:
                    numerator_p = divmod(numerator_p, gcd_p)[0]
                    denominator_p = divmod(denominator_p, gcd_p)[0]
                func_params = (numerator_p, denominator_p)
                # Skip the constant functions p_0/q_0
                if (numerator_p.degree() < 2) and (denominator_p.degree() < 2):
                    continue
                # If an equivalent function was tested already - continue
                if func_params in generated_funcs:
                    continue
                # Unclear what this is for. The zero polynomial has been filtered out at the beginning, and has a degree
                # of sympy.numbers.NegativeInfinity().
                # This is left here in case it's purpose will be cleared in the near future. Otherwise, please remove it
                # 30/01/2019
                # TODO: remove if needed
                # if numerator_p.degree() == -1 or denominator_p.degree() == -1:
                #     continue
                generated_funcs.append(func_params)
                yield RationalFuncEvaluator((numerator_p, denominator_p, 0), self.target_constant)

    # Old method, for polynomials that are lists instead of sympy.Poly
    # Delete if it's not needed in a while
    # 30/01/2019
    # TODO: Delete if needed
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