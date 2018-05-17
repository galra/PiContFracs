from basic_enum_params import BasicEnumPolyParams
import cont_fracs
from decimal import Decimal as dec
import itertools


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
        if len(pa1) > 1:
            return all([ ULCDEvaluator._is_equiv((((p1,), pb1), ulcd1_obj, post_func_ind1, convergence_info1),
                                               (((p2,), pb2), ulcd2_obj, post_func_ind2, convergence_info2))
                         for p1, p2 in zip(pa1, pa2)])
        if len(pb1) > 1:
            return all([ ULCDEvaluator._is_equiv(((pa1, (p1,)), ulcd1_obj, post_func_ind1, convergence_info1),
                                               ((pa2, (p2,)), ulcd2_obj, post_func_ind2, convergence_info2))
                         for p1, p2 in zip(pb1, pb2)])

        # if were're here, then params1 and params2 are single-element tuples
        pa1, pa2, pb1, pb2 = pa1[0], pa2[0], pb1[0], pb2[0]
        if pb1 == pb2:
            u1, l1, c1, d1 = ulcd1
            u2, l2, c2, d2 = ulcd2
            for s in [1, -1]:
                if u1 == u2 * s and l1 == l2 * s and c1 == c2 * s and d1 == -d2 * s:
                    return True
            if pa1 == pa2 or list(pa1) == [ -x for x in pa2 ]:
                ulcd_ratio = abs(d1 / d2)
                if u1 == u2 * ulcd_ratio and l1 == l2 / ulcd_ratio and c1 == c2 * ulcd_ratio:
                    return True
                if u1 == -u2 * ulcd_ratio and l1 == -l2 / ulcd_ratio and c1 == -c2 * ulcd_ratio:
                    return True
            if ulcd1 == ulcd2 and (pa1 == pa2 or list(pa1) == [ -x for x in pa2 ]):
                return True

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


class RationalFuncEvaluator:
    def __init__(self, numerator_coeffs, denominator_coeffs, target_constant):
        self.numerator_coeffs = numerator_coeffs
        self.denominator_coeffs = denominator_coeffs
        self.numerator = self.array_to_polynom(numerator_coeffs, target_constant)
        self.denominator = self.array_to_polynom(denominator_coeffs, target_constant)
        if not self.denominator.is_normal() or not self.numerator.is_normal():
            self.val = dec('nan')
        else:
            self.val = abs(self.numerator / self.denominator)

    def __str__(self):
        self.val.to_eng_string()

    def get_val(self):
        return self.val

    def get_params(self):
        return self.numerator_coeffs, self.denominator_coeffs

    @staticmethod
    def array_to_polynom(coeffs, x):
        return cont_fracs.ContFrac._array_to_polynom(coeffs, x)


class RationalFuncEnumerator(LHSEnumerator):
    def __init__(self, lhs_evaluator_params, target_value):
        super().__init__(lhs_evaluator_params, target_value)
        self._lhs_rational_numerator, self._lhs_rational_denominator = lhs_evaluator_params
        self.target_value = target_value
        self.lhs_rational_generator = BasicEnumPolyParams(avoid_int_roots=False, avoid_zero_b=False,
                                                          should_gen_contfrac=False, prec=self.prec)
        self._rational_generator, self._iter_len = lhs_rational_generator.polys_generator(range_a=[lhs_rational_numerator],
                                                                                      range_b=[lhs_rational_denominator])

    def generator(self):
        for i, (p_numerator, p_denominator) in enumerate(self.lhs_rational_generator):
            p_numerator = p_numerator[0]
            p_denominator = p_denominator[0]
            numerator = RationalFuncEvaluator.array_to_polynom(p_numerator, self.target_value)
            denominator = RationalFuncEvaluator.array_to_polynom(p_denominator, self.target_value)
            if not denominator.is_normal() or not numerator.is_normal():
                continue
            yield RationalFuncEvaluator(p_numerator, p_denominator)
