from basic_algo import set_precision
import cont_fracs
from decimal_hashtable import DecimalHashTable
from decimal import Decimal as dec
import itertools
from gen_real_consts import gen_real_pi
import csv
import sys

class BasicEnumPolyParams:
    def __init__(self, a_poly_size=3, b_poly_size=3, num_of_a_polys=1, num_of_b_polys=1, num_of_iterations=100,
                 enum_only_exp_conv = False, threshold=None, prec=100):
        self._a_poly_size = a_poly_size
        self._b_poly_size = b_poly_size
        self._num_of_a_polys = num_of_a_polys
        self._num_of_b_polys = num_of_b_polys
        self.num_of_iterations = num_of_iterations
        self._enum_only_exp_conv = enum_only_exp_conv
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
        """enum_range - a number for the range [-enum_range, enum_range] or a specific range of the format [first, last+1].
range_a/range_b - should be of the format [first, last+1].
    if two switching polynomials are used, this should be [[first1, last1+1], [first2, last2+1]]"""
        if not (enum_range or range_a or range_b):
            raise ValueError("No range was given")
        if isinstance(enum_range, list):
            # range_a = [enum_range for i in range(self.a_poly_size) ]
            range_a = [[enum_range] * self._a_poly_size] * self._num_of_a_polys
            # range_b = [enum_range for i in range(self.b_poly_size) ]
            range_b = [[enum_range] * self._b_poly_size] * self._num_of_b_polys
        elif enum_range:
            range_a = [ [[-enum_range, enum_range+1] for i in range(self._a_poly_size)] ] * self._num_of_a_polys
            range_b = [ [[-enum_range, enum_range+1] for i in range(self._b_poly_size)] ] * self._num_of_b_polys
        else:
            if self._num_of_a_polys == 1:
                range_a = [range_a]
            if self._num_of_b_polys == 1:
                range_b = [range_b]

        pi_cont_frac = cont_fracs.PiContFrac([0 for i in range_a], [0 for i in range_b])
        a_params_iterator = itertools.product(*[ itertools.product(*[ range(*r) for r in ra ]) for ra in range_a ])
        for pas in a_params_iterator:
            b_params_iterator = itertools.product(*[itertools.product(*[ range(*r) for r in rb ]) for rb in range_b ])
            for pbs in b_params_iterator:
                for pb in pbs:
                    if self._does_have_integer_roots(pb):
                        continue
                pi_cont_frac.reinitialize(pas, pbs)
                try:
                    pi_cont_frac.gen_iterations(self.num_of_iterations)
                except cont_fracs.ZeroB:
                    continue
                if len(pas) == 2 and pas[0] == pas[1]:
                    pas = (pas[0],)
                if len(pbs) == 2 and pbs[0] == pbs[1]:
                    pbs = (pbs[0],)
                if len(pas) == 1 and len(pbs) == 1 and self._polynom_degree(pbs[0]) > 2 * self._polynom_degree(pas[0]):
                    continue
                yield (pi_cont_frac, pas, pbs)

    def _does_have_integer_roots(self, poly):
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

    def _polynom_degree(self, p):
        deg = len(p)
        for i in p[::-1]:
            if i != 0:
                break
            deg -= 1
        return deg

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
    def __init__(self, target_generator=gen_real_pi, postproc_funcs=[lambda x:x], trunc_integer=True, a_poly_size=3,
                 b_poly_size=3, num_of_a_polys=1, num_of_b_polys=1, enum_only_exp_conv=True, num_of_iterations=100,
                 threshold=None, prec=50):
        self.bep = BasicEnumPolyParams(a_poly_size=a_poly_size, b_poly_size=b_poly_size, num_of_a_polys=num_of_a_polys,
                                       num_of_b_polys=num_of_b_polys, enum_only_exp_conv=enum_only_exp_conv,
                                       num_of_iterations=num_of_iterations, threshold=threshold, prec=prec)
        self.target_generator = target_generator
        self.postproc_funcs = postproc_funcs
        self.trunc_integer = trunc_integer
        self.dec_hashtable = DecimalHashTable(6)
        self.filtered_params = []
        self.uniq_params = []

    def build_hashtable(self, enum_range=None, range_a=None, range_b=None, prec=None):
        pg = self.bep.pis_generator(enum_range=enum_range, range_a=range_a, range_b=range_b, prec=prec)
        self._iter2hashtalbe(pg)

    def _iter2hashtalbe(self, itr):
        for pi_cont_frac, pa, pb in itr:
            cur_cont_frac = pi_cont_frac.get_pi()
            if not cur_cont_frac.is_normal():
                continue
            for post_func_ind, post_f in enumerate(self.postproc_funcs):
                k = post_f(cur_cont_frac)
                if not k.is_normal():
                    continue
                k = abs(k)
                if self.trunc_integer:
                    k -= int(k)
                if k not in self.dec_hashtable:
                    self.dec_hashtable[k] = []
                self.dec_hashtable[k].append(((pa, pb), post_func_ind))

    def find_clicks(self, u_range, l_range, c_range, d_range):
        if isinstance(u_range, int):
            u_range = range(-u_range, u_range+1)
        if isinstance(l_range, int):
            l_range = range(-l_range, l_range+1)
        if isinstance(c_range, int):
            c_range = range(-c_range, c_range+1)
        if isinstance(d_range, int):
            d_range = range(1, d_range+1)
        target_value = self.target_generator()
        filtered_params = []

        for u,l,c,d in itertools.product(u_range, l_range, c_range, d_range):
            if d == 0 or l == 0:
                continue
            elif abs(c) >= abs(d):
                continue
            else:
                r = abs((u/target_value + target_value/l + c) / d)
                if self.trunc_integer:
                    r -= int(r)
                if r in self.dec_hashtable:
                    filtered_params.extend([ (ab, (u,l,c,d), post_func_ind, None) for ab, post_func_ind in self.dec_hashtable[r] ])

        self.filtered_params = filtered_params

    def refine_clicks(self, accuracy=10, num_of_iterations=3000, print_clicks=False):
        refined_params = []
        target_value = self.target_generator()
        # pi_cont_frac = cont_fracs.PiContFrac()
        for ab, ulcd, post_func_ind, convergence_info in self.filtered_params:
            pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            pi_cont_frac.set_approach_type_and_params(convergence_info)
            # pi_cont_frac.reinitialize(a_coeffs=ab[0], b_coeffs=ab[1])
            pi_cont_frac.gen_iterations(num_of_iterations, dec('1E-%d' % (accuracy+3)))
            u, l, c, d = ulcd
            signed_rhs = self.postproc_funcs[post_func_ind](pi_cont_frac.get_pi())
            rhs = abs(signed_rhs)
            if not rhs.is_normal():
                continue
            signed_lhs = (u/target_value + target_value/l + c) / d
            lhs = abs(signed_lhs)
            int_rhs = int(rhs)
            int_lhs = int(lhs)
            if self.trunc_integer:
                rhs -= int_rhs
                lhs -= int_lhs
            if self.compare_dec_with_accuracy(rhs, lhs, accuracy):
                if print_clicks:
                    print(ab)
                    print(ulcd)
                    print(rhs)
                    print(lhs)
                    print('')
                u,l,c,d = ulcd
                if (signed_lhs > 0 and signed_rhs < 0) or (signed_lhs < 0 and signed_rhs > 0):
                    d *= -1
                    signed_lhs *= -1
                c += abs(d) * (int_rhs - int_lhs)
                ulcd = (u,l,c,d)
                refined_params.append((ab, ulcd, post_func_ind, convergence_info))

        self.filtered_params = refined_params

    def get_uniq_filtered_params(self):
        return self.uniq_params

    def filter_uniq_params(self):
        non_equiv_params = []
        for params in self.filtered_params:
            is_unique = True
            for uniq_params in non_equiv_params:
                if self._is_equiv_params(params, uniq_params):
                    is_unique = False
                    break
            if is_unique:
                non_equiv_params.append(params)
        self.uniq_params = non_equiv_params
        self.filtered_params = self.uniq_params

    def _is_equiv_params(self, params1, params2):
        ab1, ulcd1, post_func_ind1, convergence_info1 = params1
        ab2, ulcd2, post_func_ind2, convergence_info2 = params2
        if post_func_ind1 != 0 or post_func_ind2 != 0:
            return (ab1 == ab2 and ulcd1 == ulcd2 and post_func_ind1 == post_func_ind2)

        pa1, pb1 = ab1
        pa2, pb2 = ab2
        if len(pa1) != len(pa2) or len(pb1) != len(pb2):
            return False
        if len(pa1) > 1:
            return all([ self._is_equiv_params((((p1,), pb1), ulcd1, post_func_ind1),
                                               (((p2,), pb2), ulcd2, post_func_ind2)) for p1, p2 in zip(pa1, pa2)])
        if len(pb1) > 1:
            return all([ self._is_equiv_params(((pa1, (p1,)), ulcd1, post_func_ind1),
                                               ((pa2, (p2,)), ulcd2, post_func_ind2)) for p1, p2 in zip(pb1, pb2)])

        # if were're here, then params1 and params2 are single-element tuples
        pa1, pa2, pb1, pb2 = pa1[0], pa2[0], pb1[0], pb2[0]
        if pb1 == pb2:
            u1, l1, c1, d1 = ulcd1
            u2, l2, c2, d2 = ulcd2
            for s in [1, -1]:
                if u1 == u2 * s and l1 == l2 * s and c1 == c2 * s and d1 == -d2 * s:
                    return True
            sys.stdout.flush()
            if pa1 == pa2 or list(pa1) == [ -x for x in pa2 ]:
                ulcd_ratio = abs(d1 / d2)
                if u1 == u2 * ulcd_ratio and l1 == l2 / ulcd_ratio and c1 == c2 * ulcd_ratio:
                    return True
                if u1 == -u2 * ulcd_ratio and l1 == -l2 / ulcd_ratio and c1 == -c2 * ulcd_ratio:
                    return True
            if ulcd1 == ulcd2 and (pa1 == pa2 or list(pa1) == [ -x for x in pa2 ]):
                return True

        return False
        # TODO: Finish this

    def filter_only_exp_convergence(self):  # , filter_uniq_list=True):
        # if filter_uniq_list:
        #     params_list = self.uniq_params
        # else:
        params_list = self.filtered_params

        filtered_params_list = []
        for cf_params in params_list:
            ab, ulcd, post_func_ind, convergence_info = cf_params
            pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            if pi_cont_frac.is_convergence_exponential():
                filtered_params_list.append(cf_params)

        # if filter_uniq_list:
        #     self.uniq_params = params_list
        # else:
        self.filtered_params = filtered_params_list

    def filter_clicks_by_approach_type(self, whitelist=['exp', 'over_exp', 'fast'], blacklist=None):    # , filter_uniq_list=True):
        if not any([whitelist, blacklist]):
            raise ValueError('One is required: whitelist, blacklist')
        if whitelist and blacklist:
            raise ValueError('Only one is possible: whitelist, blacklist')
        # if filter_uniq_list:
        #     params_list = self.uniq_params
        # else:
        params_list = self.filtered_params

        filtered_params_list = []
        for cf_params in params_list:
            ab, ulcd, post_func_ind, convergence_info = cf_params
            pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            pi_cont_frac.estimate_approach_type_and_params()
            approach_type, approach_params = pi_cont_frac.get_approach_type_and_params()
            if (whitelist and approach_type in whitelist) or (blacklist and approach_type not in blacklist):
                cf_params = (ab, ulcd, post_func_ind, (approach_type, approach_params))
                filtered_params_list.append(cf_params)

        # if filter_uniq_list:
        #     self.uniq_params = params_list
        # else:
        self.filtered_params = filtered_params_list

    def fix_ulcd_constant(self, params):
        ab, ulcd, post_func_ind = params
        pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
        # pi_cont_frac.reinitialize(a_coeffs=ab[0], b_coeffs=ab[1])
        pi_cont_frac.gen_iterations(20)

    def get_filtered_params(self):
        return self.filtered_params

    def export_to_csv(self, filename, postfuncs, uniq_params=False):
        csvfile = open(filename, 'w', newline='')
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['postproc_funcs', postfuncs])
        csvwriter.writerow(['a poly [a_0, a_1, ...]', 'b poly  [b_0, b_1, ...]', 'procpost_func index',
                            'u', 'l', 'c', 'd',
                            'convergence type', 'convergence rate', 'postfunc(cont_frac)', '(u/pi+pi/l+c)/d'])
        for ab,ulcd, post_func_ind, convergence_info in [self.filtered_params, self.uniq_params][uniq_params]:
            pa, pb = ab
            u, l, c, d = ulcd
            pi_cont_frac, postproc_res, modified_pi_expression =self.build_pi_from_params((ab, ulcd, post_func_ind,
                                                                                           convergence_info))
            csvwriter.writerow([pa, pb, post_func_ind, u, l, c, d, convergence_info[0], convergence_info[1],
                                postproc_res.to_eng_string(), modified_pi_expression.to_eng_string()])
        csvfile.close()

    def build_pi_from_params(self, params, iterations=3000):
        target_value = self.target_generator()
        ab, ulcd, post_func_ind, convergence_info = params
        u,l,c,d = ulcd
        pa, pb = ab
        pi_cont_frac = cont_fracs.PiContFrac(a_coeffs=pa, b_coeffs=pb)
        pi_cont_frac.gen_iterations(iterations)
        modified_pi_expression = (u/target_value + target_value/l + c) / d
        return (pi_cont_frac, self.postproc_funcs[post_func_ind](pi_cont_frac.get_pi()), modified_pi_expression)

    @staticmethod
    def compare_dec_with_accuracy(d1, d2, accuracy):
        # +1 for decimal dot
        accuracy += 1
        return d1.to_eng_string()[:accuracy+1] == d2.to_eng_string()[:accuracy+1]

