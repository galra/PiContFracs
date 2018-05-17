import cont_fracs
from decimal_hashtable import DecimalHashTable
from decimal import Decimal as dec
from gen_real_consts import gen_real_pi
import csv
import sys
import time
from latex import latex_cont_frac
from lhs_evaluators import ULCDEnumerator, ULCDEvaluator, RationalFuncEnumerator, RationalFuncEvaluator
from postprocfuncs import POSTPROC_FUNCS_LATEX
from basic_enum_params import BasicEnumPolyParams

ENUMERATOR_TYPES = {'ulcd': ULCDEnumerator,
                    'rationalfunc': RationalFuncEnumerator}
EVALUATOR_TYPES = {'ulcd': ULCDEvaluator,
                   'rationalfunc': RationalFuncEvaluator}

EVALUATOR2TYPE = { v: k for k, v in EVALUATOR_TYPES.items() }
ENUMERATOR2TYPE = { v: k for k, v in ENUMERATOR_TYPES.items() }

# do we still want to keep the default values here? (some of them require special imports that are not needed otherwise)
class MITM:
    def __init__(self, target_generator=gen_real_pi, target_name='pi', postproc_funcs=[lambda x:x], trunc_integer=True,
                 hashtable_prec=6, a_poly_size=3, b_poly_size=3, num_of_a_polys=1, num_of_b_polys=1,
                 enum_only_exp_conv=True, num_of_iterations=100, threshold=None, prec=50):
        self.bep = BasicEnumPolyParams(a_poly_size=a_poly_size, b_poly_size=b_poly_size, num_of_a_polys=num_of_a_polys,
                                       num_of_b_polys=num_of_b_polys, enum_only_exp_conv=enum_only_exp_conv,
                                       avoid_int_roots=True, should_gen_contfrac=True, avoid_zero_b=True,
                                       num_of_iterations=num_of_iterations, threshold=threshold, prec=prec)
        self.target_generator = target_generator
        self.target_name = target_name
        self.postproc_funcs = postproc_funcs
        self.trunc_integer = trunc_integer
        self.hashtable_prec = hashtable_prec
        self.prec = prec
        self.dec_hashtable = DecimalHashTable(self.hashtable_prec)
        self.filtered_params = []

    def redefine_settings(self, target_generator=gen_real_pi, target_name='pi', postproc_funcs=[lambda x:x],
                          trunc_integer=True, hashtable_prec = 6, prec=50):
        self.target_generator = target_generator
        self.target_name = target_name
        self.postproc_funcs = postproc_funcs
        self.trunc_integer = trunc_integer
        self.hashtable_prec = hashtable_prec
        self.dec_hashtable.update_accuracy(self.hashtable_prec)
        self.filtered_params = []
        set_precision(prec)

    def build_hashtable(self, enum_range=None, range_a=None, range_b=None, prec=None):
        pg, pg_len = self.bep.polys_generator(enum_range=enum_range, range_a=range_a, range_b=range_b, prec=prec)
        self._iter2hashtalbe(pg, pg_len)

    def _iter2hashtalbe(self, itr, itr_len, print_problematic=False):
        progress_percentage=0
        # print(itr_len)
        print('0%', end='')
        sys.stdout.flush()

        for i, (cont_frac, pa, pb) in enumerate(itr):
            if int(100 * i / itr_len) > progress_percentage:
                progress_percentage = int(100 * i / itr_len)
                print('\r%d%%' % progress_percentage, end='')
                sys.stdout.flush()
            cur_cont_frac = cont_frac.get_result()
            if not cur_cont_frac.is_normal():
                continue
            # #I got trouble when having super large numbers coming going into the sin calc
            # if print_problematic and abs(cur_cont_frac > dec('1E+50')):
            #     print('problematic number')
            #     print(cur_cont_frac)
            #     continue
            # # I got trouble when having super small numbers because we use inverse and then they get large and go into the sin calc
            # if print_problematic and abs(cur_cont_frac) < dec('1E-50'):
            #     print('problematic number')
            #     print(cur_cont_frac)
            #     continue

            # this tries to replace the previous 'problematic number' check, it haven't been tested yet
            if not cur_cont_frac.is_normal():
                if print_problematic:
                    print('problematic number')
                    print(cur_cont_frac)
                continue

            for post_func_ind, post_f in enumerate(self.postproc_funcs):
                # print(cur_cont_frac)
                try:
                    k = post_f(cur_cont_frac)
                except:
                    print(k)
                    raise
                if not k.is_normal():
                    continue
                k = abs(k)
                if self.trunc_integer:
                    k -= int(k)
                if k not in self.dec_hashtable:
                    self.dec_hashtable[k] = []
                self.dec_hashtable[k].append(((pa, pb), post_func_ind))
        print()

    def find_clicks(self, lhs_type, lhs_enumerator_params):
        # if isinstance(u_range, int):
        #     u_range = range(-u_range, u_range+1)
        # if isinstance(l_range, int):
        #     l_range = range(-l_range, l_range+1)
        # if isinstance(c_range, int):
        #     c_range = range(-c_range, c_range+1)
        # if isinstance(d_range, int):
        #     d_range = range(1, d_range+1)
        filtered_params = []
        progress_percentage=0
        print('0%', end='')
        lhs_enumerator = ENUMERATOR_TYPES[lhs_type](lhs_enumerator_params, self.target_generator())
        lhs_evaluator = EVALUATOR_TYPES[lhs_type]
        lhs_generator = lhs_enumerator.generator()
        lhs_iter_len = len(lhs_enumerator)
        for i, enum_res_obj in enumerate(lhs_generator):
            r = abs(enum_res_obj.get_val())
            if self.trunc_integer:
                r -= int(r)
            if r in self.dec_hashtable or -r in self.dec_hashtable:
                filtered_params.extend([ (ab, lhs_evaluator(enum_res_obj), post_func_ind, None)
                                         for ab, post_func_ind in self.dec_hashtable[r] ])
            if int(100 * i / lhs_iter_len) > progress_percentage:
                progress_percentage = int(100 * i / lhs_iter_len)
                print('\r%d%%' % progress_percentage, end='')
        print('')
        self.filtered_params = filtered_params

    def refine_clicks(self, accuracy=10, num_of_iterations=3000, print_clicks=False):
        refined_params = []
        target_value = self.target_generator()
        # cont_frac = cont_fracs.ContFrac()
        progress_percentage = 0
        filtered_params_len = len(self.filtered_params)
        print('0%\r', end='')
        for i, (ab, lhs_res_obj, post_func_ind, convergence_info) in enumerate(self.filtered_params):
            cont_frac = cont_fracs.ContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            cont_frac.set_approach_type_and_params(convergence_info)
            # cont_frac.reinitialize(a_coeffs=ab[0], b_coeffs=ab[1])
            cont_frac.gen_iterations(num_of_iterations, dec('1E-%d' % (accuracy+100)))
            signed_rhs = self.postproc_funcs[post_func_ind](cont_frac.get_result())
            rhs = abs(signed_rhs)
            if not rhs.is_normal():
                continue
            signed_lhs = lhs_res_obj.get_val()
            lhs = abs(signed_lhs)
            int_rhs = int(rhs)
            int_lhs = int(lhs)
            if self.trunc_integer:
                rhs -= int_rhs
                lhs -= int_lhs
            if self.compare_dec_with_accuracy(rhs, lhs, accuracy):
                if print_clicks:
                    print(ab)
                    print(lhs_res_obj)
                    print(rhs)
                    print(lhs)
                    print('')
                if (signed_lhs > 0 and signed_rhs < 0) or (signed_lhs < 0 and signed_rhs > 0):
                    lhs_res_obj.flip_sign()
                    signed_lhs *= -1
                lhs_res_obj.add_int(int(signed_rhs - signed_lhs))
                refined_params.append((ab, lhs_res_obj, post_func_ind, convergence_info))
            if int(100 * i / filtered_params_len) > progress_percentage:
                progress_percentage = int(100 * i / filtered_params_len)
                print('\r%d%%' % progress_percentage, end='')
        print('')
        self.filtered_params = refined_params

    def filter_uniq_params(self):
        non_equiv_params = []
        for params in self.filtered_params:
            is_unique = True
            for uniq_params in non_equiv_params:
                _, uniq_lhs_res_obj, _, _ = uniq_params
                if uniq_lhs_res_obj.is_equiv(uniq_params, params):
                    is_unique = False
                    break
            if is_unique:
                non_equiv_params.append(params)
        self.filtered_params = non_equiv_params

    def filter_only_exp_convergence(self, print_surprising_nonexp_contfracs=False):  # , filter_uniq_list=True):
        # if filter_uniq_list:
        #     params_list = self.uniq_params
        # else:
        params_list = self.filtered_params

        progress_percentage = 0
        filtered_params_len = len(params_list)
        print('0%\r', end='')

        filtered_params_list = []
        for i, cf_params in enumerate(params_list):
            ab, lhs_res_obj, post_func_ind, convergence_info = cf_params
            cont_frac = cont_fracs.ContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            if cont_frac.is_convergence_fast():
                filtered_params_list.append(cf_params)
            elif print_surprising_nonexp_contfracs:
                print('Surprising non-exponential convergence continuous fraction:')
                print(cf_params)
                cont_frac.estimate_approach_type_and_params()
                print(cont_frac.get_approach_type_and_params())
            if int(100 * i / filtered_params_len) > progress_percentage:
                progress_percentage = int(100 * i / filtered_params_len)
                print('\r%d%%' % progress_percentage, end='')
        print('')
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

        progress_percentage = 0
        filtered_params_len = len(params_list)
        print('0%\r', end='')

        filtered_params_list = []
        for i, cf_params in enumerate(params_list):
            ab, lhs_res_obj, post_func_ind, convergence_info = cf_params
            cont_frac = cont_fracs.ContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            try:
                cont_frac.estimate_approach_type_and_params()
                approach_type, approach_params = cont_frac.get_approach_type_and_params()
            except Exception as e:
                print('Problems while estimating the following cf_params in "filter_clicks_by_approach_type"')
                print('Exception has occurred. Skipping.')
                print(cf_params)
                print(e)
                continue
            if (whitelist and approach_type in whitelist) or (blacklist and approach_type not in blacklist):
                cf_params = (ab, lhs_res_obj, post_func_ind, (approach_type, approach_params))
                filtered_params_list.append(cf_params)
            if int(100 * i / filtered_params_len) > progress_percentage:
                progress_percentage = int(100 * i / filtered_params_len)
                print('\r%d%%' % progress_percentage, end='')
        print('')
        # if filter_uniq_list:
        #     self.uniq_params = params_list
        # else:
        self.filtered_params = filtered_params_list

    def get_filtered_params(self):
        return self.filtered_params

    def delete_hashtable(self):
        del self.dec_hashtable
        self.dec_hashtable = DecimalHashTable(self.hashtable_prec)

    def get_results_as_eqns(self, postfuncs):
        eqns = []
        eval_poly = cont_fracs.ContFrac._array_to_polynom
        known_targets = {'pi': '\pi',
                        'phi': r'\varphi'}

        for ab, lhs_res_obj, post_func_ind, convergence_info in self.filtered_params:
            pa, pb = ab
            cont_frac, postproc_res, lhs_res_obj = self.build_contfrac_from_params((ab, lhs_res_obj, post_func_ind,
                                                                                                   convergence_info))

            depth = 5
            a = [eval_poly(pa[0], i) for i in range(depth)]
            b = [eval_poly(pb[0], i) for i in range(depth)]

            if self.target_name in known_targets:
                target_name = known_targets[self.target_name]
            else:
                target_name = self.target_name

            # Creates the equation object
            lhs = lhs_res_obj.get_latex_exp(target_name)
            rhs = latex_cont_frac(a, b)
            lhs, rhs = POSTPROC_FUNCS_LATEX[postfuncs[post_func_ind]](lhs, rhs)
            eqn = r'{0} = {1}'.format(lhs, rhs)

            # Appends equation
            eqns.append(eqn)

        return eqns

    def export_to_csv(self, filename, postfuncs):
        csvfile = open(filename, 'w', newline='')
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['postproc_funcs', postfuncs, 'target_name', self.target_name])
        csvwriter.writerow(['a poly [a_0, a_1, ...]', 'b poly  [b_0, b_1, ...]', 'procpost_func index',
                            'LHS type', 'LHS params',
                            'convergence type', 'convergence rate', 'postfunc(cont_frac)',
                            'LHS val'])
        for ab, lhs_res_obj, post_func_ind, convergence_info in self.filtered_params:
            pa, pb = ab
            lhs_type = EVALUATOR2TYPE[type(lhs_res_obj)]
            cont_frac, postproc_res, lhs_res_obj = self.build_contfrac_from_params((ab, lhs_res_obj,
                                                                                                   post_func_ind,
                                                                                                   convergence_info))
            csvwriter.writerow([pa, pb, post_func_ind, lhs_type, lhs_res_obj.get_params(),
                                convergence_info[0], convergence_info[1],
                                postproc_res.to_eng_string(), lhs_res_obj.get_val().to_eng_string()])
        csvfile.close()

    def build_contfrac_from_params(self, params, iterations=3000):
        ab, lhs_res_obj, post_func_ind, convergence_info = params
        pa, pb = ab
        cont_frac = cont_fracs.ContFrac(a_coeffs=pa, b_coeffs=pb)
        cont_frac.gen_iterations(iterations)
        return (cont_frac, self.postproc_funcs[post_func_ind](cont_frac.get_result()), lhs_res_obj)

    @staticmethod
    def compare_dec_with_accuracy(d1, d2, accuracy):
        # +1 for decimal dot
        accuracy += 1
        return d1.to_eng_string()[:accuracy+1] == d2.to_eng_string()[:accuracy+1]

