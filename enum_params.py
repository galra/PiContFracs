import cont_fracs
from decimal_hashtable import DecimalHashTable
from decimal import Decimal as dec
from gen_consts import gen_pi_const
import csv
import sys
from latex import latex_cont_frac
from lhs_evaluators import ULCDEnumerator, ULCDEvaluator, RationalFuncEnumerator, RationalFuncEvaluator
from postprocfuncs import POSTPROC_FUNCS ,POSTPROC_FUNCS_LATEX
from basic_enum_params import BasicEnumPolyParams, set_precision
import io
import shutil
import numpy

ENUMERATOR_TYPES = {'ulcd': ULCDEnumerator,
                    'rationalfunc': RationalFuncEnumerator}
EVALUATOR_TYPES = {'ulcd': ULCDEvaluator,
                   'rationalfunc': RationalFuncEvaluator}

EVALUATOR2TYPE = { v: k for k, v in EVALUATOR_TYPES.items() }
ENUMERATOR2TYPE = { v: k for k, v in ENUMERATOR_TYPES.items() }

# do we still want to keep the default values here? (some of them require special imports that are not needed otherwise)
class MITM:
    def __init__(self, target_generator=gen_pi_const, target_name='pi', postproc_funcs=[lambda x: x],
                 postproc_funcs_filter=[], trunc_integer=True, hashtable_prec=15, ab_poly_class=BasicEnumPolyParams,
                 enum_only_exp_conv=True, num_of_iterations=300, threshold=None, prec=50):
        # rhs_polys_enumer = Basic Enum Params
        self.rhs_polys_enumer = ab_poly_class(num_of_iterations=num_of_iterations,
                                              enum_only_exp_conv=enum_only_exp_conv, avoid_int_roots=True,
                                              should_gen_contfrac=True, avoid_zero_b=True,
                                              threshold=threshold, prec=prec)
        set_precision(prec)
        self.target_generator = target_generator
        self.target_name = target_name
        self.postproc_funcs = postproc_funcs
        self.postproc_funcs_filter = postproc_funcs_filter
        self.trunc_integer = trunc_integer
        self.hashtable_prec = hashtable_prec
        self.prec = prec
        self.dec_hashtable = DecimalHashTable(self.hashtable_prec)
        self.filtered_params = []


    def redefine_settings(self, target_generator=gen_pi_const, target_name='pi', postproc_funcs=[lambda x:x],
                          postproc_funcs_filter=[], trunc_integer=True, ab_poly_class=BasicEnumPolyParams,
                          hashtable_prec = 6, prec=50):
        if not isinstance(self.rhs_polys_enumer, ab_poly_class):
            self.rhs_polys_enumer = ab_poly_class(prec=self.prec, enum_only_exp_conv=True, avoid_int_roots=True,
                                                  should_gen_contfrac=True, num_of_iterations=100, threshold=None)
        set_precision(prec)
        self.target_generator = target_generator
        self.target_name = target_name
        self.postproc_funcs = postproc_funcs
        self.postproc_funcs_filter = postproc_funcs_filter
        self.trunc_integer = trunc_integer
        self.hashtable_prec = hashtable_prec
        self.dec_hashtable.update_accuracy(self.hashtable_prec)
        self.filtered_params = []


    def build_hashtable(self, range_a=None, range_b=None, prec=None):
        pg, pg_len = self.rhs_polys_enumer.polys_generator(range_a=range_a, range_b=range_b, prec=prec)
        self._iter2hashtalbe(pg, pg_len)

    def _iter2hashtalbe(self, itr, itr_len, print_problematic=False):
        progress_percentage=0
        # print(itr_len)
        print('0%', end='')
        sys.stdout.flush()

        filtered_postproc_funcs = [(i, f) for i,f in enumerate(self.postproc_funcs) if i in self.postproc_funcs_filter ]

        for i, (cont_frac, pa, pb) in enumerate(itr):
            if int(100 * i / itr_len) > progress_percentage:
                progress_percentage = int(100 * i / itr_len)
                print('\r%d%%' % progress_percentage, end='')
                sys.stdout.flush()
            cur_cont_frac = cont_frac.get_result()
            if not cur_cont_frac.is_normal():
                continue

            # this tries to replace the previous 'problematic number' check, it haven't been tested yet
            if not cur_cont_frac.is_normal():
                if print_problematic:
                    print('problematic number')
                    print(cur_cont_frac)
                continue
            for post_func_ind, post_f in filtered_postproc_funcs:
                # print(cur_cont_frac)
                if post_func_ind not in self.postproc_funcs_filter:
                    pass
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
                if ((pa, pb), post_func_ind) not in self.dec_hashtable[k]:
                    self.dec_hashtable[k].append(((pa, pb), post_func_ind))
        print()

    def find_clicks(self, lhs_type, lhs_enumerator_params):
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
                sys.stdout.flush()
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
                    # print('flipping sign')
                    lhs_res_obj.flip_sign()
                    signed_lhs *= -1
                try:
                    lhs_res_obj.add_int(int(signed_rhs - signed_lhs))
                except:
                    print(signed_rhs)
                    print(signed_lhs)
                    print(int(signed_rhs - signed_lhs))
                    exit()
                lhs_res_obj.canonalize_params()
                refined_params.append((ab, lhs_res_obj, post_func_ind, convergence_info))
            if int(100 * i / filtered_params_len) > progress_percentage:
                progress_percentage = int(100 * i / filtered_params_len)
                print('\r%d%%' % progress_percentage, end='')
                sys.stdout.flush()
        print('')
        self.filtered_params = refined_params

    def filter_uniq_params(self):
        print('0%\r', end='')
        progress_percentage = 0
        filtered_params_len = len(self.filtered_params)
        non_equiv_params = []
        for i, params in enumerate(self.filtered_params):
            is_unique = True
            for uniq_params in non_equiv_params:
                uniq_lhs_res_obj= uniq_params[1]
                if uniq_lhs_res_obj.is_equiv(uniq_params, params):
                    is_unique = False
                    break
            if is_unique:
                non_equiv_params.append(params)
            if int(100 * i / filtered_params_len) > progress_percentage:
                progress_percentage = int(100 * i / filtered_params_len)
                print('\r%2d%%, num. of filtered unique params: %d' % (progress_percentage, len(non_equiv_params)),
                      end='')
                sys.stdout.flush()
        print('')
        self.filtered_params = non_equiv_params

    def filter_integer_roots_numerators(self):
        valid_params = []
        for params in self.filtered_params:
            b = params[0][1]
            for b_p in b:
                b_p_roots = numpy.roots(b_p)
                b_p_real_roots = [ r for r in b_p_roots if abs(r.imag) < 0.001 ]
                b_p_int_roots = [ r for r in b_p_real_roots if abs(r - round(r)) < 0.001 ]
            if not b_p_int_roots:
                valid_params.append(params)
        self.filtered_params = valid_params

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
                sys.stdout.flush()
        print('')
        self.filtered_params = filtered_params_list

    def filter_clicks_by_approach_type(self, whitelist=['exp', 'super_exp', 'fast'], blacklist=None):     # , filter_uniq_list=True):
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
                sys.stdout.flush()
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

    def get_results_as_eqns(self, postfuncs, ignore_zerob_exceptions=True):
        eqns = []
        eval_poly = cont_fracs.ContFrac._array_to_polynom
        known_targets = {'pi': '\pi',
                        'phi': r'\varphi'}

        for ab, lhs_res_obj, post_func_ind, convergence_info in self.filtered_params:
            pa, pb = ab
            try:
                cont_frac, postproc_res, lhs_res_obj = self.build_contfrac_from_params((ab, lhs_res_obj, post_func_ind,
                                                                                        convergence_info))
            except cont_fracs.ZeroB:
                if ignore_zerob_exceptions:
                    continue
                else:
                    raise

            depth = 5
            a = [eval_poly(pa[i % len(pa)], i) for i in range(depth)]
            b = [eval_poly(pb[i % len(pb)], i) for i in range(depth)]

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

    def export_to_csv(self, filename):
        with io.StringIO(newline='') as csvbuffer:
            csvwriter = csv.writer(csvbuffer)
            # csvwriter.writerow(['postproc_funcs', postfuncs, 'target_name', self.target_name])
            csvwriter.writerow(['target_name', self.target_name])
            csvwriter.writerow(['a poly [a_0, a_1, ...]', 'b poly  [b_0, b_1, ...]', 'postproc_func',
                                'LHS type', 'LHS params',
                                'convergence type', 'convergence rate', 'postfunc(cont_frac)',
                                'LHS val'])
            for ab, lhs_res_obj, post_func_ind, convergence_info in self.filtered_params:
                pa, pb = ab
                lhs_type = EVALUATOR2TYPE[type(lhs_res_obj)]
                cont_frac, postproc_res, lhs_res_obj = self.build_contfrac_from_params((ab, lhs_res_obj,
                                                                                        post_func_ind,
                                                                                        convergence_info))
                # csvwriter.writerow([pa, pb, post_func_ind, lhs_type, lhs_res_obj.get_params(),
                #                     convergence_info[0], convergence_info[1],
                #                     postproc_res.to_eng_string(), lhs_res_obj.get_val().to_eng_string()])
                csvwriter.writerow([pa, pb, POSTPROC_FUNCS[post_func_ind], lhs_type, lhs_res_obj.get_params(),
                                    convergence_info[0], convergence_info[1],
                                    postproc_res.to_eng_string(), lhs_res_obj.get_val().to_eng_string()])
            with open(filename, 'w', newline='') as csvfile:
                csvbuffer.seek(0)
                shutil.copyfileobj(csvbuffer, csvfile)

    def build_contfrac_from_params(self, params, iterations=3000):
        ab, lhs_res_obj, post_func_ind, convergence_info = params
        pa, pb = ab
        cont_frac = cont_fracs.ContFrac(a_coeffs=pa, b_coeffs=pb)
        cont_frac.gen_iterations(iterations)
        return cont_frac, self.postproc_funcs[post_func_ind](cont_frac.get_result()), lhs_res_obj

    @staticmethod
    def compare_dec_with_accuracy(d1, d2, accuracy):
        # +1 for decimal dot
        accuracy += 1
        return d1.to_eng_string()[:accuracy+1] == d2.to_eng_string()[:accuracy+1]

