import numpy as np
import basic_algo
from decimal import Decimal as dec
import decimal
from gen_real_consts import gen_real_pi
import math
import scipy.stats

# used to avoid numpy casting int into numpy.int32/64, which doesn't handle bignums
# class myint(int):
#     pass


class ZeroB(Exception):
    pass


# The default diff_mat and first_diff_mat compute the following:
#   The first diff matrix is a_0 = a_0 / 1
#   then the i'th iteration (a_0 is the 0'th iteration) is a_0 + b_1 / (a1 + b_2 / (... + b_i / a_i))
class PiContFrac(basic_algo.PiBasicAlgo):
    # if b_n = 0 then we get a rational number, and should skip this results. This is what avoid_zero_b for.
    def __init__(self, a_coeffs, b_coeffs, avoid_zero_b=True, check_b_threshold=False, b_threshold=10**10, diff_mat_gen=None,
                 first_diff_mat=None, target_val=None, logging=False, dtype=int):
        if diff_mat_gen is None:
            diff_mat_gen = self._default_diff_mat_gen
        if target_val is None:
            target_val = gen_real_pi()

        self._dtype_0 = dtype(0)
        self._dtype_1 = dtype(1)

        if not isinstance(a_coeffs[0], list) and not isinstance(a_coeffs[0], tuple):
            a_coeffs = [a_coeffs]
        if not isinstance(b_coeffs[0], list) and not isinstance(b_coeffs[0], tuple):
            b_coeffs = [b_coeffs]

        params = {'a_coeffs': a_coeffs, 'b_coeffs': b_coeffs, 'pi': 1}
        if first_diff_mat is None:
            self._auto_first_diff_mat = True
            first_diff_mat = self._autogen_first_diff_mat(params, dtype)
        else:
            self._auto_first_diff_mat = False

        self._avoid_zero_b = avoid_zero_b
        self._check_b_threshold = check_b_threshold
        self._b_threshold = dtype(b_threshold)
        self._approach_type, self._approach_params = 'undefined', 0
        super().__init__(params, diff_mat_gen, first_diff_mat, target_val, logging, dtype)

    def reinitialize(self, a_coeffs=None, b_coeffs=None, first_diff_mat=None, **kwargs):
        self._approach_type, self._approach_params = 'undefined', 0
        params = kwargs
        if a_coeffs:
            params['a_coeffs'] = a_coeffs
        else:
            params['a_coeffs'] = self._params['a_coeffs']
        if b_coeffs:
            params['b_coeffs'] = b_coeffs
        else:
            params['b_coeffs'] = self._params['b_coeffs']
        if first_diff_mat is None and self._auto_first_diff_mat:
            first_diff_mat = self._autogen_first_diff_mat(params)

        super().reinitialize(params, first_diff_mat)

    def iteration_algorithm(self, params, diff_mat, i):
        # diff_mat is actually diff_vec in this case
        # it is the vector (p_n, p_{n-1}, dp_n/dx_m, dp_{n-1}/dx_m)
        # and similarly for (q_n, q_{n-1}, dq_n/dx_m, dq_{n-1}/dx_m)
        i += 1
        a_coeffs = params['a_coeffs']
        b_coeffs = params['b_coeffs']
        # ts = time.time()
        # for k in range(100):
        diff_mats_p, diff_mats_q = self._diff_mat_gen(i, a_coeffs, b_coeffs)
        # print((time.time()-ts)/100)
        diff_vec_p, diff_vec_q = diff_mat

        # new_vec_p = [ new_mat * diff_vec for new_mat, diff_vec in zip(diff_mats_p, diff_vec_p) ]
        # new_vec_q = [ new_mat * diff_vec for new_mat, diff_vec in zip(diff_mats_q, diff_vec_q) ]
        # TODO: if calculations are slow, once in a while (i % 100 == 0?) check if gcd(p_i,p_{i-1}) != 0
        # TODO: and take it out as a factor (add p_gcd, q_gcd params)
        new_vec_p = (sum([l*m for l,m in zip(diff_vec_p, diff_mats_p)]), diff_vec_p[0])
        new_vec_q = (sum([l*m for l,m in zip(diff_vec_q, diff_mats_q)]), diff_vec_q[0])

        # new_vec_p is a list of vectors, which fit for dp_i/dx_j. The first element in each vector is p_n
        # Because it will be the same p_n for each dp_i/dx_j (only the derivatives change), we take arbitrarily
        # the first one. Same for new_vec_q
        # p_i = new_vec_p[0][0,0]
        # q_i = new_vec_q[0][0,0]
        if self._logging:
            p_i = dec(new_vec_p[0])
            q_i = dec(new_vec_q[0])
            if q_i.is_zero():
                pi_i = dec('NaN')
            else:
                pi_i = p_i / q_i
            params['pi'] = pi_i
        new_diff_mat = (new_vec_p, new_vec_q)

        # OLD COMMENTS - MAY BE IRRELEVANT:
        # the vectors are (p_n, p_{n-1}, dp_n/dx_m, dp_{n-1}/dx_m)
        # so for x_j we take the j'th vector which suits the differentiation by x_j
        # and we take the 3rd element (which is [2]) for the derivative
        return (params, new_diff_mat)

    def finalize_iterations(self, params, diff_mat, n):
        diff_vec_p, diff_vec_q = diff_mat
        p = dec(diff_vec_p[0])
        q = dec(diff_vec_q[0])
        if q.is_zero():
            pi = dec('NaN')
        else:
            pi = p / q
        params['pi'] = pi
        return (params, diff_mat)

    def is_accuracy_met(self, accuracy, params, diff_mat, i):
        # error < gcd^2/q^2 < accuracy
        # q^2 > gcd^2 * ceil(1 / accuracy) > gcd^2 / accuracy
        # # diff_vec_p, diff_vec_q = diff_mat
        # # p = diff_vec_p[0]
        # # q = diff_vec_q[0]
        # # pq_gcd = math.gcd(p, q)
        #
        # return q**2 > pq_gcd**2 * math.ceil(1/accuracy)
        if i < 1:
            return

    def _default_diff_mat_gen(self, i, a_coeffs, b_coeffs):
        a2p = PiContFrac._array_to_polynom
        # This supports only 1 or 2 polynomials
        a_i = a2p(a_coeffs[i % len(a_coeffs)], (i + len(a_coeffs)-1) >> (len(a_coeffs)-1))
        b_i = a2p(b_coeffs[(i-1) % len(b_coeffs)], (i + len(b_coeffs)-1) >> (len(b_coeffs)-1))
        if self._avoid_zero_b and b_i == 0:
            raise ZeroB()
        if self._check_b_threshold and b_i / self.gcd(a_i, b_i) > self._b_threshold:
            # raise self.SufficientAccuracy()
            pass
        # mat_p = [np.matrix(( (a_i, b_i),
        #                         (1, 0)), self._dtype)]
        # mat_q = [np.matrix(((a_i, b_i),
        #                         (1, 0),), self._dtype)]
        mat_p = (a_i, b_i)
        mat_q = mat_p
        return (mat_p, mat_q)

    def _autogen_first_diff_mat(self, params, dtype=None):
        if dtype is None:
            dtype = self._dtype
        if hasattr(self, '_dtype') and dtype == self._dtype:
            dtype_0 = self._dtype_0
            dtype_1 = self._dtype_1
        else:
            dtype_0 = dtype(0)
            dtype_1 = dtype(1)
        a_coeffs = params['a_coeffs']
        # b_coeffs = params['b_coeffs']
        # first_diff_mat = ([], [])
        # first_diff_mat[0].append(np.matrix(((dtype(a_coeffs[0]),), (dtype_1,),), dtype=dtype))
        # first_diff_mat[1].append(np.matrix(((dtype_1,), (dtype_0,)), dtype=dtype))
        # return first_diff_mat
        return (dtype(a_coeffs[0][0]), dtype_1), (dtype_1, dtype_0)

    def compare_result(self, target_val=None):
        if target_val is not None:
            comp_val = target_val
        elif self._target_val is not None:
            comp_val = self._target_val
        else:
            raise ValueError('No target value to compare to.')

        try:
            r = self.params_log['pi'][-1] - comp_val
            # if not r.is_normal():
            return abs(r)
            # return abs(r - round(r)
        except:
            raise RuntimeError('Run gen_iterations first. No PI was generated!')

    def get_pi(self):
        return self.params_log['pi'][-1]

    def get_p_q(self):
        diff_vec_p, diff_vec_q = self.diff_mats[-1]
        return diff_vec_p[0], diff_vec_q[0]

    def is_pi_valid(self):
        return self.params_log['pi'][-1].is_normal()

    def estimate_approach_type_and_params(self):
        iters = 5000
        initial_cutoff = 1500
        iters_step = 500
        approach_type, approach_params = self._estimate_approach_type_and_params_inner_alg(find_poly_parameter=True,
                                                                                           iters=iters,
                                                                                           initial_cutoff=initial_cutoff,
                                                                                           iters_step=iters_step)
        while approach_type == 'fast' and initial_cutoff > 0:
            initial_cutoff >>= 1
            iters = max(50, iters >> 1)
            iters_step = max(int((iters-initial_cutoff) / 30), iters_step >> 1)
            approach_type, approach_params = self._estimate_approach_type_and_params_inner_alg(find_poly_parameter=True,
                                                                                               iters=iters,
                                                                                               initial_cutoff=initial_cutoff,
                                                                                               iters_step=iters_step)
        if initial_cutoff == 0:
            iters = 50
            iters_step = 1
            while approach_type == 'fast' and iters > 10:
                initial_cutoff >>= 1
                iters -= 5
                approach_type, approach_params = self._estimate_approach_type_and_params_inner_alg(find_poly_parameter=True,
                                                                                                   iters=iters,
                                                                                                   initial_cutoff=initial_cutoff,
                                                                                                   iters_step=iters_step)
        self._approach_type = approach_type
        self._approach_params = approach_params
        return

    def get_approach_type_and_params(self):
        return (self._approach_type, self._approach_params)

    def set_approach_type_and_params(self, convergence_info):
        self._approach_type, self._approach_params = convergence_info

    def is_convergence_exponential(self, find_poly_parameter=False, iters=5000, initial_cutoff=1500,
                                                     iters_step=500, exponential_threshold=1.1):
        """Returns true if the convergence type is exponential or over exponential.
False if it's sub exponential (e.g. linear)."""
        if iters_step < 6:
            ValueError('iters_step should be at least 4')

        self.gen_iterations(initial_cutoff, exec_finalize=False)

        return_val = True
        for i in range(initial_cutoff+iters_step, iters+1, iters_step):
            p_0, q_0 = self.get_p_q()
            self.add_iterations(1, exec_finalize=False)
            p_1, q_1 = self.get_p_q()
            self.add_iterations(1, exec_finalize=False)
            p_2, q_2 = self.get_p_q()
            self.add_iterations(1, exec_finalize=False)
            p_3, q_3 = self.get_p_q()
            self.add_iterations(1, exec_finalize=False)
            p_4, q_4 = self.get_p_q()
            self.add_iterations(1, exec_finalize=False)
            p_5, q_5 = self.get_p_q()

            # q_4(p_2q_0-p_0q_2) > q_0(p_4q_2-p_2q_4)
            lhs_pair = abs(q_4 * (p_2 * q_0 - p_0 * q_2))
            rhs_pair = abs(q_0 * (p_4 * q_2 - p_2 * q_4))
            lhs_odd = abs(q_5 * (p_3 * q_1 - p_1 * q_3))
            rhs_odd = abs(q_1 * (p_5 * q_3 - p_3 * q_5))
            if (((lhs_pair - rhs_pair) <= (rhs_pair >> 4) ) or
                    ((lhs_odd - rhs_odd) <= (rhs_odd >> 4))):
                return_val = False
                break
            # -3 for the iterations of res_1, res_2, res_3 that were already executed
            self.add_iterations(iters_step - 5, exec_finalize=False)

        self.gen_iterations(0)
        return return_val

    def _estimate_approach_type_and_params_inner_alg(self, find_poly_parameter=False, iters=5000, initial_cutoff=1500,
                                          iters_step=500, exponential_threshold=1.1):
        """Returns 'exp', 'over_exp', 'poly', 'undefined', 'fast' and 'mixed', as a tuple of (string,num): (approach_type, approach_parameter)
    or ('poly', (approach_parameter, R**2))."""
        if iters_step < 6:
            ValueError('iters_step should be at least 4')

        approach_type = None
        approach_parameter = 0

        delta_pair = []
        delta_odd = []
        self.gen_iterations(initial_cutoff)
        res_0 = self.get_pi()
        # if res_0.is_nan():
        #     print('res_0 nan')
        #     print(self.get_p_q())
        self.add_iterations(1)
        res_1 = self.get_pi()
        # if res_1.is_nan():
        #     print('res_1 nan')
        #     print(self.get_p_q())
        self.add_iterations(1)
        res_2 = self.get_pi()
        # if res_2.is_nan():
        #     print('res_2 nan')
        #     print(self.get_p_q())
        self.add_iterations(1)
        res_3 = self.get_pi()
        # if res_3.is_nan():
        #     print('res_3 nan')
        #     print(self.get_p_q())
        self.add_iterations(1)
        res_4 = self.get_pi()
        # if res_4.is_nan():
        #     print('res_4 nan')
        #     print(self.get_p_q())
        self.add_iterations(1)
        res_5 = self.get_pi()
        # if res_5.is_nan():
        #     print('res_5 nan')
        #     print(self.get_p_q())
        delta_pair.append((initial_cutoff, abs(res_2 - res_0)))
        # if not delta_pair[-1][1].is_normal() and not delta_pair[-1][1].is_zero():
        #     print('first, res_0: %s' % res_0.to_eng_string())
        #     print('first, res_2: %s' % res_2.to_eng_string())
        delta_pair.append((initial_cutoff + 2, abs(res_4 - res_2)))
        # if not delta_pair[-1][1].is_normal() and not delta_pair[-1][1].is_zero():
        #     print('first, res_2: %s' % res_2.to_eng_string())
        #     print('first, res_4: %s' % res_4.to_eng_string())
        delta_odd.append((initial_cutoff + 1, abs(res_3 - res_1)))
        # if not delta_odd[-1][1].is_normal() and not delta_odd[-1][1].is_zero():
        #     print('first, res_1: %s' % res_1.to_eng_string())
        #     print('first, res_3: %s' % res_3.to_eng_string())
        delta_odd.append((initial_cutoff + 3, abs(res_5 - res_3)))
        # if not delta_odd[-1][1].is_normal() and not delta_odd[-1][1].is_zero():
        #     print('first, res_3: %s' % res_3.to_eng_string())
        #     print('first, res_5: %s' % res_5.to_eng_string())

        for i in range(initial_cutoff+iters_step, iters+1, iters_step):
            # -3 for the iterations of res_1, res_2, res_3 that were already executed
            self.add_iterations(iters_step - 5)
            res_0 = self.get_pi()
            self.add_iterations(1)
            res_1 = self.get_pi()
            self.add_iterations(1)
            res_2 = self.get_pi()
            self.add_iterations(1)
            res_3 = self.get_pi()
            self.add_iterations(1)
            res_4 = self.get_pi()
            self.add_iterations(1)
            res_5 = self.get_pi()
            delta_pair.append((i, abs(res_2 - res_0)))
            # if not delta_pair[-1][1].is_normal() and not delta_pair[-1][1].is_zero():
            #     print('res_0: %s' % res_0.to_eng_string())
            #     print('res_2: %s' % res_2.to_eng_string())
            delta_pair.append((i + 2, abs(res_4 - res_2)))
            # if not delta_pair[-1][1].is_normal() and not delta_pair[-1][1].is_zero():
            #     print('res_2: %s' % res_2.to_eng_string())
            #     print('res_4: %s' % res_4.to_eng_string())
            delta_odd.append((i + 1, abs(res_3 - res_1)))
            # if not delta_odd[-1][1].is_normal() and not delta_odd[-1][1].is_zero():
            #     print('res_1: %s' % res_1.to_eng_string())
            #     print('res_3: %s' % res_3.to_eng_string())
            delta_odd.append((i + 3, abs(res_5 - res_3)))
            # if not delta_odd[-1][1].is_normal() and not delta_odd[-1][1].is_zero():
            #     print('res_3: %s' % res_3.to_eng_string())
            #     print('res_5: %s' % res_5.to_eng_string())
            # if show_progress and i % 500 == 0:
            #     print('\r%d' % i, end='')
        # return (delta_pair, delta_odd)
        pair_diminish = False
        odd_diminish = False
        if len(delta_pair) > 3 and all([ p[1].is_zero() for p in delta_pair[-3:] ]):
            pair_diminish = True
        if len(delta_odd) > 3 and all([ p[1].is_zero() for p in delta_odd[-3:] ]):
            odd_diminish = True
        # if one diminishes and the other isn't, return 'undefined'
        if pair_diminish ^ odd_diminish:
            approach_type = 'undefined'
        elif pair_diminish and odd_diminish:
            approach_type = 'fast'

        # if approach_type:
        #     return (approach_type, approach_parameter)

        pair_ratio = [ (delta_pair[i][0], delta_pair[i][1] / delta_pair[i+1][1])
                       for i in range(0, len(delta_pair), 2) if delta_pair[i][1] != 0 and delta_pair[i+1][1] != 0 and
                      not delta_pair[i][1].is_nan() and not delta_pair[i+1][1].is_nan() ]
        odd_ratio = [ (delta_odd[i][0], delta_odd[i][1] / delta_odd[i+1][1])
                      for i in range(0, len(delta_odd), 2) if delta_odd[i][1] != 0 and delta_odd[i+1][1] != 0 and
                      not delta_odd[i][1].is_nan() and not delta_odd[i+1][1].is_nan() ]

        if len(pair_ratio) < 6:
            return (approach_type, approach_parameter)

        mean_pair_ratio = sum([ p for i, p in pair_ratio] ) / len(pair_ratio)
        mean_pair_ratio_avg_square_error = sum([ (r-mean_pair_ratio)**2 for i, r in pair_ratio ]).sqrt() / len(pair_ratio)
        mean_odd_ratio = sum([ p for i, p in odd_ratio ]) / len(odd_ratio)
        mean_odd_ratio_avg_square_error = sum([ (r-mean_odd_ratio)**2 for i, r in odd_ratio ]).sqrt() / len(odd_ratio)
        relative_pair_sq_err = mean_pair_ratio_avg_square_error / mean_pair_ratio
        relative_odd_sq_err = mean_odd_ratio_avg_square_error / mean_odd_ratio
        if relative_pair_sq_err > 0.5 or relative_odd_sq_err > 0.5:
            if all([ i[1] > 2 for i in (pair_ratio[3*int(len(pair_ratio)/4):] +
                                         odd_ratio[3*int(len(odd_ratio)/4):]) ]):

                approach_type = 'over_exp'
            else:
                approach_type = 'undefined'
        if (relative_odd_sq_err <= 0.5 and relative_pair_sq_err <= 0.5) or approach_type == 'over_exp':
            is_pair_exp = mean_pair_ratio > exponential_threshold
            is_odd_exp = mean_odd_ratio > exponential_threshold
            # in case one is exponential and the other isn't return 'mixed'
            if is_pair_exp ^ is_odd_exp:
                approach_type = 'mixed'
            elif is_pair_exp and is_odd_exp:
                if approach_type != 'over_exp':
                    approach_type = 'exp'
                approach_parameter_pair = mean_pair_ratio**type(mean_pair_ratio)(0.5)
                approach_parameter_odd = mean_odd_ratio**type(mean_odd_ratio)(0.5)
                approach_parameter = min(approach_parameter_pair, approach_parameter_odd)
                approach_coeff_pair = [ abs(delta_pair[i][1] * approach_parameter**(delta_pair[i][0]) /
                                        (1 - approach_parameter**(-2))) for i in range(0, len(delta_pair))
                                         if delta_pair[i][1] != 0 ]
                approach_coeff_pair = sum(approach_coeff_pair) / len(approach_coeff_pair)
                approach_coeff_odd = [ abs(delta_odd[i][1] * approach_parameter**(delta_odd[i][0]) /
                                       (1 - approach_parameter**(-2))) for i in range(0, len(delta_odd))
                                       if delta_odd[i][1] != 0 ]
                approach_coeff_odd = sum(approach_coeff_odd) / len(approach_coeff_odd)
                approach_coeff = min(approach_coeff_pair, approach_coeff_odd)
                approach_parameter = (approach_parameter, approach_coeff)
            else:
                approach_type = 'poly'

        if approach_type != 'poly' or not find_poly_parameter:
            return (approach_type, approach_parameter)

    #     We're requested to find the poly parameter
        log_x_pair = [ math.log(i) for i, d in delta_pair ]
        log_y_pair = [ math.log(d) for i, d in delta_pair ]
        slope_pair, intercept_pair, r_value_pair, p_value_pair, std_err_pair = scipy.stats.linregress(log_x_pair,
                                                                                                      log_y_pair)
        log_x_odd = [ math.log(i) for i, d in delta_odd ]
        log_y_odd = [ math.log(d) for i, d in delta_odd ]
        slope_odd, intercept_odd, r_value_odd, p_value_odd, std_err_odd = scipy.stats.linregress(log_x_odd, log_y_odd)

        approach_parameter = (min(abs(slope_pair), abs(slope_odd))-1, min(intercept_pair, intercept_odd),
                              min(r_value_pair**2, r_value_odd**2))
        return (approach_type, approach_parameter)

    @staticmethod
    def _array_to_polynom(coeffs, x):
        return sum([coeffs[j] * x ** j for j in range(len(coeffs))])

    @staticmethod
    def gcd(a, b):
        """Calculate the Greatest Common Divisor of a and b.

        Unless b==0, the result will have the same sign as b (so that when
        b is divided by it, the result comes out positive).
        """
        while b:
            a, b = b, a%b
        return a