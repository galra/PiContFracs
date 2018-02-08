import numpy as np
import basic_algo
from decimal import Decimal as dec
import decimal
from gen_real_pi import gen_real_pi


class PiContFrac(basic_algo.PiBasicAlgo):
    def __init__(self, a_coeffs, b_coeffs, diff_mat_gen=None, first_diff_mat=None, target_val=None, dtype=dec):
        if first_diff_mat is None:
            first_diff_mat = ([], [])
            first_diff_mat[0].extend([ np.matrix(((dec(a_coeffs[0]),), (dec(1),), (dec(a_coeffs[i] * (i==0)),), (dec(0),)), dtype=dtype)
                                for i in range(len(a_coeffs)) ])
            first_diff_mat[0].extend([ np.matrix(((dec(a_coeffs[0]),), (dec(1),), (dec(0),), (dec(0),)), dtype=dtype)
                                for i in range(len(a_coeffs)) ])
            first_diff_mat[1].extend([ np.matrix(((dec(1),), (dec(0),), (dec(0),), (dec(0),)), dtype=dtype)
                                 for i in range(len(b_coeffs) + len(a_coeffs)) ])
        if diff_mat_gen is None:
            diff_mat_gen = self._default_diff_mat_gen
        if target_val is None:
            target_val = gen_real_pi()

        params = {'a_coeffs': a_coeffs, 'b_coeffs': b_coeffs, 'pi': 1}
        super().__init__(params, diff_mat_gen, first_diff_mat, target_val, dtype)

    def reinitialize(self, a_coeffs=None, b_coeffs=None, first_diff_mat=None, **kwargs):
        params = kwargs
        if a_coeffs:
            params['a_coeffs'] = a_coeffs
        if b_coeffs:
            params['b_coeffs'] = b_coeffs
        super().reinitialize(params, first_diff_mat)

    def iteration_algorithm(self, params, diff_mat, i):
        # diff_mat is actually diff_vec in this case
        # it is the vector (p_n, p_{n-1}, dp_n/dx_m, dp_{n-1}/dx_m)
        # and similarly for (q_n, q_{n-1}, dq_n/dx_m, dq_{n-1}/dx_m)

        i += 1
        a_coeffs = params['a_coeffs']
        b_coeffs = params['b_coeffs']
        diff_mats_p, diff_mats_q = self._diff_mat_gen(i, a_coeffs, b_coeffs)
        diff_vec_p, diff_vec_q = diff_mat
        new_vec_p = [ new_mat * diff_vec for new_mat, diff_vec in zip(diff_mats_p, diff_vec_p) ]
        new_vec_q = [ new_mat * diff_vec for new_mat, diff_vec in zip(diff_mats_q, diff_vec_q) ]
        # new_vec_p is a list of vectors, which fit for dp_i/dx_j. The first element in each vector is p_n
        # Because it will be the same p_n for each dp_i/dx_j (only the derivatives change), we take arbitrarily
        # the first one. Same for new_vec_q
        p_i = new_vec_p[0][0,0]
        q_i = new_vec_q[0][0,0]
        try:
            pi_i = p_i/q_i
        except:
            pi_i = dec('NaN')
        params['pi'] = pi_i
        new_diff_mat = (new_vec_p, new_vec_q)

        # the vectors are (p_n, p_{n-1}, dp_n/dx_m, dp_{n-1}/dx_m)
        # so for x_j we take the j'th vector which suits the differentiation by x_j
        # and we take the 3rd element (which is [2]) for the derivative
        grad = []
        if not pi_i.is_nan():
            grad.append([ new_vec_p[j][2,0] / q_i - new_vec_q[j][2,0] * p_i/q_i**2 for j in range(len(a_coeffs)) ])
            grad.append([ new_vec_p[j][2,0] / q_i - new_vec_q[j][2,0] * p_i/q_i**2 for j in
                          range(len(a_coeffs), len(a_coeffs) + len(b_coeffs)) ])

        return (params, new_diff_mat, grad)

    def _default_diff_mat_gen(self, i, a_coeffs, b_coeffs):
        a2p = PiContFrac._array_to_polynom
        a_i = a2p(a_coeffs, i)
        b_i = a2p(b_coeffs, i)
        mat_p = []
        mat_q = []
        for j in range(len(a_coeffs)):
            mat_p.append(np.matrix(((a_i, b_i, 0, 0),), self._dtype))
        for j in range(len(b_coeffs)):
            mat_p.append(np.matrix(( (a_i, b_i, 0, 0),), self._dtype))
        for j in range(len(a_coeffs)):
            mat_q.append(np.matrix(((a_i, b_i, 0, 0),
                                (1, 0, 0, 0),
                                (0, 0, a_i, b_i),
                                (0, 0, 1, 0)), self._dtype))
        for j in range(len(b_coeffs)):
            mat_q.append(np.matrix(((a_i, b_i, 0, 0),
                                (1, 0, 0, 0),
                                (0, i**j, a_i, b_i),
                                (0, 0, 1, 0)), self._dtype))
        return (mat_p, mat_q)

    # def get_derivative(self):
    #     return self.compare_result() * self.gradient_vectors[-1]

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

    def is_pi_valid(self):
        return self.params_log['pi'][-1].is_normal()


    @staticmethod
    def _array_to_polynom(coeffs, i):
        return sum([coeffs[j] * i**j for j in range(len(coeffs))])