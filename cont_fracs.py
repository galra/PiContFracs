import numpy as np
import basic_algo
from decimal import Decimal as dec
import decimal
from gen_real_pi import gen_real_pi
import time

# used to avoid numpy casting int into numpy.int32/64, which doesn't handle bignums
# class myint(int):
#     pass


class ZeroB(Exception):
    pass


class PiContFrac(basic_algo.PiBasicAlgo):
    def __init__(self, a_coeffs, b_coeffs, avoid_zero_b=True, check_b_threshold=False, b_threshold=10**10, diff_mat_gen=None,
                 first_diff_mat=None, with_grad=False, target_val=None, dtype=dec):
        if diff_mat_gen is None and with_grad:
            diff_mat_gen = self._default_diff_mat_with_grad_gen
        elif diff_mat_gen is None and not with_grad:
            diff_mat_gen = self._default_diff_mat_no_grad_gen
        if target_val is None:
            target_val = gen_real_pi()

        self._dtype_0 = dtype(0)
        self._dtype_1 = dtype(1)
        self._with_grad = with_grad
        params = {'a_coeffs': a_coeffs, 'b_coeffs': b_coeffs, 'pi': 1}
        if first_diff_mat is None:
            self._auto_first_diff_mat = True
            if self._with_grad:
                first_diff_mat = self._autogen_first_diff_mat_with_grad(params, dtype)
            else:
                first_diff_mat = self._autogen_first_diff_mat_no_grad(params, dtype)
        else:
            self._auto_first_diff_mat = False

        self._avoid_zero_b = avoid_zero_b
        self._check_b_threshold = check_b_threshold
        self._b_threshold = dtype(b_threshold)
        super().__init__(params, diff_mat_gen, first_diff_mat, target_val, dtype)

    def reinitialize(self, a_coeffs=None, b_coeffs=None, first_diff_mat=None, **kwargs):
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
            if self._with_grad:
                first_diff_mat = self._autogen_first_diff_mat_with_grad(params)
            else:
                first_diff_mat = self._autogen_first_diff_mat_no_grad(params)

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
        new_vec_p = (sum([l*m for l,m in zip(diff_vec_p, diff_mats_p)]), diff_vec_p[0])
        new_vec_q = (sum([l*m for l,m in zip(diff_vec_q, diff_mats_q)]), diff_vec_q[0])

        # new_vec_p is a list of vectors, which fit for dp_i/dx_j. The first element in each vector is p_n
        # Because it will be the same p_n for each dp_i/dx_j (only the derivatives change), we take arbitrarily
        # the first one. Same for new_vec_q
        # p_i = new_vec_p[0][0,0]
        # q_i = new_vec_q[0][0,0]
        p_i = new_vec_p[0]
        q_i = new_vec_q[0]
        if q_i == 0: # q_i.is_zero():
            pi_i = dec('NaN')
        else:
            pi_i = p_i / q_i
        params['pi'] = pi_i
        new_diff_mat = (new_vec_p, new_vec_q)

        # the vectors are (p_n, p_{n-1}, dp_n/dx_m, dp_{n-1}/dx_m)
        # so for x_j we take the j'th vector which suits the differentiation by x_j
        # and we take the 3rd element (which is [2]) for the derivative
        grad = []
        if not self._with_grad:
            return (params, new_diff_mat, grad)

        if not pi_i.is_nan():
            grad.append([ new_vec_p[j][2,0] / q_i - new_vec_q[j][2,0] * p_i/q_i**2 for j in range(len(a_coeffs)) ])
            grad.append([ new_vec_p[j][2,0] / q_i - new_vec_q[j][2,0] * p_i/q_i**2 for j in
                          range(len(a_coeffs), len(a_coeffs) + len(b_coeffs)) ])

        return (params, new_diff_mat, grad)

    def _default_diff_mat_with_grad_gen(self, i, a_coeffs, b_coeffs):
        a2p = PiContFrac._array_to_polynom
        a_i = a2p(a_coeffs, i)
        b_i = a2p(b_coeffs, i)
        if self._avoid_zero_b and b_i == 0:
            raise ZeroB()
        if self._avoid_zero_b and b_i == 0:
            raise ZeroB()
        if self._check_b_threshold and b_i / self.gcd(a_i, b_i) > self._b_threshold:
            raise self.SufficientAccuracy()
        mat_p = []
        mat_q = []
        for j in range(len(a_coeffs)):
            mat_p.append(np.matrix(( (a_i, b_i, 0, 0),
                                (1, 0, 0, 0),
                                (i**j, 0, a_i, b_i),
                                (0, 0, 1, 0)), self._dtype))
        for j in range(len(b_coeffs)):
            mat_p.append(np.matrix(( (a_i, b_i, 0, 0),
                                (1, 0, 0, 0),
                                (0, 0, a_i, b_i),
                                (0, 0, 1, 0)), self._dtype))
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

    def _default_diff_mat_no_grad_gen(self, i, a_coeffs, b_coeffs):
        a2p = PiContFrac._array_to_polynom
        a_i = a2p(a_coeffs, i)
        b_i = a2p(b_coeffs, i)
        if self._avoid_zero_b and b_i == 0:
            raise ZeroB()
        if self._check_b_threshold and b_i / self.gcd(a_i, b_i) > self._b_threshold:
            raise self.SufficientAccuracy()
        # mat_p = [np.matrix(( (a_i, b_i),
        #                         (1, 0)), self._dtype)]
        # mat_q = [np.matrix(((a_i, b_i),
        #                         (1, 0),), self._dtype)]
        mat_p = (a_i, b_i)
        mat_q = mat_p
        return (mat_p, mat_q)

    def _autogen_first_diff_mat_with_grad(self, params, dtype=None):
        if dtype is None:
            dtype = self._dtype
        if hasattr(self, '_dtype') and dtype == self._dtype:
            dtype_0 = self._dtype_0
            dtype_1 = self._dtype_1
        else:
            dtype_0 = dtype(0)
            dtype_1 = dtype(1)
        a_coeffs = params['a_coeffs']
        b_coeffs = params['b_coeffs']
        first_diff_mat = ([], [])
        first_diff_mat[0].extend([np.matrix(((dtype(a_coeffs[0]),), (dtype_1,),
                                             (dtype(a_coeffs[i] * (i == 0)),), (dtype_0,)), dtype=dtype)
                                  for i in range(len(a_coeffs))])
        first_diff_mat[0].extend([np.matrix(((dtype(a_coeffs[0]),), (dtype_1,),
                                             (dtype_0,), (dtype_0,)), dtype=dtype)
                                  for i in range(len(a_coeffs))])
        first_diff_mat[1].extend([np.matrix(((dtype_1,), (dtype_0,),
                                             (dtype_0,), (dtype_0,)), dtype=dtype)
                                  for i in range(len(b_coeffs) + len(a_coeffs))])
        return first_diff_mat

    def _autogen_first_diff_mat_no_grad(self, params, dtype=None):
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
        return (dtype(a_coeffs[0]), dtype_1), (dtype_1, dtype_0)
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