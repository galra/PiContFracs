from decimal import Decimal as dec
import decimal
from gen_real_pi import gen_real_pi
import numpy as np


def set_precision(prec):
    decimal.getcontext().prec=prec


class RNNAlgo:
    def __init__(self, params, diff_mat_gen, first_diff_mat, target_val=None, dtype=dec):
        self._params = params
        self._diff_mat_gen = diff_mat_gen
        self._first_diff_mat = first_diff_mat
        self._target_val = target_val
        self._dtype = dtype

        self.gradient_vectors = []
        self.diff_mats = []
        self.params_log = {p:[] for p in self._params}

    def reinitialize(self, params=None, first_diff_mat=None, **kwargs):
        if params is None:
            params = kwargs
        for p in params:
            self._params[p] = params[p]

        if first_diff_mat is not None:
            self._first_diff_mat = first_diff_mat

    def gen_iterations(self, num_of_iters, iteration_algorithm):
        self.params_log = {p : [self._params[p]] for p in self._params}
        self.gradient_vectors = []
        self.diff_mats = []

        params = self._params
        diff_mat = self._first_diff_mat
        for i in range(num_of_iters):
            params, diff_mat, gradient_vec = iteration_algorithm(params, diff_mat, i)
            for p in params:
                self.params_log[p].append(params[p])
            self.diff_mats.append(diff_mat)
            self.gradient_vectors.append(gradient_vec)

    def get_derivative(self):
        return np.multiply(self.compare_result(), self.gradient_vectors[-1])

    def compare_result(self, parameter, target_val=None):
        if target_val is not None:
            comp_val = target_val
        elif self._target_val is not None:
            comp_val = self._target_val
        else:
            raise ValueError('No target value to compare to.')

        try:
            return self.params_log[parameter][-1] - comp_val
        except:
            raise RuntimeError('Run gen_iterations first. No PI was generated!')


class PiBasicAlgo(RNNAlgo):
    def __init__(self, params, diff_mat_gen=None, first_diff_mat=None, target_val=None, dtype=dec):
        if target_val is None:
            target_val = gen_real_pi()
        super().__init__(params, diff_mat_gen, first_diff_mat, target_val, dtype)

    def gen_iterations(self, num_of_iters):
        super().gen_iterations(num_of_iters, self.iteration_algorithm)

    def compare_result(self, real_pi=None):
        return super().compare_result('pi', real_pi)

    def iteration_algorithm(self, params, diff_mat, i):
        raise NotImplementedError()