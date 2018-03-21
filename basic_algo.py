from decimal import Decimal as dec
import decimal
from gen_real_pi import gen_real_pi
import math
# import numpy as np


def set_precision(prec):
    decimal.getcontext().prec=prec


class RNNAlgo:
    class SufficientAccuracy(Exception):
        pass

    def __init__(self, params, diff_mat_gen, first_diff_mat, target_val=None, logging=False, dtype=dec):
        self._params = params
        self._diff_mat_gen = diff_mat_gen
        self._first_diff_mat = first_diff_mat
        self._target_val = target_val
        self._logging = logging
        self._dtype = dtype

        self.diff_mats = []
        self.params_log = {p:[] for p in self._params}
        self.num_of_executed_iters = 0

    def reinitialize(self, params=None, first_diff_mat=None, **kwargs):
        if params is None:
            params = kwargs
        for p in params:
            self._params[p] = params[p]

        if first_diff_mat is not None:
            self._first_diff_mat = first_diff_mat

    def gen_iterations(self, num_of_iters, accuracy, iteration_algorithm):
        self.params_log = {p : [self._params[p]] for p in self._params}
        self.diff_mats = []

        params = self._params
        diff_mat = self._first_diff_mat
        num_of_iters_div_100 = math.ceil(num_of_iters/100)
        num_of_iters = num_of_iters_div_100 * 100
        for i in range(num_of_iters_div_100):
            old_params = params
            old_diff_mat = diff_mat
            for j in range(100):
                try:
                    params, diff_mat = iteration_algorithm(params, diff_mat, i*100+j)
                    for p in params:
                        if self._logging:
                            self.params_log[p].append(params[p])
                    if self._logging:
                        self.diff_mats.append(diff_mat)
                except self.SufficientAccuracy:
                    break
            # if self.is_accuracy_met(accuracy, old_params, old_diff_mat, i):
            #     break
        if hasattr(self, 'finalize_iterations'):
            params, diff_mat = self.finalize_iterations(params, diff_mat, num_of_iters)

        if not self._logging:
            for p in params:
                self.params_log[p] = [params[p]]
            self.diff_mats = [diff_mat]
        self.num_of_executed_iters = i*100+j

    def get_num_of_executed_iters(self):
        return self.num_of_executed_iters

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
    def __init__(self, params, diff_mat_gen=None, first_diff_mat=None, target_val=None, logging=False, dtype=dec):
        if target_val is None:
            target_val = gen_real_pi()
        super().__init__(params, diff_mat_gen, first_diff_mat, target_val, logging, dtype)

    def gen_iterations(self, num_of_iters, accuracy):
        super().gen_iterations(num_of_iters, accuracy, self.iteration_algorithm)

    def compare_result(self, real_pi=None):
        return super().compare_result('pi', real_pi)

    def iteration_algorithm(self, params, diff_mat, i):
        raise NotImplementedError()