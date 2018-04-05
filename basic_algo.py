from decimal import Decimal as dec
import decimal
from gen_real_consts import gen_real_pi
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

        self.params_log = {p : [self._params[p]] for p in self._params}
        self.num_of_executed_iters = 0
        self.diff_mats = [self._first_diff_mat]

    def reinitialize(self, params=None, first_diff_mat=None, **kwargs):
        if params is None:
            params = kwargs
        for p in params:
            self._params[p] = params[p]

        if first_diff_mat is not None:
            self._first_diff_mat = first_diff_mat

    def gen_iterations(self, num_of_iters, iteration_algorithm, exec_finalize=True):
        self.params_log = {p : [self._params[p]] for p in self._params}
        self.num_of_executed_iters = 0
        self.diff_mats = [self._first_diff_mat]

        self.add_iterations(num_of_iters=num_of_iters, iteration_algorithm=iteration_algorithm,
                            exec_finalize=exec_finalize)

    def add_iterations(self, num_of_iters, iteration_algorithm, exec_finalize=True):
        params = dict([ (p,self.params_log[p][-1]) for p in self.params_log ])
        diff_mat = self.diff_mats[-1]
        for i in range(num_of_iters):
            try:
                params, diff_mat = iteration_algorithm(params, diff_mat, i + self.num_of_executed_iters)
                for p in params:
                    if self._logging:
                        self.params_log[p].append(params[p])
                if self._logging:
                    self.diff_mats.append(diff_mat)
            except self.SufficientAccuracy:
                break
        if hasattr(self, 'finalize_iterations') and exec_finalize:
            params, diff_mat = self.finalize_iterations(params, diff_mat, num_of_iters)

        if not self._logging:
            for p in params:
                self.params_log[p] = [params[p]]
            self.diff_mats = [diff_mat]
        self.num_of_executed_iters += num_of_iters

    def get_num_of_executed_iters(self):
        return self.num_of_executed_iters

    def compare_result(self, parameter, target_val=None, abs=True):
        if target_val is not None:
            comp_val = target_val
        elif self._target_val is not None:
            comp_val = self._target_val
        else:
            raise ValueError('No target value to compare to.')

        try:
            if abs:
                return abs(self.params_log[parameter][-1] - comp_val)
            else:
                return self.params_log[parameter][-1] - comp_val
        except:
            raise RuntimeError('Run gen_iterations first. No PI was generated!')


class PiBasicAlgo(RNNAlgo):
    def __init__(self, params, diff_mat_gen=None, first_diff_mat=None, target_val=None, logging=False, dtype=dec):
        if target_val is None:
            target_val = gen_real_pi()
        self._approach_type, self._approach_params = None, None
        super().__init__(params, diff_mat_gen, first_diff_mat, target_val, logging, dtype)

    def gen_iterations(self, num_of_iters, accuracy=0, exec_finalize=True):
        if self._approach_type in ['exp', 'over_exp'] and self._approach_params and accuracy:
            # print(self._approach_type)
            power_base, coeff = self._approach_params
            required_num_of_iters = math.ceil(abs((accuracy / power_base).ln() / power_base.ln()))
            # print('num_of_iters', required_num_of_iters)
            num_of_iters = min(num_of_iters, required_num_of_iters)
        super().gen_iterations(num_of_iters, self.iteration_algorithm, exec_finalize)

    def add_iterations(self, num_of_iters, iteration_algorithm=None, exec_finalize=True):
        if not iteration_algorithm:
            iteration_algorithm = self.iteration_algorithm
        super().add_iterations(num_of_iters, iteration_algorithm=iteration_algorithm,
                               exec_finalize=exec_finalize)

    def compare_result(self, real_pi=None, abs=True):
        return super().compare_result('pi', real_pi)

    def iteration_algorithm(self, params, diff_mat, i):
        raise NotImplementedError()