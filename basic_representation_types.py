import math


class BasicRecursiveConstRepresenter:
    """Base class for all recursive representation methods, i.e. continued fractions.
    The most basic implementation should supply only 'iteration_algorithm' to the iterations' functions.
    If 'finalize_iterations' is implemented, it's invoked at the end of gen_iterations/add_iterations."""
    # class SufficientAccuracy(Exception):
    #     pass

    def __init__(self, params, iter_calc_matrix_generator, first_iter_calc_matrix, target_val=None, logging=False):
        """params - required parameters for the representation algorithm.
        iter_calc_matrix_generator - a function that generates the matrix to step/calculate the next iteration.
        first_iter_calc_matrix - the first matrix to be used to step/calcualte the next iteration.
        target_val - a value which is the wished target to compare to.
        logging - whether all the iterations mid-results are recorded and logged (increasing runtime) or not."""
        self._params = params
        self._iter_calc_matrix_generator = iter_calc_matrix_generator
        self._first_iter_calc_matrix = first_iter_calc_matrix
        self._target_val = target_val
        self._logging = logging

        # The results will be saved here. If logging is enabled, the logging will be saved here too.
        self.params_log = {p : [self._params[p]] for p in self._params}
        self.num_of_executed_iters = 0
        self.iter_calc_matrices = [self._first_iter_calc_matrix]

    def reinitialize(self, params=None, first_iter_calc_matrix=None, **kwargs):
        """Reinitialize the instance with new parameters. Useful mostly to avoid intensive instances creation."""
        if params is None:
            params = kwargs
        for p in params:
            self._params[p] = params[p]

        if first_iter_calc_matrix is not None:
            self._first_iter_calc_matrix = first_iter_calc_matrix

    def gen_iterations(self, num_of_iters, iteration_algorithm, exec_finalize=True):
        """Resets the current state to 0 iterations and generates 'num_of_iters' iterations FROM SCRATCH,
        using 'iteration_algorithm' method for each iteration.
        If 'exec_finalize' is True and self.finalize_iterations method exists, it is being invoked at the end."""
        self.params_log = {p : [self._params[p]] for p in self._params}
        self.num_of_executed_iters = 0
        self.iter_calc_matrices = [self._first_iter_calc_matrix]

        self.add_iterations(num_of_iters=num_of_iters, iteration_algorithm=iteration_algorithm,
                            exec_finalize=exec_finalize)

    def add_iterations(self, num_of_iters, iteration_algorithm, exec_finalize=True):
        """Adds 'num_of_iters' iterations to the current state, using 'iteration_algorithm' method.
        If 'exec_finalize' is True and self.finalize_iterations method exists, it is being invoked at the end."""
        params = dict([ (p, self.params_log[p][-1]) for p in self.params_log ])
        diff_mat = self.iter_calc_matrices[-1]
        for i in range(num_of_iters):
            # try:
            params, diff_mat = iteration_algorithm(params, diff_mat, i + self.num_of_executed_iters)
            if self._logging:
                for p in params:
                    self.params_log[p].append(params[p])
                self.iter_calc_matrices.append(diff_mat)
                # except self.SufficientAccuracy:
                #     break
        if hasattr(self, 'finalize_iterations') and exec_finalize:
            params, diff_mat = self.finalize_iterations(params, diff_mat, num_of_iters)

        if not self._logging:
            for p in params:
                self.params_log[p] = [params[p]]
            self.iter_calc_matrices = [diff_mat]
        self.num_of_executed_iters += num_of_iters

    def get_num_of_executed_iters(self):
        """Returns the number of executed iterations for the current state."""
        return self.num_of_executed_iters

    def compare_result(self, parameter, target_val=None, abs=True):
        """Compares the result of 'parameter' (which should be a key to a value in self.params_log) to 'target_val'.
        In:
        parameter - a key to self.params_log to compare.
        target_val - a value to compare to. Is None and self._target_val isn't
                     (if it was defined at the instance generation), self._target_val used instead.
                     If both are None, excpetion is raised.
        abs - return absolute value of difference instead of difference
        Out: (parameter value-target_val) or absolute value of it.
        """
        if target_val is not None:
            comp_val = target_val
        elif self._target_val is not None:
            comp_val = self._target_val
        else:
            raise ValueError('No target value to compare to.')

        try:
            diff = self.params_log[parameter][-1] - comp_val
            if abs:
                return abs(diff)
            else:
                return diff
        except:
            raise RuntimeError('Run gen_iterations first. No result was generated!')


class BasicContFracAlgo(BasicRecursiveConstRepresenter):
    """This is an extension of BasicRecursiveConstRepresenter that adds:
    * support for an approach type/parameters infrastructure (poly, exp, super_exp)
    * feeds automatically the iterations' functions with self.iteration_algorithm method
    * compare_result compares to the 'contfrac_res' parameter automatically"""
    def __init__(self, params, iter_calc_matrix_generator=None, first_iter_calc_matrix=None, target_val=None,
                 logging=False):
        """params - required parameters for the representation algorithm.
        iter_calc_matrix_generator - a function that generates the matrix to step/calculate the next iteration.
        first_iter_calc_matrix - the first matrix to be used to step/calcualte the next iteration.
        target_val - a value which is the wished target to compare to.
        logging - whether all the iterations mid-results are recorded and logged (increasing runtime) or not.
        """
        # These will be used to save the type and parameters of approach to convergence value, if evaluated.
        self._approach_type, self._approach_params = None, None
        super().__init__(params, iter_calc_matrix_generator, first_iter_calc_matrix, target_val, logging)

    def gen_iterations(self, num_of_iters, accuracy=0, exec_finalize=True):
        """Resets the current state to 0 iterations and generates 'num_of_iters' iterations FROM SCRATCH.
        If the approach type is evaluated and supported, and less then 'num_of_iters' interations are required to
        achieve an error of at most 'accuracy', if 'accuracy' isn't 0.
        If 'exec_finalize' is True and self.finalize_iterations method exists, it is being invoked at the end."""
        # Calculates the required number of iterations if required, for fast convergence.
        if self._approach_type in ['exp', 'super_exp'] and self._approach_params and accuracy:
            power_base, coeff = self._approach_params
            required_num_of_iters = math.ceil(abs((accuracy / power_base).ln() / power_base.ln()))
            num_of_iters = min(num_of_iters, required_num_of_iters)
        super().gen_iterations(num_of_iters, self.iteration_algorithm, exec_finalize)

    def add_iterations(self, num_of_iters, iteration_algorithm=None, exec_finalize=True):
        """Adds 'num_of_iters' iterations to the current state, using 'iteration_algorithm' method. If no iteration
        algorithm is supplied, self.iteration_algorithm is used as default.
        If 'exec_finalize' is True and self.finalize_iterations method exists, it is being invoked at the end."""
        if not iteration_algorithm:
            iteration_algorithm = self.iteration_algorithm
        super().add_iterations(num_of_iters, iteration_algorithm=iteration_algorithm,
                               exec_finalize=exec_finalize)

    def compare_result(self, target_val=None, abs=True):
        """Compare the last result (from the last evaluated iteration) to 'target_val' and return the difference.
        If 'abs' is True, returns the absolute value of the difference."""
        return super().compare_result('contfrac_res', target_val)

    def iteration_algorithm(self, params, promo_mat, i):
        """params - inner algorithm required parameters.
        promo_mat - promotion matrix (or similar, in the broad sense), needed to generate the next iteration.
        i - the index of the generated iteration."""
        raise NotImplementedError()