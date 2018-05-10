import basic_algo
from cont_fracs import ContFrac
from decimal import Decimal as dec
import numpy as np
import random
from math import floor, ceil
import itertools

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class GradientDescentBasicAlgo:
    def __init__(self, enforce_Z=False, max_num_of_iterations = 7000, threshold=None, step_size='dynamic'):
        # lattice parameters
        self.enforce_Z = enforce_Z

        # Grad-Des parameters
        if threshold is None:
            threshold = dec(10)**dec(-4)
        self.threshold = threshold
        if step_size is None:
            step_size = [dec(0.001)] * 2
        self.step_size = step_size
        self.max_num_of_iterations = max_num_of_iterations

        # set decimal precision
        basic_algo.set_precision(100)

    def find_params(self, a_coeffs, b_coeffs, show_progress=True):
        """Runs gradient-descent with the given parameters. step-size can be either dynamic of a number.
        dynamic is good"""
        # init initial values
        ITERATIONS_CHUNCK_SIZE = 50

        round_lottery = lambda : [round, floor, ceil][random.random() < 0.0]

        a_coeffs = [ dec(m) for m in a_coeffs ]
        b_coeffs = [ dec(m) for m in b_coeffs ]
        if self.enforce_Z:
            a_coeffs_param = [ dec(round_lottery()(m)) for m in a_coeffs ]
            b_coeffs_param = [ dec(round_lottery()(m)) for m in b_coeffs ]
        else:
            a_coeffs_param = a_coeffs
            b_coeffs_param = b_coeffs

        pcf = ContFrac(a_coeffs_param, b_coeffs_param)
        pcf.gen_iterations(ITERATIONS_CHUNCK_SIZE)
        iter_num = 1
        dec_0 = dec(0)

        while (abs(pcf.compare_result()) > self.threshold) and (iter_num < self.max_num_of_iterations):
            grad = pcf.get_derivative()
            if self.step_size == 'dynamic':
                    dx = dec(0.1) * np.array((max(abs(grad[0])).log10().max(dec(1)),
                                                  max(abs(grad[1])).log10().max(dec(1))))
            dx = dx.tolist()
            # minus sign because we want to go downhill instead of uphill
            a_coeffs_diff = np.multiply(grad[0], dx[0])
            b_coeffs_diff = np.multiply(grad[1], dx[1])
            a_coeffs += a_coeffs_diff
            b_coeffs += b_coeffs_diff
            if self.enforce_Z:
                a_coeffs_param = [ dec(round_lottery()(m)) for m in a_coeffs ]
                b_coeffs_param = [ dec(round_lottery()(m)) for m in b_coeffs ]
            else:
                a_coeffs_param = [ dec(m) for m in a_coeffs ]
                b_coeffs_param = [ dec(m) for m in b_coeffs ]
            pcf.reinitialize(a_coeffs_param, b_coeffs_param)
            pcf.gen_iterations(ITERATIONS_CHUNCK_SIZE)
            iter_num += 1
            if iter_num % 100 == 0 and show_progress:
                print('\r%d, error=%f' % (iter_num, pcf.compare_result()), end='')
        if show_progress:
            print('')

        if iter_num >= self.max_num_of_iterations and show_progress:
            print('Iterations limit reached. Aborting.')
        elif show_progress:
            print('Result distance: %s' % abs(pcf.compare_result()))

        if iter_num >= self.max_num_of_iterations:
            return
        return (a_coeffs_param, b_coeffs_param)
#
# def gen_pi_map(x1_range, y1_range, pi0=None, resolution=100, iterations=10):
#     """Plots a map of the function iterations results over a square range"""
#     map_data = {'x1': [], 'y1': [], 'pi': []}
#     if pi0 is None:
#         pi0 = dec(2) + dec(2).sqrt()
#     ba = basic_algo.PiBasicAlgo(x1=dec(1), y1=dec(1), pi0=dec(1))
#
#     if x1_range[0] < 0 or y1_range[0] < 0:
#         raise ValueError('Range most be of positive values!')
#
#     if x1_range[0] >= x1_range[1] or y1_range[0] >= y1_range[1]:
#         raise ValueError('Range boundaries are disordered!')
#
#     x_space = np.linspace(x1_range[0], x1_range[1], resolution)
#     y_space = np.linspace(y1_range[0], y1_range[1], resolution)
#     map_data['x1'], map_data['y1'] = np.meshgrid(x_space, y_space)
#     map_data['pi'] = np.array(map_data['x1'])
#     for x1_i in range(len(x_space)):
#         for y1_i in range(len(y_space)):
#             ba.reinitialize(dec(x_space[x1_i]), dec(y_space[y1_i]), dec(pi0))
#             ba.gen_iterations(iterations)
#             map_data['pi'][x1_i,y1_i] = (ba.compare_result()**dec(2)).log10()
#
#     return map_data
#
# def gen_pi_path(x1, y1, pi0=None, iterations=10):
#     """Calculates the whole path of x,y,pi during the iterations process
#     :param x1:
#     :param y1:
#     :param pi0:
#     :param iterations:
#     :return:
#     """
#     if pi0 is None:
#         pi0 = dec(2) + dec(2).sqrt()
#     ba = basic_algo.PiBasicAlgo(x1=x1, y1=y1, pi0=pi0)
#     ba.gen_iterations(iterations)
#     return {'x': ba.xs, 'y': ba.ys, 'pi': ba.pis}
#
# def draw_map(map_data, fig):
#     """Plots a map to a given figure"""
#     # ax = fig.gca(projection='3d')
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot_wireframe(map_data['x1'], map_data['y1'], map_data['pi']) #, rstride=1, cstride=1)
#
# def draw_path(path_data, fig):
#     """Plots a given path of x,y,pi (path_data is a dictionary) to a given figure"""
#     ax = fig.gca(projection='3d')
#     ax.plot(path_data['x'], path_data['y'], path_data['pi']) #, label='parametric curve')
#
# def draw_all(x_range=[0.1, 2], y_range=[0.1, 2], starting_points=None,
#              pi0=None, resolution=100, iterations=10, zero_threshold=0.001):
#     if pi0 is None:
#         pi0 = dec(2) + dec(2).sqrt()
#     if starting_points is None:
#         x_space = np.linspace(x_range[0], x_range[1], resolution)
#         y_space = np.linspace(y_range[0], y_range[1], resolution)
#         starting_points = itertools.product(x_space, y_space)
#
#     fig = plt.figure()
#     map_data = gen_pi_map(x_range, y_range, pi0, resolution, iterations)
#     draw_map(map_data, fig)
#
#     # find the path coordinates
#     path_xy = [ (x_space[i], y_space[j])
#                 for i,j in itertools.product(range(len(x_space)), range(len(y_space)))
#                 if map_data['pi'][i,j] < -3.5 ]
#
#     # zero_path = [ (x_space[i], y_space[j], map_data['pi'][i,j])
#     #               for i,j in i,j in itertools(range(len(x_space)), range(len(y_space)))
#     #               if map_data['pi'][i,j]) < zero_threshold ]
#     path_x, path_y = [ x for x,y in path_xy ], [ y for x,y in path_xy ]
#     fig_path = plt.figure()
#     ax_path = fig_path.gca()
#     ax_path.plot(path_x, path_y, 'b.')
#
#     fig.show()
#     fig_path.show()
#     return (path_x, path_y)
#
# def array_to_matlab(array_name, array):
#     """Converts an array a to string to paste to matlab"""
#     return '%s = [%s];' % (array_name, ' '.join([ i.__repr__() for i in array ]))
#
