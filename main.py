#!/usr/bin/env python3

# TODO: check for integer roots for p \in Z_3,4[x]
# TODO: fix the ulcd constant
# TODO: fix the convergence params (usecase: the e contfracs: 'results_5303_050418.csv'
# this is a temporary execution file
import enum_params
from decimal import Decimal as dec
import time
import datetime
from gen_real_consts import gen_real_pi, gen_real_e

class MeasureRuntime():
    def __init__(self):
        self.times = []

    def start_measure(self):
        self.times.append(time.time())

    def measure_time(self):
        self.times.append(time.time())
        return self.times[-1] - self.times[-2]

# example - this is the standard continuous fraction with sixes. it works.
# mitm = enum_params.MITM(postproc_func=lambda x:x)
# mitm.build_hashtable(range_a = [[6,7], [0,1], [0,1]], range_b=[[1,2],[-4,-3],[4,5]])
# mitm.find_clicks(u_range=[0], l_range=[1], c_range=range(-4,4), d_range=[1])
# mitm.refine_clicks()
# print(mitm.filtered_params)

# example - this is for testing a specific set of parameters
# bep = enum_params.BasicEnumPolyParams()
# bep.pis_generator(range_a=[[6,7], [-1,0], [1,2]], range_b=[[5,6],[4,5],[-3,-2]])


# default is a,b in Z_2[x]
def safe_inverse(x):
    if x.is_zero():
        return dec('inf')
    else:
        return 1/x

def main():
    postproc_funcs = ['safe_inverse', 'lambda x: x', 'lambda x: x**2']

    evaluated_postproc_funcs = [ eval(ppf) for ppf in postproc_funcs ]
    measure_runtime = MeasureRuntime()
    measure_runtime.start_measure()
    mitm = enum_params.MITM(target_generator=gen_real_e, postproc_funcs=evaluated_postproc_funcs)
    print('Finished creating mitm object. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
    # a,b polynoms coefficients will be enumerated in [-2,2]
    # one can either set enum_range to set a uniform boundary to all the coefficients,
    # or set a different range to the a's coefficients and b's coefficients.
    # the given value should be either int (then the range will be [-a,a], enumeration includes both edges), or a 2-elements tuple/list
    # of the form [a,b] where a<b. enumeration includes only lower edge (b isn't included)
    mitm.build_hashtable(enum_range=2)
    print('Finished building hashtable. Runtime: %s ' % str(datetime.timedelta(seconds=measure_runtime.measure_time())))
    # for finding clicks, we enumerate u,l,c,d: (u/pi+pi/l+c)*1/d
    # TODO: add n/d instead of 1/d? equivalent to k*pi/l, technically
    # here a range should e either an int (then the enumeration is over [-i,i]), or an iterable of any type
    # (e.g. list, range object etc.)
    mitm.find_clicks(u_range=2, l_range=2, c_range=2, d_range=2)
    mitm.filter_uniq_params()
    print('Finished finding clicks. Number of clicks: %d. Runtime: %s ' %
          (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    mitm.filter_only_exp_convergence()
    print('Finished fast filtering exponential convergence. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    mitm.filter_clicks_by_approach_type()
    print('Finished full filtering exponential convergence. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    mitm.refine_clicks(accuracy=10, num_of_iterations=4000, print_clicks=False)
    print('Finished refining clicks, 8 digits accuracy, 2000 iterations. Number of clicks left: %d. Runtime: %s ' %
          (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    # mitm.filter_uniq_params()
    # print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
    #       (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    # mitm.refine_clicks(accuracy=12, num_of_iterations=10000)
    # print('Finished refining clicks, 12 digits accuracy, 10000 iterations. Number of clicks left: %d. Runtime: %s ' %
    #       (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    # mitm.filter_uniq_params()
    # print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
    #       (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    print('---REFINING CAN BE CANCELLED NOW---')
    try:
        mitm.refine_clicks(accuracy=15, num_of_iterations=2000, print_clicks=False)
        print('Finished refining clicks, 14 digits accuracy, 40000 iterations. Number of clicks left: %d. Runtime: %s ' %
              (len(mitm.get_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
        mitm.filter_uniq_params()
        print('Finished filtering unique parameters. Number of unique parameters: %d. Runtime: %s ' %
              (len(mitm.get_uniq_filtered_params()), str(datetime.timedelta(seconds=measure_runtime.measure_time()))))
    except KeyboardInterrupt as e:
        print('Canceled refining clicks, 14 digits accuracy, 40000 iterations. Runtime: %s ' %
              (str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))

    export_filename = 'results_%s.csv' % time.strftime('%M%H_%d%m%y')
    mitm.export_to_csv(export_filename, postproc_funcs)
    print('Finished saving results. Filename: %s. Runtime: %s ' %
          (export_filename, str(datetime.timedelta(seconds=measure_runtime.measure_time())) ))
    # the above enumeration is of complexity: 2**6 + 4**4, approximately 8 bits. Should be fine.

if __name__ == '__main__':
    main()