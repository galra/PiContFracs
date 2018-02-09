#!/usr/bin/env python3

# this is a temporary execution file
import enum_params
from decimal import Decimal as dec

# example - this is the standard continuous fraction with sixes. it works.
# mitm = enum_params.MITM(postproc_func=lambda x:x)
# mitm.build_hashtable(range_a = [[6,7], [0,1], [0,1]], range_b=[[1,2],[-4,-3],[4,5]])
# mitm.find_clicks(u_range=[0], l_range=[1], c_range=range(-4,4), d_range=[1])
# mitm.refine_clicks()
# print(mitm.filtered_params)

# default is a,b in Z_2[x]
def safe_inverse(x):
    if x.is_zero():
        return dec('inf')
    else:
        return 1/x

mitm = enum_params.MITM(postproc_func=safe_inverse)
# a,b polynoms coefficients will be enumerated in [-2,2]
# one can either set enum_range to set a uniform boundary to all the coefficients,
# or set a different range to the a's coefficients and b's coefficients.
# the given value should be either int (then the range will be [-a,a], enumeration includes both edges), or a 2-elements tuple/list
# of the form [a,b] where a<b. enumeration includes only lower edge (b isn't included)
mitm.build_hashtable(enum_range=2)
# for finding clicks, we enumerate u,l,c,d: (u/pi+pi/l+c)*1/d
# TODO: add n/d instead of 1/d? equivalent to k*pi/l, technically
# here a range should e either an int (then the enumeration is over [-i,i]), or an iterable of any type
# (e.g. list, range object etc.)
mitm.find_clicks(u_range=4, l_range=4, c_range=4, d_range=4)
mitm.refine_clicks()
print(mitm.filtered_params)

# the above enumeration is of complexity: 2**6 + 4**4, approximately 8 bits. Should be fine.