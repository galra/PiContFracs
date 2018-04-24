import numpy as np
from pylatex import Document, Section, Math, Alignat
import os


def generate_latex(eqns=[]):
    doc = Document()

    with doc.create(Section('Automatic Conjectures')):
        doc.append('These are the conjectures detected by the algorithm.')

        for eqn in eqns:
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(eqn)

    doc.generate_pdf('results', clean_tex=False)


def latex_cont_frac(a, b, current_iteration=''):
    len_a = len(a)
    len_b = len(b)

    if current_iteration == '':
        current_iteration = str(a[len_a - 1]) + ' + \dots'

    if len_a > 1:
        new_iteration = r'{0} + \frac{{ {1} }} {{ {2} }}'.format(a[len_a - 2], b[len_b - 1], current_iteration)
        return latex_cont_frac(a[:len_a - 1], b[:len_b - 1], new_iteration)
    else:
        return current_iteration
    

if __name__ == '__main__':
    # generate_latex([r'1=2', r'2=3'])
    print(latex_cont_frac([2, 4, 1, 5], [3, 2, 7, 6]))
