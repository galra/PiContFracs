import numpy as np

from pylatex import Document, Section, Math, Alignat
import os


def generate_latex(eqns=[]):
    doc = Document()

    with doc.create(Section('Potential identities')):
        doc.append('These are potential identities detected by the algorithm.')

        u, l, c, d = -2, 1, 3, -4
        a0, a1, a2, a3 = 2, 2, 2, 2
        b0, b1, b2 = 1, 1, 1
        # (u/target_value + target_value/l + c) / d

        for eqn in eqns:
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                # lhs = r'\frac{{ \frac{{ {0} }}{{\pi}} + \frac{{ \pi }} {{ {1} }} + {2} }} {{ {3} }}'.format(u, l, c, d)
                # rhs = r'{0} + \frac{{ {1} }} {{ {2} + \frac{{ {3} }} {{ {4} + \frac{{ {5} }} {{ {6} + \dots }} }} }}'.format(a0, b1, a1, b2, a2, b2, a3)
                # eqn = lhs + r'&=' + rhs + r' \\'
                agn.append(eqn)

    doc.generate_pdf('full', clean_tex=False)


if __name__ == '__main__':
    generate_latex([r'1=2', r'2=3'])
