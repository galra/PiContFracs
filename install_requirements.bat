@echo off
pip install dill 1> install_requirements.log 2>&1 && echo Successfully installed dill || echo Error while instally dill
pip install numpy 1>> install_requirements.log 2>&1 && echo Successfully installed numpy || echo Error while instally numpy
pip install ordered-set 1>> install_requirements.log 2>&1 && echo Successfully installed ordered-set || echo Error while instally ordered-set
pip install PyLaTex 1>> install_requirements.log 2>&1 && echo Successfully installed PyLaTex || echo Error while instally PyLaTex
pip install scipy 1>> install_requirements.log 2>&1 && echo Successfully installed scipy || echo Error while instally scipy