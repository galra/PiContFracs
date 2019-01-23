@echo off
rem pip install dill 1> install_requirements.log 2>&1 && echo Successfully installed dill || echo Error while instally dill
rem pip install numpy 1>> install_requirements.log 2>&1 && echo Successfully installed numpy || echo Error while instally numpy
rem pip install ordered-set 1>> install_requirements.log 2>&1 && echo Successfully installed ordered-set || echo Error while instally ordered-set
rem pip install PyLaTex 1>> install_requirements.log 2>&1 && echo Successfully installed PyLaTex || echo Error while instally PyLaTex
rem pip install scipy 1>> install_requirements.log 2>&1 && echo Successfully installed scipy || echo Error while instally scipy
rem pip install decmath 1>> install_requirements.log 2>&1 && echo Successfully installed decmath || echo Error while instally decmath
rem pip install python-flint 1>> install_requirements.log 2>&1 && echo Successfully installed flint || echo Error while instally flint
pip install -r requirements.txt
echo download perl, add it to PATH