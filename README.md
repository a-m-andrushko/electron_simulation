This repository contains a simulation of a quantum particle (in this example, an electron) in a 1D box with an oscillating electric field.

A detailed description on how it was built is laid out in _instructions.pdf_ (in Polish).

_electron_simulation.cpp_ is the main file where all calculations are done, it can be launched using any up-to-date C++ compiler and needs to be located in the same folder as _parameters.txt_; the latter file consists of examplary parameters for a basic simulation to be performed. Running _electron_simulation.cpp_ will produce 3 text files -- with wave function states (_wave_function.txt_); particle's norm, average position, and energy (_quantities.txt_); and particle's probability density (_density.txt_). Those are necessary for further analysis.

In order to visualise evolution of probability density via _gnuplot_, _animate.zip_ needs to be unpacked. In _gnuplot.rot_, the name of the processed file needs to be adjusted. In _animate.dem_, number of iterations must be the same as the number of iterations in the evolution loop in _electron_simulation.cpp_. By typing in terminal _gnuplot animate.dem_, the simulation will be launced.

Finally, _electron_simulation.ipynb_ is a Python code containing analysis of the programme based on the points listed in _instructions.pdf_.
