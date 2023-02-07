.. ASSIST documentation master file, created by
   sphinx-quickstart on Mon Feb  6 13:12:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome
=======

ASSIST is a software package for ephemeris-quality integrations of test particles. 
ASSIST is an extension of the REBOUND framework and makes use of its IAS15 integrator 
to integrate test particle trajectories in the field of the Sun, Moon, planets, and 16 
massive asteroids, with the positions of the masses extracted from the JPL DE441
ephemeris and its associated asteroid perturber file. The package incorporates the most
significant gravitational harmonics and general relativistic corrections. ASSIST also
accounts for position- and velocity-dependent non-gravitational effects according to
the Marsden (1973) model. All components in the equations of motion have been verified 
to machine precision in a term-by-term comparison with output from JPL's small body 
integrator. The first order variational equations are included for all terms to support 
orbit fitting and covariance mapping. This framework is meant to provide an open-source 
package written in a modern language to enable high-precision orbital analysis and 
science by the small body community.

If you are using ASSIST for the first time check out the quickstart guides below with installation instructions.

If you clone the repository at `https://github.com/mholman/assist <https://github.com/mholman/assist>`_ you can 
load and run all the jupyter notebook examples locally (under assist/jupyter_examples) 
as well as the C examples 
(under assist/examples. In the terminal you can run the example in each folder with ``make clean && make && ./rebound``). 

Get Started!
------------

ASSIST is written in C, but we also provide a convenient Python wrapper.

:ref:`python_quickstart`

:ref:`c_quickstart`

.. toctree::
    :numbered:
    :hidden:

    python_quickstart
    c_quickstart

