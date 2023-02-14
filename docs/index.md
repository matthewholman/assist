# Welcome to ASSIST

ASSIST is a software package for ephemeris-quality integrations of test particles. ASSIST is an extension of the [REBOUND framework](https://github.com/hannorein/rebound) and makes use of its [IAS15 integrator](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.1424R/abstract) to integrate test particle trajectories in the field of the Sun, Moon, planets, and 16 massive asteroids, with the positions of the masses coming from the JPL DE441 ephemeris and its associated asteroid perturber file. The package incorporates the most significant gravitational harmonics and general relativistic corrections. ASSIST also accounts for position- and velocity-dependent non-gravitational effects according to the [Marsden (1973) model](https://ui.adsabs.harvard.edu/abs/1973AJ.....78..211M/abstract). All components in the equations of motion have been verified to machine precision in a term-by-term comparison with output from JPL's small body integrator. The first order variational equations are included for all terms to support orbit fitting and covariance mapping. This framework is meant to provide an open-source package written in a modern language to enable high-precision orbital analysis and science by the small body community.

## License
ASSIST is open source, freely distributed under the [GNU General Public license, version 3](https://github.com/matthewholman/blob/main/LICENSE).

## Contributors

* Matthew J. Holman, Center for Astrophysics | Harvard & Smithsonian, <mholman@cfa.harvard.edu>
* Arya Akmal, Montgomery College, Rockville
* Davide Farnocchia, Jet Propulsion Laboratory, California Institute of Technology 
* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Matthew J. Payne, Center for Astrophysics | Harvard & Smithsonian
* Robert Weryk, University of Western Ontario
* Dan Tamayo, Harvey Mudd College, <dtamayo@hmc.edu>
* David M. Hernandez, Center for Astrophysics | Harvard & Smithsonian 


