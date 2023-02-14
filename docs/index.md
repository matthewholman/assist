# ASSIST

ASSIST is a software package for ephemeris-quality integrations of test particles. ASSIST is an extension of the [REBOUND framework](https://github.com/hannorein/rebound) and makes use of its [IAS15 integrator](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.1424R/abstract) to integrate test particle trajectories in the field of the Sun, Moon, planets, and 16 massive asteroids, with the positions of the masses coming from the JPL DE441 ephemeris and its associated asteroid perturber file. The package incorporates the most significant gravitational harmonics and general relativistic corrections. ASSIST also accounts for position- and velocity-dependent non-gravitational effects according to the [Marsden (1973) model](https://ui.adsabs.harvard.edu/abs/1973AJ.....78..211M/abstract). All components in the equations of motion have been verified to machine precision in a term-by-term comparison with output from JPL's small body integrator. The first order variational equations are included for all terms to support orbit fitting and covariance mapping. This framework is meant to provide an open-source package written in a modern language to enable high-precision orbital analysis and science by the small body community.



## Installation

=== "Python"
    It's easiest to install ASSIST into a python virtual environment. If you already have a virtual environment or do not want to use one, you can skip this step. Otherwise, run the following command in an empty directory. They will setup and activate a new virtual environment in a directory. 
    ```bash
        python3 -m venv venv
        source venv/bin/activate
    ```
    Now we can install numpy, REBOUND, and ASSIST:
    ```bash
        pip install numpy
        pip install rebound 
        pip install assist
    ```

=== "C"
    To use the C version of ASSIST, first clone the REBOUND and then the ASSIST repository. In an empty directory, run:
    ```bash
        git clone https://github.com/hannorein/rebound.git
        git clone https://github.com/matthewholman/assist.git
    ```bash


To use ASSIST, you also need to download ephemeris data files. One file for planet ephemeris and another suplementary file for asteroid ephemeris. The following commands download these files with curl. You can also manually download them using your browser. Note that these are large files, almost 1GB in size. You can store these files anywhere. If you're installing REBOUND from the repository, a good place to put them is in a new directory with the name `data`.

```bash
    mkdir data
    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440
    curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o data/sb441-n16.bsp
```
    
For some of the examples, you will also need the planet ephemeris file with an extended coverage.

!!! Warning

    This file is 2.6GB in size.

```bash
    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441 -o assist/data/linux_m13000p17000.441
```

Now you can try out if assist works.

=== "Python"
    Start a python interpreter, or open a new notebook if you're planning to use Jupyter notebooks.
    ```bash
        python3
    ```
    Then run the following code:
    ```python
        import assist
        ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
        print(ephem.jd_ref)
        ephem.get_particle("Earth", 0)
    ```
    You should see the default reference Julian date (2451545.0) and the position of the Earth at that time printed on the screen.

=== "C"
    Go to one of the example directories and compile the problem file. This will also trigger the installation of the REBOUND and ASSIST shared libraries.

    ```bash
        cd assist/examples/plain_interface
        make
    ```

    Now, you're ready to run the example with:

    ```bash
        ./rebound
    ```


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

