[![Version](https://img.shields.io/badge/assist-v1.0.1b7-green.svg?style=flat)](https://assist.readthedocs.org)

# ASSIST

ASSIST is a package supporting high-accuracy ephemeric calculations for small Solar System bodies in REBOUND.


## Installation (Python)

It's easiest to install ASSIST into a python virtual environment. If you already have a virtual environment or do not want to use one, you can skip this step. Otherwise, run the following command in an empty directory. They will setup and activate a new virtual environment in a directory. 

    python3 -m venv venv
    source venv/bin/activate

Now we can install numpy, REBOUND, and ASSIST:

    pip install numpy
    pip install rebound 
    pip install assist

To use use ASSIST, you also need to download ephemeris data files. One file for planet ephemeris and another suplementary file for asteroid ephemeris. The following commands download these files with curl. You can also manually download them using your browser. Note that these are large files, almost 1GB in size.

    mkdir data
    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440
    curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o data/sb441-n16.bsp

Now you can try out if assist works.

    python3

    >>> import assist
    >>> ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
    >>> print(ephem.jd_ref)
    >>> ephem.get_particle("Earth", 0)

You should see the default reference Julian date (2451545.0) and the position of the Earth at that time printed on the screen.

## Installation (C)

To install the C version of ASSIST, first clone the REBOUND and then the ASSIST repositories. In an empty directory, run:

    git clone https://github.com/hannorein/rebound.git
    git clone https://github.com/matthewholman/assist.git

To use use ASSIST, you also need to download ephemeris data files. One file for planet ephemeris and another suplementary file for asteroid ephemeris. The following commands download these files with curl. You can also manually download them using your browser. Note that these are large files, almost 1GB in size.

    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o assist/data/linux_p1550p2650.440
    curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o assist/data/sb441-n16.bsp

For some of the examples, you will also need the planet ephemeris file with an extended coverage. Note that this file is 2.6GB in size.

    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441 -o assist/data/linux_m13000p17000.441

Next, go to one of the example directories and compile the problem file. This will also trigger the installation of the REBOUND and ASSIST shared libraries.

    cd assist/examples/plain_interface
    make

Now, you're ready to run the example with:

    ./rebound


