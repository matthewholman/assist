.. image:: https://img.shields.io/badge/assist-v1.0.0-green.svg?style=flat
    :target: https://assist.readthedocs.org

# ASSIST

ASSIST is a package supporting high-accuracy ephemeric calculations for small Solar System bodies in REBOUND.


## Installation (Python)

It's easiest to install ASSIST into a python virtual environment. The following bash commands install and activate a new virtual environment in a directory. If you already have a virtual environment or do not want to use one, you can skip this step.

    python3 -m venv venv
    source venv/bin/activate

Now we can install numpy, REBOUND, and ASSIST:

    pip install numpy
    pip install rebound 
    pip install assist

To use use ASSIST, you also need to download ephemeris data file, specifically one file for planet ephemeris (DE440) and a suplementary file for asteroid ephemeris. The following commands download these files with curl. You can also manually download them using your browser. Note that these are large files, almost 1GB in size.

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
