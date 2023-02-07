.. _python_quickstart:

Quickstart (Python)
===================

Installation
------------
It's easiest to install ASSIST into a python virtual environment. If you already have a virtual environment 
or do not want to use one, you can skip this step. Otherwise, run the following command in an empty directory. 
They will setup and activate a new virtual environment in a directory.

.. code-block:: console

	python3 -m venv venv
	source venv/bin/activate

Now we can install numpy, REBOUND, and ASSIST:

.. code-block:: console	

	pip install numpy
	pip install rebound 
	pip install assist

To use use ASSIST, you also need to download ephemeris data files. One file 
for planet ephemeris and another suplementary file for asteroid ephemeris. 
The following commands download these files with curl. You can also manually 
download them using your browser. Note that these are large files, almost 1GB in size.

.. code-block:: console

	mkdir data
	curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440
	curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o data/sb441-n16.bsp

Now you can try out if assist works.

.. code-block:: console

	python3

	>>> import assist
	>>> ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")
	>>> print(ephem.jd_ref)
	>>> ephem.get_particle("Earth", 0)

You should see the default reference Julian date (2451545.0) and the position of the Earth at that time printed on the screen.

At this point you are done and can skip to the Quick Start Guide below.

For a more complete installation, i.e., if you want any of the following: 

* Source code
* The example files so that you can modify them locally.
* To also use the C version
 
First follow the installation instructions for the C version in :ref:`c_quickstart`.
Then, to install the Python version from this repository, navigate to the `assist` directory and
(you'd also do this to install the Python version after modifying any of the C code)::

    pip install -e ./

.. _python_qs:

Quick Start Guide
-----------------

(See also: ipython notebook: assist/jupyter_examples/Getting started.ipynb on github)

You begin by importing rebound and assist


.. code:: python

    import rebound
    import assist

You then import the downloaded ephemeris data, using the path to the data directory


.. code:: python

    ephem = assist.Ephem("data/linux_p1550p2650.440", "data/sb441-n16.bsp")

Then  set up your REBOUND simulation with initial positions in AU and initial velocities in AU/day 
(Here we use the initial conditions for asteroid 3666 Holman) 

.. code:: python

    import rebound
    sim = rebound.Simulation()
    holman_initial = rebound.Particle(
    x=3.338875348598862E+00, y=-9.176518412197102E-01, z=-5.038590741719294E-01, 
    vx=2.805663364339457E-03, vy=7.550408665778840E-03, vz=2.980028207875623E-03)
    sim.add(holman_initial)

Next, set the initial simulation time corresponding to the initial conditions above. 
The above initial conditions are valid at 2458849.5 Julian Days (2020-Jan-01). In ASSIST, 
we measure time relative to the jd_ref parameter in the ephemeris structure.

.. code:: python

    sim.t = 2458849.5 - ephem.jd_ref

attach assist to the simulation:


.. code:: python

    ax = assist.Extras(sim, ephem)

and integrate forward to a desired final time (here, 10000 days):

.. code:: python

    t_final = sim.t + 100000
    ax.integrate_or_interpolate(t_final)

display the courdinates of the asteroid at t_final:

.. code:: python
 
   sim.particles[0].xyz
