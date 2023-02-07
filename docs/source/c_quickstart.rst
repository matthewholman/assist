.. _c_quickstart:

Quickstart (C)
==============

Installation
------------

To install the C version of ASSIST, first clone the REBOUND and then the ASSIST repositories. In an empty directory, run:

.. code-block:: console

    git clone https://github.com/hannorein/rebound.git
    git clone https://github.com/matthewholman/assist.git

To use use ASSIST, you also need to download ephemeris data files. 
One file for planet ephemeris and another suplementary file for asteroid ephemeris. 
The following commands download these files with curl. 
You can also manually download them using your browser. Note that these are large files, almost 1GB in size.

.. code-block:: console

	curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o assist/data/linux_p1550p2650.440
	curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o assist/data/sb441-n16.bsp

For some of the examples, you will also need the planet ephemeris file with an extended coverage. Note that this file is 2.6GB in size.

.. code-block:: console

	curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441 -o assist/data/linux_m13000p17000.441


.. _c_qs:

Quick Start Guide
-----------------

Once you have installed REBOUND and ASSIST and downloaded the ephemeris files to assist/data,
go to one of the example directories and compile the problem file in the directory. 
(This will also trigger the installation of the REBOUND and ASSIST shared libraries.)

.. code-block:: console

	cd assist/examples/asteroid
	make

And run the example with:

.. code-block:: console

	./rebound
