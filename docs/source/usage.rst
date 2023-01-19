Usage
=====

Installation
------------

Python:

To use Assist, first install it using pip:

.. code-block:: console

   (.venv) $ pip install assist

If installation fails, make sure that wheel is installed in your virtual environment:

.. code-block:: console

   (.venv) $ pip install wheel

C:

Navigate to the parent directory that holds the rebound folder (to install in a custom folder, set the REB_DIR environment variable as shown below). Then in a terminal:


.. code-block:: console

	git clone https://github.com/mholman/assist.git

(install git if you don’t have it). 

By default, the code looks for the ephemeris files in the assist/data directory, unless the environment variables are set differently, as below.

You can download the DE441 ephemeris file here:

https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441

You can find the SPK files and documentation here:

ftp://ssd.jpl.nasa.gov/pub/eph/small_bodies/asteroids_de441/sb441-n16.bsp

For custom install locations, set the following environment variables as appropriate:

.. code-block:: console

	export REB_DIR=/Users/mholman/rebound

	export ASSIST_DIR=/Users/mholman/assist

You access the ephemeris files adding equivalent versions of these lines to your .bash_profile .

.. code-block:: console

	export JPL_PLANET_EPHEM=/Users/mholman/assist/data/linux_m13000p17000.441

	export JPL_SB_EPHEM=/Users/mholman/assist/data/sb441-n16.bsp




