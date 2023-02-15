# Forces

ASSIST includes several forces when calculating trajectories. By default the gravitational forces from the Sun, Moon, planets, and 16 massive asteroids, forces from the most significant gravitational harmonics, general relativistic corrections, and forces which account for position- and velocity-dependent non-gravitational effects are included. You turn each force contribution on or off. This can be used to determine how important each effect is and to increase efficiency.

## Configuring forces

The configuration of the forces happens in the assist extras structure after you attach it to a simulation (see the examples on how to do that). 
The following code removes the higher order gravitational harmonics of the Earth and Sun from the force calculation.


=== "Python"

    ```python
    extras = assist.Extras(sim, ephem)
    forces = extras.forces
    forces.remove("EARTH_HARMOMNICS")
    forces.remove("SUN_HARMOMNICS")
    extras.forces = forces
    ```
    
    !!! Warning
        `forces` is just a list of strings. 
        If you want to remove a force, make sure you first make a copy of the list, remove an item, and then set the copy of the list to the property. 
        Doing something like `extras.forces.remove("EARTH_HARMONICS")` does not work.


=== "C"

    ```c
    struct assist_extras* ax = ...
    ax->forces ^= ASSIST_FORCE_EARTH_HARMONICS;
    ax->forces ^= ASSIST_FORCE_SUN_HARMONICS;
    ```
    !!! Info 
        Rather than having one parameter for each force, the `forces` parameter is an integer which acts as a bitfield. 
        For example, by default the bit that corresponds to earth harmonics in `forces` is turned on (set to 1). 
        The `^=` operator is an in-place bitwise xor operator which flips the earth harmonics bit to 0.


## Non-gravitational forces

We use the model by [Marsen 1973](https://ui.adsabs.harvard.edu/abs/1991AJ....102.1539M) for non-gravitational forces such as radiation pressure, Yarkovsky effect, and outgassing. 
For this force to have any effect on your particles, you need to provide a list of three parameters \(A_1\), \(A_2\), and \(A_3\) for each particle. These are the coefficients, in units of \(au/d^2\), for the radial, tangential, and normal components of the non-gravitational forces. The following code demonstrates how this can be done.

=== "Python"
    
    !!! Note inline end
        This has to be a numpy array.

    ```python
    import numpy as np
    ...
    extras = assist.Extras(sim, ephem)
    ex.particle_params = np.array([4.99e-13, -2.90e-14, 0.0])
    ```
    

=== "C"

    ```c
    struct assist_extras* ax = ...
    double params[] = {4.99e-13, -2.90e-14, 0.0};
    ax->particle_params = params;
    ```

## Einstein, Infeld, Hoffman GR treatment

The Einstein, Infeld, Hoffman GR treatment (EIH) is very accurate up to order \((v/c)^2\), but the implementation is slow and computationally expensive.
By default only the Sun is considered as a source for general relativistic correction.
However, for some cases such as when an asteroid has a close encounter with a planet, it can be important to include general relativistic corrections for other 
planets as well.
The following code turns on EIH GR corrections for the Sun and all planets.

!!! Info inline end
    The numerical value is 11 because there are 9 planets (including Pluto), the Sun, and the Moon.

=== "Python"

    ```python
    extras = assist.Extras(sim, ephem)
    extras.gr_eih_sources = 11
    ```
    
=== "C"

    ```c
    struct assist_extras* ax = ...
    ax->gr_eih_sources = 11
    ```

## Available force routines

The following table provides a summary of all available force routines in ASSIST.

Name (C)   | Name (python) | Numerical value (hex) | Default | Description
---------- | ------------- | --------------------- | ------- | -----------
ASSIST_FORCE_SUN               | SUN               | 0x01  | On  | Gravitational interaction (\(1/r\)) due to the Sun 
ASSIST_FORCE_PLANETS           | PLANETS           | 0x02  | On  | Gravitational interaction (\(1/r\)) due to all planets
ASSIST_FORCE_ASTEROIDS         | ASTEROIDS         | 0x04  | On  | Gravitational interaction (\(1/r\)) due to 16 massive asteroids
ASSIST_FORCE_NON_GRAVITATIONAL | NON_GRAVITATIONAL | 0x08  | On  | Non-gravitational effects from radiation pressure, outgassing, etc. Using the Marsden (1973) model. 
ASSIST_FORCE_EARTH_HARMONICS   | EARTH_HARMONICS   | 0x10  | On  | Earth's \(J_2\), \(J_3\), \(J_4\), \(J_5\) zonal harmonics  
ASSIST_FORCE_SUN_HARMONICS     | SUN_HARMONICS     | 0x20  | On  | \(J_2\) zonal harmonic of the Sun
ASSIST_FORCE_GR_EIH            | GR_EIH            | 0x40  | On  | Einstein, Infeld, Hoffman GR treatment (accurate but expensive)
ASSIST_FORCE_GR_SIMPLE         | GR_SIMPLE         | 0x80  | Off | Nobili and Roxburgh GR treatment 
ASSIST_FORCE_GR_POTENTIAL      | GR_POTENTIAL      | 0x100 | Off | Damour and Deruelle GR treatment
