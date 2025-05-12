
## N-Body Gravitational Simulation

This project simulates the gravitational interactions between multiple celestial bodies using Newtonian mechanics and the Velocity Verlet integration algorithm. The core numerical computations are written in C++, with a Python wrapper for interfacing and visualizing results.

Key outputs of the simulation model include: 
- 2D animations of the orbital motions 
- Calculation of the orbital eccentricities and periods
- Analysis of energy conservation and integration error 
- Parametric studies (e.g., Moon mass variation and solar mass effects)
The code is currently set up for a 3-body Earth–Moon–Sun system.

### Project structure

The file `orbit_simul.cpp` includes the C++ implementation of the Velocity Verlet algorithm and the energy and eccentricity calculations. 

The Python script `orbit_simul_wrapper.py` is the main interface which: 
- Loads and runs the C++ shared library using `ctypes`
- Sets up initial conditions for different systems
- Calls the C++ functions to perform simulations
- Analyses and plots simulation results by accessing the plotting functions that are stored in `orbit_plot_module.py`.


### How to run

To compile the C++ code into a shared library (so the python wrapper can use it):

        g++ -fPIC -shared orbit_simul.cpp -o orbit_simul.so

To run the python wrapper:

        python3 orbit_simul_wrapper.py


### Example 3-Body Earth, Sun and Moon system

<img src="https://github.com/nicolaegues/orbit-cpp-project/blob/main/my_orbits.gif" width="50%" height="50%">
