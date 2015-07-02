Solia is a collection of utilities for exploring the space of possible planetary systems. More specifically, it was inspired by and focused on a need to investigate the stability of co-orbital systems as exemplified by Saturn's moons Janus & Epimetheus. It consists of numerical integrators for n-body gravitational simulations and a set of graphing and data analysis tools.

Data Format
=====
The different components of Solia all communicate through a common text-based data format for simulation data. Data files consist of a series of snapshots which have the following format:

    N T
    m_1 x_1 y_1 vx_1 vy_1
    ...
    m_n x_n y_n vx_n vy_n
    r_1,2 ... r_(n-1),n

where N is the total number of particles in the simulation, T is the current time, m_1 through m_n are the particles masses, x_i and y_i are the coordinates for particle i, vx_i and vy_i are the velocity components for particle i, and r_i,j is the distance between particles i and j. All values are in MKS units. The last line could be recalculated from other information rather stored in the data file, but is included for the sake of efficiency- simulators have to calculate it anyway, so we might as well make things easier on the analysis programs.

A future version may extend this to arbitrary numbers of dimensions. Altering the code to handle 3D simulations is pretty simple, but currently has to be done by hand.

Simulators
=====

Solia includes two numerical n-body simulators:

* NBody.py is a very simple second-order simulator, useful for playing around but not terribly fast nor terribly accurate.
* nbody.cpp is a much faster and more accurate 4th-order simulator with variable global timestep based on the Hermite integrator starter code by [Piet, Makino, and McMillan](https://www.ids.ias.edu/~piet/act/comp/algorithms/codes.html). It is built with the included Makefile.

While there is plenty of other n-body simulation software out there, I wrote these because I could not find any that fit all three of the following criteria:

1. Takes input and produces output in standard units. E.g., the original code on which nbody.cpp was based is unitless (implicitly using natural units), which introduces unnecessary confusion in how to translate the representation of a given system into a form that the integrator can use, and then translate the results back into something easy for a person to understand.
2. Has a sufficiently high degree efficiency and accuracy over several years of simulation time.
3. Makes it easy to input data for arbitrary theoretical systems. Many high-accuracy simulators that produce output in human-readable units are, for example, hard-coded for solar system simulation, and make it difficult or impossible to do arbitrary experiments.

Both simulators read input from stdin and write to stdout. The first snapshot in an input data file is used to provide initial conditions to start a simulation; any additional snapshots are ignored. nbody.cpp also writes diagnostic information (primarily energy conservation info) to stderr.

NBody.py takes three optional command-line arguments:

    --step [seconds]: the timestep size
    --end [seconds]: the ending time of the simulation
    --out [seconds]: the interval at which to produce output snapshots

nbody.cpp takes five optional command-line arguments:

    -a [float]: accuracy parameter, used to control the size of the variable timestep
    -d [seconds]: diagnostic interval, the simulation time between diagnostic output
    -x: eXtra diagnostics; setting this flag causes diagnostics to output all internal particle data
    -o [seconds]: output interval, the simulation time between output snapshots
    -t [seconds]: total simulation duration; this is slightly different from "--end" for NBody.py, in that it specifies duration after the initial start time read from the input file, which may be greater than zero, not a hard end time.

Note that, due to the variable timestep, output times and total duration may not match the provided parameters exactly, but output will occur as close as soon as possible after each scheduled interval.

Graphing Tools
=====

* GraphOrbits.py produces line graphs showing the varying distances over time between every pair of particles. Command-line options allow for averaging over a rolling time window (to smooth out regular oscillations due to orbital eccentricity, for example), selecting a specific time range out of a larger data set, altering the time-axis scale (so you can view years instead of megaseconds, for example), and saving graphs to a file vs. displaying in a window.

* GraphMotion.py produces animations of the evolution of a system in 2D space over time. It is particularly useful for debugging errors in initial conditions (like a planet moving in the wrong direction), which can be difficult to see on the distance graphs.

Other
=====

* GenerateSystems.py produces an ensemble of systems with different initial conditions and produces a DOS Batch file to run simulations on each one and save graphs of the results.
* OrbitData.py is a Python module for reading and writing the common data format which is used by the other Python scripts.