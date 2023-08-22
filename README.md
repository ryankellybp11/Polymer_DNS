# dns-code
In-house direct numerical simulation code for channel/boundary layer flows with an external force field, and optional passive scalar and polymer additives.

The solver is based on [*Kim, J., Moin, P., & Moser, R. (1987). Turbulence statistics in fully developed channel flow at low Reynolds number. Journal of Fluid Mechanics, 177, 133-166.*](https://doi.org/10.1017/S0022112087000892)

The external force field is based on [*Goldstein, D., Handler, R., Sirovich, L. (1993). Modeling a No-Slip Flow Boundary with an External Force Field. Journal of Computational Physics, Volume 105, Issue 2, Pages 354-366.*](http://cfpl.ae.utexas.edu/wp-content/uploads/1993/04/goldstein_handler_sirovich_1993.pdf)

This branch has a self-contained case: all the necessary code, binaries, and files are located within the same directory. This should make setting up the code easier.

## Code structure:


    .
    ├── code                             # Contains source code, compiled code, binariy files, and geometry file
    │   ├── bin                          # Contains all binary files
    │   │   ├── dns                      # Main executable - runs DNS
    │   │   ├── geometry                 # Contains geometry generator and geometry file
    │   │   │   ├── Geom                 # Compiled geometry generator
    │   │   │   └── grid_size.mod        # Binary module for grid size data
    │   │   ├── grid_size.mod            # Binary module for grid size data
    │   └── src                          # Contains all the source files
    │       ├── dns                      # Contains DNS source files
    │       │   ├── dns.f90              # Main DNS source file (with polymer)
    │       │   ├── dns_no_poly.f90      # Main DNS source file (without polymer)
    │       │   ├── ffts.f90             # Useful fft functions for dns.f90
    │       │   ├── grid_size.f90        # Module with all grid size data
    │       │   ├── init_flow.f90        # Defines flow type (Channel, Couette, BL, etc.)
    │       │   ├── part_track.f90       # Particle tracking integrator
    │       │   └── polymer_routines.f90 # Additional subroutines for polymer effects
    │       ├── geometry                 # Contains geometry generator source code
    │       │   └── Geom.f90             # Geometry generator source code
    │       └── preconfigure.py          # preconfigure tool: compiles all the necessary files. Run with Python >=2.7.13 
    │── setup                            # Contains all the setup files that are inputs to the DNS. 
    │   ├── dns.config                   # Main DNS input file
    │   ├── vort.config                  # Contains input information for generating vortices (soon to be vortex ring as well, but as of now, that's in dns.config)
    │   ├── restart                      # Restart file (outputs/last-restart file of previous simulation)
    │   ├── c-restart                    # Scalar/polymer restart file 
    │   └── particles			         # Particle setup files
    │       ├── generator.f90            # Generates random particles in flow domain 
    │       ├── particles.dat            # Contains particle trajectory info (ASCII)
    │       └── pgen					 # Binary particle generator file   
    ├── outputs                          # Contails all the outputs produced by the DNS
    │   ├── flowfield                    # Contains tecplot-readable binary data files of the flow field
    │   │   ├── grid.plt                 # Grid file
    │   │   ├── time-00000.plt           # Flow field at time step 0 (in Tecplot binary format)
    │   │   └── ...                      # etc.
    │   ├── particles                    # Contains particle data
    │   │   ├── part-t-00000.dat         # Particles at time step 0 (in Tecplot ASCII format)
    │   │   └── ...                      # etc.
    │   ├── last-restart                 # Latest restart file. Do `mv outputs/last-restart setup/restart` if you want to restart latest simulation.
    │   ├── c-last-restart               # Latest polymer restart file. `mv outputs/c-last-restart setup/c-restart` if you want to restart latest simulation.
    │   ├── err                          # Runtime errors will appear here
    │   ├── log                          # Runtime outputs
    │   └── KE,enstrophy,...             # Channel flow output files for comparing different flow variables
    ├── jobcode                          # File for submitting a job to TACC. Usage: `sbatch jobcode`
    ├── clean                            # File for cleaning outputs directory. Usage: `. clean'
    ├── run                              # File for running the code locally (or on a work node). Usage: `. run'
    ├── hushrun                          # Same as `run' but runs the code in the background
    └── err_comp.txt                     # Output file for any compilation errors (empty if no errors)



## Prerequisites (OPTIONAL)
The code can export the data in three different formats:
1. ASCII (but tecplot-readable)
2. ASCII and then converted to Tecplot binary (.plt) using preplot
3. SZL Tecplot binary (.szplt) directly.

The third option is MUCH faster, but it requires installing Tecplot and the TecIO libraries as follows:

For fast data export, the TecIO library from Tecplot 360 is used. This library comes preinstalled with tecplot. 
### To install Tecplot from scratch (if not already installed):
1. Download the latest tecplot version:
```console
wget https://tecplot.azureedge.net/products/360/2021r1/tecplot360ex2021r1_linux64.sh
```
2. Install tecplot:
```console
chmod +x tecplot360ex2021r1_linux64.sh
./tecplot360ex2021r1_linux64.sh
```
On prompt, choose a *local* installation directory (e.g. `~/applications/tecplot`)

### Add relevant library paths to .bashrc:
1. Find where tecplot was installed. If tecplot is already installed, then:

```console
which tec360
```

will give you something like `path_to_tecplot/bin/tec360`. We need only the `path_to_tecplot`.

2. Replace `path_to_tecplot` with your installation directory and add the following lines at the bottom of `~/.bashrc`:

```console
export TECPLOT_360_INCLUDE=path_to_tecplot/include/
export TECPLOT_360_LIB=path_to_tecplot/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TECPLOT_360_LIB
```

3. Source `~/.bashrc`:
```console
source ~/.bashrc
```

## Installation:
1. Add the following line to `~/.bashrc`:
```console
alias preconfigure-dns='python code/src/preconfigure.py'
```
2. Source `~/.bashrc`:
```console
source ~/.bashrc
```

## Compile:
Run `preconfigure-dns` from the case directory after each change in the grid size or geometry.

## Run:
There are a few ways to run the code:
- Run the DNS code in a single core, locally:
```console
./code/bin/dns
```
- Run the DNS code with multipe cores, locally:
```console
./omprun
```
- Run the DNS code as a job in TACC:
```console
sbatch jobcode
```
where `jobcode` has been **edited** (job name, number of nores/mpi tasks, time, email, etc) according to the user's needs.

## Using Github in TACC:
TACC does not allow one to create ssh keys within a cloud node. In order to clone/pull/push this repo on a TACC server, one has to use Personal Access Tokens: https://docs.github.com/en/github/authenticating-to-github/keeping-your-account-and-data-secure/creating-a-personal-access-token

In order to save a personal access token on TACC (to avoid entering the token every time), do the following:
```console
git config credential.helper store
git clone [LINK TO THIS REPO]
```
or
```console
git config credential.helper store
git pull
```
from within this repo, and enter your username and **access token instead of your password**.

## Known issues with current version of code:
- Polymer targeting unexpectedly causes CFL failure
- Blasius BL flow initialization is extremely slow
- Misc. unused variables, inappropriately allocated arrays
