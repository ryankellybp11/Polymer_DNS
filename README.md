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
    │   │   ├── ***.mod                  # Binary module files for use in main DNS code
    │   └── src                          # Contains all the source files
    │       ├── dns                      # Contains DNS source files
    │       │   ├── derivatives.f90      # Contains subroutines related to gradients and derivatives
    │       │   ├── grid_size.f90        # Module with all grid size data
    │       │   ├── helpers.f90          # Contains various helpful functions and subroutines
    │       │   ├── init_flow.f90        # Defines flow type (Channel, Couette, BL, etc.)
    │       │   ├── main_dns.f90         # Main DNS source file 
    │       │   ├── part_track.f90       # Particle tracking integrator
    │       │   └── solvers.f90          # Contains the solvers used in the code, including vcw3dp
    │       ├── geometry                 # Contains geometry generator source code
    │       │   └── Geom.f90             # Geometry generator source code
    │       └── preconfigure.py          # preconfigure tool: compiles grid and geometry files. Run with Python >=2.7.13 
    ├── outputs                          # Contails all the outputs produced by the DNS
    │   ├── flowfield                    # Contains tecplot-readable binary data files of the flow field
    │   │   ├── grid.plt                 # Grid file
    │   │   ├── time-000001.szplt           # Flow field at time step 1 (in Tecplot binary format)
    │   │   └── ...                      # etc.
    │   ├── particles                    # Contains particle data
    │   │   ├── part-t-000001.dat         # Particles at time step 1 (in Tecplot ASCII format)
    │   │   └── ...                      # etc.
    │   ├── last-restart                 # Latest restart file. Copy to setup/restart to restart latest simulation
    │   ├── c-last-restart               # Latest polymer restart file. Copy to setup/c-restart to restart latest simulation
    │   ├── err                          # Runtime errors will appear here
    │   ├── log                          # Runtime outputs
    │   └── KE,enstrophy,...             # Output files for comparing different flow variables
    │── setup                            # Contains all the setup files that are inputs to the DNS. 
    │   ├── dns.config                   # Main DNS input file
    │   ├── vort.config                  # Contains input information for generating vortices 
    │   ├── restart                      # Restart file (outputs/last-restart file of previous simulation)
    │   ├── c-restart                    # Scalar/polymer restart file 
    │   └── particles			         # Particle setup files
    │       ├── generator.f90            # Generates random particles in flow domain 
    │       ├── particles.dat            # Contains particle trajectory info (ASCII)
    │       └── pgen					 # Binary particle generator file   
    ├── bash_clean                       # File for cleaning outputs directory. Usage: `. bash_clean'
    ├── err_comp.txt                     # Output file for any compilation errors (empty if no errors)
    ├── hushrun                          # Same as `run' but runs the code in the background
    ├── jobcode                          # File for submitting a job to TACC. Usage: `sbatch jobcode`
    ├── Makefile                         # Makefile for compiling DNS code - uses preconfigure.py
    └── run                              # File for running the code locally (or on a work node). Usage: `. run'


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
- Run the DNS code locally or on a compute node on TACC:
```console
. run
```

- Run the DNS code in the background locally or on a compute node on TACC:
```console
. hushrun
```

- Run the DNS code as a batch job on TACC:
```console
sbatch jobcode
```
where `jobcode` has been **edited** (job name, number of nores/mpi tasks, time, email, etc) according to the user's needs.

You can also clean up the output files by running the script
```console
. clean
```

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
