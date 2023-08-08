# dns-code
In-house direct numerical simulation code for channel/boundary layer flows with an external force field.

The solver is based on [*Kim, J., Moin, P., & Moser, R. (1987). Turbulence statistics in fully developed channel flow at low Reynolds number. Journal of Fluid Mechanics, 177, 133-166.*](https://doi.org/10.1017/S0022112087000892)

The external force field is based on [*Goldstein, D., Handler, R., Sirovich, L. (1993). Modeling a No-Slip Flow Boundary with an External Force Field. Journal of Computational Physics, Volume 105, Issue 2, Pages 354-366.*](http://cfpl.ae.utexas.edu/wp-content/uploads/1993/04/goldstein_handler_sirovich_1993.pdf)

This branch has a self-contained case: all the necessary code, binaries, and files are located within the same directory. This should make setting up the code easier.

## Code structure:

    .
    ├── code                          # Contains source code, compiled code, binariy files, and geometry file
    │   ├── bin                       # Contains all binary files
    │   │   ├── dns                   # Main executable - runs DNS
    │   │   ├── geometry              # Contains geometry generator and geometry file
    │   │   │   ├── Geom              # Compiled geometry generator
    │   │   │   └── geometry          # Geometry file
    │   │   ├── grid_size.mod         # Binary module for grid size data
    │   │   ├── jet_data.mod          # Binary module for jet data
    │   │   ├── snap_data.mod         # Binary module for snapshot data
    │   │   └── tecplot               # Contains preplot tool
    │   │       ├── preplot           # Preplot tool: converts ascii data to tecplot-readable binaries
    │   │       ├── preplot-all       # Obsolete
    │   │       └── preplot-bin       # Called from preplot
    │   └── src                       # Contains all the source files
    │       ├── dns                   # Contains DNS source files
    │       │   ├── dns.f90           # Main DNS source file
    │       │   ├── ffts.f90          # Useful fft functions for dns.f90
    │       │   └── grid_size.f90     # Module with all grid size data
    │       ├── geometry              # Contains geometry generator source code
    │       │   └── Geom.f90          # Geometry generator source code
    │       └── preconfigure.py       # preconfigure tool: compiles all the necessary files. Run with Python >=2.7.13 
    ├── jobcode                       # File for submitting a job to TACC. Usage: `sbatch jobcode`
    ├── omprun                        # File for running the code locally, but in parallel (OMP)
    ├── outputs                       # Contails all the outputs produced by the DNS
    │   ├── err                       # Runtime errors will appear here
    │   ├── flowfield                 # Contains tecplot-readable binary data files of the flow field
    │   │   ├── grid.plt              # Grid file
    │   │   ├── time-00000.plt        # Flow field at time step 0 (in Tecplot binary format)
    │   │   └── ...                   # etc.
    │   ├── last-restart              # Latest restart file. Do `mv outputs/last-restart setup/restart` if you want to restart latest simulation.
    │   ├── log                       # Runtime outputs
    │   ├── particles                 # Contains particle data
    │   │   ├── particles-t-00000.dat # Particles at time step 0 (in Tecplot ASCII format)
    │   │   └── ...                   # etc.
    │   ├── snapshots                 # Contains snapshot data
    │   │   ├── gird.dat              # Grid of snapshot data
    │   │   ├── input.dat             # Control inputs at each snapshot time step
    │   │   ├── u.dat                 # u-velocity shapshots
    │   │   ├── v.dat                 # v-velocity shapshots
    │   │   ├── w.dat                 # w-velocity shapshots
    │   │   └── wy.dat                # wy-vorticity shapshots
    │   └── statistics                # Contains statistics data (currently not used)
    ├── README.md                     # Readme file
    ├── run                           # File 
    └── setup                         # Contains all the setup files that are inputs to the DNS. 
        ├── dns.config                # Main DNS input file
        ├── jet.config                # Jet setup file
        ├── jet_input_signal.dat      # Jet input signal controlling the magnitude of the jet
        ├── restart                   # Restart file (outputs/last-restart file of previous simulation)
        ├── snapshots.config          # Snapshots setup file
        └── tracer_particles.config   # Tracer particles configuration file

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

