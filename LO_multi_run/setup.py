'''
This script is designed to set up multiple Lamb-Oseen vortex runs by copying the source code into several directories within this folder and altering desired parameters (by default, this will be vR)
'''
import numpy as np
import os

### CREATE COMPARISON CASES ###
num_runs = 8
for n in range(num_runs):
    os.system(f'cp -r ../LambOseen case{n}')

# Define grid/domain parameters for comparison runs
nsteps = 30000
iprnfrq = 1000
dt = 0.0005
Nx, Ny, Nz = 2, 256, 256
Lx, Ly, Lz = 5/64, 4.0, 4.0

# Define constant Polymer/Scalar Parameters
tpoly = 0.1
Sc_p = 10 # DIFFPOLY in DNS
Sc   = 10 # DIFF in DNS
deltaT = 10 
sigmap = np.pi/10

# Define constant vortex parameters for comparison runs
Gamma, sigma = 4.0, 0.1
vY, vZ = 0.0, 0.0

# Set variable to compare
vRlist = [sigma*0.25*(i+1) for i in range(num_runs)]

### BUILD SETUP FILES ###
# Read dns.config and vort.config for writing later
dnsin = open('case0/setup/dns.config','r') # Choose the first one since they're all identical
dnslines = dnsin.readlines()
dnsin.close()

vortin = open('case0/setup/vort.config','r')
vortlines = vortin.readlines()
vortin.close()

# Update lines with specified parameters above for all cases
for n in range(num_runs):
    # dns.config lines
    linetxt = f'    {nsteps}        ! NSTEPS      ( number of steps to run before stopping )\n'
    dnslines[3] = linetxt
    linetxt = f'    {iprnfrq}       ! IPRNFRQ     ( number of steps between each output of the flowfield )\n'
    dnslines[4] = linetxt
    linetxt = f'    {dt}            ! DT          ( time step: 0.00125/0.000625/0.0003125/0.00015625 )\n'
    dnslines[5] = linetxt
    linetxt = f'    {Ny}            ! NY1         ( number of grid cells in the y-direction - defined in grid_size module )\n'
    dnslines[9] = linetxt
    linetxt = f'    {Nz}            ! NZ1         ( number of grid cells in the z-direction on the regular grid - defined in grid_size module )\n'
    dnslines[10] = linetxt
    linetxt = f'    {Nx}            ! NX1         ( number of grid cells in the x-direction on the regular grid - defined in grid_size module )\n'
    dnslines[11] = linetxt
    linetxt = f'    {Ly}            ! YL          ( length of the computational domain in the wall-normal direction )\n'
    dnslines[13] = linetxt
    linetxt = f'    {Lz}            ! ZL          ( length of the computational domain in the spanwise direction ) 3.1415927\n'
    dnslines[14] = linetxt
    linetxt = f'    {Lx}            ! XL          ( length of the computational domain in the streamwise direction ) 12.566371\n'
    dnslines[15] = linetxt
    linetxt = f'    {tpoly}                          ! tpoly - polymer relaxation time\n'
    dnslines[55] = linetxt
    linetxt = f'    {Sc_p}                           ! DIFFPOLY - Schmidt number of polymer (numerical stabilizer)\n'
    dnslines[57] = linetxt
    linetxt = f'    {Sc}                             ! DIFF - Prandtl/Schmidt number of scalar\n'
    dnslines[62] = linetxt
    linetxt = f'    {deltaT}                         ! deltaT - Max concentration (PPM) or scalar emission rate (PPM/(time step)/particle)\n'
    dnslines[63] = linetxt
    linetxt = f'    100.0, {sigmap}, {sigmap}          ! sigmax,sigmay,sigmaz\n'
    dnslines[64] = linetxt

    # vort.config lines
    linetxt = f'    {Gamma}, {sigma}   ! vortGamma, vortSigma\n'
    vortlines[13] = linetxt
    linetxt = f'    {vY}, {vZ}      ! vortY, vortZ      \n'
    vortlines[14] = linetxt
    linetxt = f'    {vRlist[n]}     ! vR - Radial distance of scalar from vortex center\n'
    vortlines[17] = linetxt

    # Write output files
    dnsout = open(f'case{n}/setup/dns.config','w+')
    dnsout.writelines(dnslines)
    dnsout.close()

    vortout = open(f'case{n}/setup/vort.config','w+')
    vortout.writelines(vortlines)
    vortout.close()

# Create commandlines file
here = os.getcwd()
cmdfile = open('commandlines','w')
cmdlines = []
for n in range(num_runs):
    linetxt = f'cd {here}/case{n} && make polymer_dns && . run \n'
    cmdlines.append(linetxt)

cmdfile.writelines(cmdlines)
cmdfile.close()