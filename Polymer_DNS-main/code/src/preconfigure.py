# Preconfigure dns code before running

from os.path import expanduser
import subprocess
import os

home = expanduser('~')
compiler = ['ifort']


print('\n')
print('____________________________________________\n')
print('         Preconfiguring DNS code')
print('    ____________________________________    \n')


# 1. Read grid size and generate grid_size.f90 file (module grid_size). Save grid_size.f90 to setup/bin/
print('    1. Setting up the grid_size module... ')

# Read file with grid data
setup_dns = open('setup/dns.config', 'r')
setup_data = list([line.strip().split() for line in setup_dns])
setup_dns.close()

# Remember that Python indexing starts at 0
ny = int(setup_data[9][0])
nz = int(setup_data[10][0])
nx = int(setup_data[11][0])
bftail = int(setup_data[12][0])

output_format = int(setup_data[33][0]) # For writing Tecplot files

ipoly = int(setup_data[40][0])
iscl = int(setup_data[41][0])
npart = int(setup_data[47][0])
new_part = int(setup_data[48][0])


grid_size = open('code/src/dns/grid_size.f90', 'r')
lines = grid_size.readlines()
lines[1] = '    integer, parameter :: nx_ = %s, ny_ = %s, nz_ = %s, bftail_ = %s, npart = %s \n' % (nx, ny, nz, bftail, npart)
grid_size.close()

grid_size_updated = open('code/src/dns/grid_size.f90', 'w+')
grid_size_updated.writelines(lines)
grid_size_updated.close()


print('       Done!\n')


# 2. Make geometry

print('    2. Generating geometry file... ')

os.chdir('code/bin/geometry/')
comp_out, comp_err = subprocess.Popen(compiler + ['-g', '-traceback','../../src/dns/grid_size.f90', '../../src/geometry/Geom.f90', '-o', 'Geom'], stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()
# os.chmod 755 'Geom'
exec_err = subprocess.call('./Geom', stderr=subprocess.PIPE)

if not comp_err and not exec_err:
    print('       Done!\n')
elif comp_err:
    print('Generating geometry file failed!')
    print('Compilation error:')
    print(comp_err)
else:
    print('Generating geometry file failed!')
    print('Runtime error:')
    print(exec_err)

os.chdir('../../..')


# 3. Generate particles
if new_part > 0:
	if npart == 1:
		print('    3. Generating 1 particle... ')
	elif new_part == 1:	
		print('    3. Generating',npart,'random particles... ')
	elif new_part == 2:	
		print('    3. Generating',npart,'random particles in a sphere... ')
	elif new_part == 3:	
		print('    3. Generating',npart,'particles in a circle... ')
	elif new_part == 4:	
		print('    3. Generating',npart,'random particles in a sheet... ')
	elif new_part == 5:	
		print('    3. Generating',npart,'random particles near the walls... ')
	elif new_part == 6:	
		print('    3. Generating',npart,'random particles in a cylinder... ') 
	elif new_part == 7:	
		print('    3. Generating',npart,'structured particles in a y-z plane... ') 
	elif new_part == 10:	
		print('    3. Generating',npart,'particles to be released regularly... ') 
	elif new_part == 11:	
		print('    3. Generating',npart,'structured particles throughout the domain... ') 
	else:
		print('    3. Generating',npart,'particles... ') 
	
	os.chdir('setup/particles/')
	comp_out, comp_err = subprocess.Popen(compiler + ['generator.f90', '-o', 'pgen'], stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	exec_err = subprocess.call('./pgen', stderr=subprocess.PIPE)
	
	if not comp_err and not exec_err:
	    print('       Done!\n')
	elif comp_err:
	    print('Generating particle file failed!')
	    print('Compilation error:')
	    print(comp_err)
	else:
	    print('Generating particle file failed!')
	    print('Runtime error:')
	    print(exec_err)
	
	os.chdir('../..')
	
	
	# 4. Compile code from the source code in code_dir and save executable to setup/bin
	
	if (ipoly != 0) or (iscl != 0):
		print('    4. Compiling DNS code with polymer... ')
	else:
		print('    4. Compiling DNS code (no polymer)... ')
else:
	if (ipoly != 0) or (iscl != 0):
		print('    3. Compiling DNS code with polymer... ')
	else:
		print('    3. Compiling DNS code (no polymer)... ')


#compiler_options = ['-mavx', '-r8', '-fr', '-132', '-w', '-O2', '-unroll', '-qopenmp', '-parallel', '-mcmodel', 'medium', '-fpp', '-shared-intel','-g','-traceback','-MD']
# More debug options
compiler_options = ['-mavx', '-r8', '-fr', '-132', '-w', '-O3', '-unroll', '-qopenmp', '-parallel', '-mcmodel', 'medium', '-fpp', '-shared-intel','-g','-traceback','-Wall','-check all']

if (ipoly != 0) and (iscl != 0):
    # DNS w/ Polymer
    source_files = ['../src/dns/grid_size.f90', '../src/dns/dns.f90', '../src/dns/ffts.f90', '../src/dns/init_flow.f90', '../src/dns/part_track.f90', '../src/dns/polymer_routines.f90']
else:    
    # No poly dns
    source_files = ['../src/dns/grid_size.f90', '../src/dns/dns_no_poly.f90', '../src/dns/ffts.f90', '../src/dns/init_flow.f90', '../src/dns/part_track.f90']

other_options = ['-cref', '-wl', '-mkl', '-o']
exec_name = ['dns']

if output_format == 3:
	current_env = os.environ.copy()
	TECPLOT_360_INCLUDE = current_env["TECPLOT_360_INCLUDE"]
	TECPLOT_360_LIB = current_env["TECPLOT_360_LIB"]
	link_tecio = ['-lirc', '-I', TECPLOT_360_INCLUDE, '-L', TECPLOT_360_LIB, '-ltecio', '-DOUTPUTFORM=3']
	print('       (Using tecio)')
else:
	print('       (Not using tecio)')
	link_tecio = []


os.chdir('code/bin/')
#subprocess.call(['cp', code_dir + 'dns.f90', '.']) # Just a copy of the code - not used in compilation
comp_out, comp_err = subprocess.Popen(compiler + compiler_options + link_tecio + source_files + other_options + exec_name, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

err_str = comp_err.decode("utf-8")
with open('../../err_comp.txt','w') as f:
    f.write(err_str)


if not comp_err:
    print('       Done!\n')
elif 'warning' in err_str:
	print('       Warning:')
	print('       '+err_str)
	print('')
	print('       Code compilation is otherwise successful')
else:	
    print('compilation failed!\n')
    print('Check error in err_comp.txt\n')

os.chdir('../..')


if not comp_err:
    print('    ____________________________________    \n')
    print('         DNS code is ready to run!')
    print('____________________________________________\n')
elif 'warning' in err_str:
    print('    ____________________________________    \n')
    print('         DNS code is ready to run!')
    print('____________________________________________\n')
else:
    print('    ************************************    \n')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

