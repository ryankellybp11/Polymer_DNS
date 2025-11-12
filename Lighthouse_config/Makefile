#COMP := ifx
COMP := ifort
LIB := -L$(FFTW_LIB) -lfftw3 -lm -L$(TEC360_LIB) -lirc -ltecio -L$(HDF5_LIB) -lhdf5_fortran -lhdf5
INC := -I$(FFTW_INCLUDE) -I$(TEC360_INC) -DOUTPUTFORM=3 -I$(HDF5_INCLUDE) 
FLAGS := -r8 -fr -fpp -132 -O1 -qmkl -parallel -march=native -mtune=native -axicelake-client -unroll -shared-intel -mcmodel=medium
#FLAGS := -r8 -fr -fpp -132 -O0 -mcmodel=medium -shared-intel -xhost -qmkl -parallel -g -traceback -debug full -warn all -check all -check noarg_temp_created -nogen-interfaces
OMPFLAGS := -qopenmp -lfftw3_omp -qopenmp-simd
DIR := code/src/dns/
POLY := -DSCALAR -DPOLYMER
SCL := -DSCALAR


temp: temp.f90 Makefile
	ifort temp.f90 $(FLAGS) -o temp
	clear

dns:
	python3 code/src/preconfigure.py
	$(COMP) $(DIR)grid_size.f90 $(DIR)derivatives.f90 $(DIR)helpers.f90 $(DIR)solvers.f90 $(DIR)main_dns.f90 $(DIR)init_flow.f90 $(DIR)part_track.f90 -o code/bin/dns $(LIB) $(INC) $(FLAGS) $(OMPFLAGS) > outputs/log 2> err_comp.txt 
#	mv *genmod* code/bin/.
	mv *.mod code/bin/.

scalar_dns:
	python3 code/src/preconfigure.py
	$(COMP) $(DIR)grid_size.f90 $(DIR)derivatives.f90 $(DIR)helpers.f90 $(DIR)solvers.f90 $(DIR)main_dns.f90 $(DIR)init_flow.f90 $(DIR)part_track.f90 -o code/bin/dns $(LIB) $(INC) $(FLAGS) $(OMPFLAGS) $(SCL) > outputs/log 2> err_comp.txt 
#	mv *genmod* code/bin/.
	mv *.mod code/bin/.

polymer_dns:
	python3 code/src/preconfigure.py
	$(COMP) $(DIR)grid_size.f90 $(DIR)derivatives.f90 $(DIR)helpers.f90 $(DIR)solvers.f90 $(DIR)main_dns.f90 $(DIR)init_flow.f90 $(DIR)part_track.f90 -o code/bin/dns $(LIB) $(INC) $(FLAGS) $(OMPFLAGS) $(POLY) > outputs/log 2> err_comp.txt 
#	mv *genmod* code/bin/.
	mv *.mod code/bin/.

clean: 
	rm -rf temp fort.* a.out *.mod
	rm -rf *_genmod* *.opt* r00*
	clear
