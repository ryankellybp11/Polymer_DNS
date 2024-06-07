#COMP := ifx
COMP := ifort
LIB := -L$(TACC_FFTW3_LIB) -lfftw3 -lfftw3f -lm -L$(TECPLOT_360_LIB) -lirc -ltecio
INC := -I$(TACC_FFTW3_INC) -I$(TECPLOT_360_INCLUDE) -DOUTPUTFORM=3
#FLAGS := -r8 -fr -fpp -132 -O1 -qmkl -parallel 
FLAGS := -r8 -fr -fpp -132 -O3 -qmkl -parallel -march=native -mtune=native -axicelake-client -unroll -shared-intel -g -fp-model fast=2
#FLAGS := -r8 -fr -fpp -mavx -132 -O0 -mcmodel=medium -shared-intel -xhost -qmkl -parallel -g -traceback -debug full -warn all -check all -check noarg_temp_created -nogen-interfaces
OMPFLAGS := -qopenmp -lfftw3_omp -qopenmp-simd
#DIR := main_test_code/src/dns/
DIR := code/src/dns/
POLY := -DSCALAR -DPOLYMER
SCL := -DSCALAR


old_dns: code/src/dns/* Makefile setup/dns.config
	python3 old_source_code/src/preconfigure.py
	
fftw_test: fftw_test.f90 Makefile
	ifort code/src/dns/grid_size.f90 fftw_test.f90 -o fftw_test $(LIB) $(INC) $(FLAGS) 
	clear

omp_compare_fft: omp_compare_ffts.f90 Makefile
	mpiifort code/src/dns/grid_size.f90 code/src/dns/ffts.f90 omp_compare_ffts.f90 -o omp_compare_ffts $(LIB) $(INC) $(FLAGS) $(OMPFLAGS)
	clear

compare_fft: compare_ffts.f90 Makefile
	mpiifx code/src/dns/grid_size.f90 code/src/dns/ffts.f90 compare_ffts.f90 -o compare_ffts $(LIB) $(INC) $(FLAGS)
	clear

yfft_test: yfft_test.f90 Makefile
	ifx code/src/dns/grid_size.f90 yfft_test.f90 -o yfft_test $(LIB) $(INC) $(FLAGS)
	clear

temp: temp.f90 Makefile
	ifort temp.f90 $(FLAGS) -o temp
	clear

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
dns:
	python3 code/src/preconfigure.py
	$(COMP) $(DIR)grid_size.f90 $(DIR)derivatives.f90 $(DIR)helpers.f90 $(DIR)solvers.f90 $(DIR)main_dns.f90 $(DIR)init_flow.f90 $(DIR)part_track.f90 -o code/bin/dns $(LIB) $(INC) $(FLAGS) $(OMPFLAGS) > outputs/log 2> err_comp.txt 
#	mv *genmod* code/bin/.
	mv *.mod code/bin/.
	

clean: 
	rm -rf fftw_test temp fort.* a.out compare_ffts yfft_test omp_compare_ffts
	rm -rf *_genmod* *.opt* r00*
	clear
