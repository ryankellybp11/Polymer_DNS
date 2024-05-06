LIB := -L$(TACC_FFTW3_LIB) -lfftw3 -lm
INC := -I$(TACC_FFTW3_INC)
FLAGS := -r8 -mavx -132

test: fftw_test.f90 Makefile
	ifort code/src/dns/grid_size.f90 fftw_test.f90 -o test $(LIB) $(INC) $(FLAGS)
	clear

clean: 
	rm -f test fort.*
	clear
