FORTRAN = gfortran

all: atomic serial tpriv

serial: serial_hw4.f95
	$(FORTRAN) -o serial -fopenmp serial_hw4.f95
atomic: atomic_hw4.f95
	$(FORTRAN) -o atomic -fopenmp atomic_hw4.f95
tpriv: tpriv_hw4.f95
	$(FORTRAN) -o tpriv -fopenmp tpriv_hw4.f95
clean:
	rm -f *.a *.o a.out core* atomic serial tpriv
	
