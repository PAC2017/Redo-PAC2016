
FC=ifort
LD= ifort
FFLAGS=-qopenmp -mcmodel=large

main: main.o Covariance1.o Covariance.o check.o
	${FC} ${FFLAGS} main.o Covariance1.o Covariance.o check.o -o main
main.o: main.f90
	${FC} ${FFLAGS} -c main.f90
Covariance1.o: Covariance1.f90 
	${FC} ${FFLAGS} -c Covariance1.f90
#Covariance.o: Covariance.f90
#	${FC} ${FFLAGS} -c Covariance.f90
Covariance.o: Covariance.c
	icc -fopenmp -std=c99 -axCOMMON-AVX512 -O3 -ipo -xMIC-AVX512 -funroll-loops -fast -march=native -mtune=core-avx2 -xHost -no-prec-div -fomit-frame-pointer -qopt-streaming-stores=always -c Covariance.c
check.o: check.f90
	${FC} ${FFLAGS} -c check.f90

clean:
	rm -rf *.o *.i main
