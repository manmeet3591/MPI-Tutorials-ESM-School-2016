CC = mpicc
FC = mpiifort
all:
	$(FC) trapezoid.f90 -o trapezoid.x

serial:
#	gcc -lm serial_trapezoid.c -o serial_trapezoid.x
	ifort serial_trapezoid.f90 -o serial_trapezoid.x
