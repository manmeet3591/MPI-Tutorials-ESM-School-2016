CC = mpicc
FC = mpiifort
all:
	$(FC) dotproduct.f90 -o dotproduct.x

extra:
	$(CC) extra_dotproduct.c -o extra_dotproduct.x

serial:	
#	gcc serial_dotproduct.c -o serial_dotproduct.x
#	gcc -lm poisson_serial.c -o poisson_serial.x
#	ifort serial_dotproduct.f90 -o serial_dotproduct.x
	ifort poisson_serial_jacobi.f90 -o poisson_serial_jacobi.x
